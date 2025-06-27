#include <iostream>
#include <vector>
#include <atomic> // For thread-safe counter
#include <thread> // For hardware_concurrency
#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStopwatch.h" // For timing
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

using namespace std;
using namespace ROOT;

// Function to compute momentum from px, py, pz
float computeMomentum(const RVec<float>& px, const RVec<float>& py, const RVec<float>& pz, size_t index) {
    return sqrt(px[index] * px[index] + py[index] * py[index] + pz[index] * pz[index]);
}

int main() {
    // Determine the number of threads supported by the CPU
    unsigned int nThreads = std::thread::hardware_concurrency();
    if (nThreads == 0) {
        cout << "Could not determine the number of threads; defaulting to 4." << endl;
        nThreads = 4; // Fallback value if hardware_concurrency fails
    }
    cout << "Number of threads supported by CPU: " << nThreads << endl;

    // Enable multi-threading in ROOT using the detected number of threads
    ROOT::EnableImplicitMT(nThreads);
    cout << "Multi-threading enabled for RDataFrame processing with " << nThreads << " threads." << endl;

    // Open the ROOT file and create an RDataFrame from the TTree
    TFile* inputFile = new TFile("output/particle_data.root", "READ");
    if (!inputFile->IsOpen()) {
        cerr << "Error: Could not open output/particle_data.root!" << endl;
        return 1;
    }

    // Create RDataFrame from the TTree
    ROOT::RDataFrame df("ParticleTree", inputFile);

    // Get the total number of events for progress tracking
    Long64_t totalEvents = *df.Count();
    cout << "Total number of events to process: " << totalEvents << endl;

    // Create the 2D histograms for pions (PID 211) and protons (PID 2212)
    TH2F hChi2PIDvsPPID211("hChi2PIDvsPPID211", "Chi2PID vs. Momentum for Positive Pions (PID 211);Momentum p (GeV);Chi2PID", 100, 0, 10, 100, -5, 5);
    TH2F hChi2PIDvsPPID2212("hChi2PIDvsPPID2212", "Chi2PID vs. Momentum for Protons (PID 2212);Momentum p (GeV);Chi2PID", 100, 0, 10, 100, -5, 5);

    // Thread-safe counter for processed events
    std::atomic<Long64_t> processedEvents(0);
    Long64_t progressInterval = totalEvents / 20; // Update progress every 5%
    if (progressInterval == 0) progressInterval = 1; // Avoid division by zero for small datasets

    // Set up a timer to measure the processing time
    TStopwatch timer;
    timer.Start();

    // Process the data in parallel using RDataFrame with Foreach
    df.Foreach([&](const RVec<int>& pid, const RVec<float>& chi2pid, const RVec<float>& px, const RVec<float>& py, const RVec<float>& pz) {
        // Increment the counter for each event (thread-safe)
        Long64_t currentEvent = processedEvents.fetch_add(1, std::memory_order_relaxed);

        // Print progress at regular intervals
        if (currentEvent % progressInterval == 0 || currentEvent == totalEvents - 1) {
            double percent = (currentEvent + 1) * 100.0 / totalEvents;
            cout << "\rProgress: " << fixed << setprecision(1) << percent << "%" << flush;
        }

        // Process the event for both pions and protons
        for (size_t i = 0; i < pid.size(); ++i) {
            if (chi2pid[i] == 9999.0) continue; // Skip invalid chi2pid

            float p = computeMomentum(px, py, pz, i);
            if (pid[i] == 211) { // Pions
                hChi2PIDvsPPID211.Fill(p, chi2pid[i]);
            } else if (pid[i] == 2212) { // Protons
                hChi2PIDvsPPID2212.Fill(p, chi2pid[i]);
            }
        }
    }, {"pid", "chi2pid", "px", "py", "pz"});

    // Ensure the progress line ends with a newline
    cout << endl;

    // Stop the timer and print the elapsed time in minutes
    timer.Stop();
    double timeInMinutes = timer.RealTime() / 60.0; // Convert seconds to minutes
    cout << "Time taken for RDataFrame processing: " << timeInMinutes << " minutes" << endl;

    // Debug: Check entries in the histograms
    cout << "Entries in hChi2PIDvsPPID211 (Pions): " << hChi2PIDvsPPID211.GetEntries() << endl;
    cout << "Entries in hChi2PIDvsPPID2212 (Protons): " << hChi2PIDvsPPID2212.GetEntries() << endl;

    // Plot: Chi2PID vs. Momentum for Positive Pions (PID 211) with logarithmic z-axis
    TCanvas* canvasChi2PIDvsPPion = new TCanvas("canvasChi2PIDvsPPion", "Chi2PID vs. Momentum for Positive Pions", 800, 600);
    if (hChi2PIDvsPPID211.GetEntries() > 0) {
        hChi2PIDvsPPID211.SetStats(0);
        hChi2PIDvsPPID211.SetMinimum(1e-5); // Set minimum for log scale
        hChi2PIDvsPPID211.Draw("COLZ");
        gPad->SetLogz(); // Set logarithmic z-axis
    } else {
        cout << "Warning: hChi2PIDvsPPID211 (Pions) is empty, plot will be empty." << endl;
    }
    canvasChi2PIDvsPPion->Print("output_mt/chi2pid_vs_p_pion.pdf");

    // Plot: Chi2PID vs. Momentum for Protons (PID 2212) with logarithmic z-axis
    TCanvas* canvasChi2PIDvsPProton = new TCanvas("canvasChi2PIDvsPProton", "Chi2PID vs. Momentum for Protons", 800, 600);
    if (hChi2PIDvsPPID2212.GetEntries() > 0) {
        hChi2PIDvsPPID2212.SetStats(0);
        hChi2PIDvsPPID2212.SetMinimum(1e-5); // Set minimum for log scale
        hChi2PIDvsPPID2212.Draw("COLZ");
        gPad->SetLogz(); // Set logarithmic z-axis
    } else {
        cout << "Warning: hChi2PIDvsPPID2212 (Protons) is empty, plot will be empty." << endl;
    }
    canvasChi2PIDvsPProton->Print("output_mt/chi2pid_vs_p_proton.pdf");

    // Clean up
    delete canvasChi2PIDvsPPion;
    delete canvasChi2PIDvsPProton;
    inputFile->Close();
    delete inputFile;

    // Disable multi-threading (good practice)
    ROOT::DisableImplicitMT();

    cout << "Program completed successfully." << endl;
    return 0;
}