#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

using namespace std;
using namespace ROOT;

// Function to compute momentum from px, py, pz (for RDataFrame Define)
float computeMomentum(const RVec<float>& px, const RVec<float>& py, const RVec<float>& pz, size_t index) {
    return sqrt(px[index] * px[index] + py[index] * py[index] + pz[index] * pz[index]);
}

int main() {
    // Enable multi-threading in ROOT (set number of threads, e.g., 4)
    ROOT::EnableImplicitMT(4); // Adjust the number of threads based on your CPU (e.g., 4, 8, etc.)
    cout << "Multi-threading enabled for RDataFrame processing." << endl;

    // Open the ROOT file and create an RDataFrame from the TTree
    TFile* inputFile = new TFile("output/particle_data.root", "READ");
    if (!inputFile->IsOpen()) {
        cerr << "Error: Could not open output/particle_data.root!" << endl;
        return 1;
    }

    // Create RDataFrame from the TTree
    ROOT::RDataFrame df("ParticleTree", inputFile);

    // Define the model for the 2D histogram
    TH2F hChi2PIDvsPPID211("hChi2PIDvsPPID211", "Chi2PID vs. Momentum for Positive Pions (PID 211);Momentum p (GeV);Chi2PID", 100, 0, 10, 100, -5, 5);

    // Process the data in parallel using RDataFrame
    // Step 1: Filter for PID 211 and valid chi2pid values
    auto filtered_df = df.Filter([](const RVec<int>& pid, const RVec<float>& chi2pid) {
        for (size_t i = 0; i < pid.size(); ++i) {
            if (pid[i] == 211 && chi2pid[i] != 9999.0) {
                return true; // At least one particle in the event is a valid pion
            }
        }
        return false;
    }, {"pid", "chi2pid"});

    // Step 2: Flatten the RVecs into individual entries and compute momentum
    auto flattened_df = filtered_df.Define("flat_idx", [](const RVec<int>& pid) {
        vector<size_t> indices;
        for (size_t i = 0; i < pid.size(); ++i) {
            if (pid[i] == 211) {
                indices.push_back(i);
            }
        }
        return indices;
    }, {"pid"})
    .Define("flat_chi2pid", [](const RVec<float>& chi2pid, const RVec<size_t>& indices) {
        vector<float> values;
        for (const auto& idx : indices) {
            if (chi2pid[idx] != 9999.0) {
                values.push_back(chi2pid[idx]);
            }
        }
        return values;
    }, {"chi2pid", "flat_idx"})
    .Define("flat_p", [](const RVec<float>& px, const RVec<float>& py, const RVec<float>& pz, const RVec<size_t>& indices) {
        vector<float> momenta;
        for (const auto& idx : indices) {
            momenta.push_back(computeMomentum(px, py, pz, idx));
        }
        return values;
    }, {"px", "py", "pz", "flat_idx"});

    // Step 3: Create the 2D histogram using Histo2D
    auto histo = flattened_df.Histo2D(
        {&hChi2PIDvsPPID211, "hChi2PIDvsPPID211"},
        "flat_p",
        "flat_chi2pid",
        "flat_chi2pid" // Use chi2pid as weight (ensures correct filling for flattened entries)
    );

    // Debug: Check entries in the histogram
    cout << "Entries in hChi2PIDvsPPID211: " << histo->GetEntries() << endl;

    // Plot: Chi2PID vs. Momentum for Positive Pions (PID 211) with logarithmic z-axis
    TCanvas* canvasChi2PIDvsP = new TCanvas("canvasChi2PIDvsP", "Chi2PID vs. Momentum for Positive Pions", 800, 600);
    if (histo->GetEntries() > 0) {
        histo->SetStats(0);
        histo->SetMinimum(1e-5); // Set minimum for log scale
        histo->Draw("COLZ");
        gPad->SetLogz(); // Set logarithmic z-axis
    } else {
        cout << "Warning: hChi2PIDvsPPID211 is empty, plot will be empty." << endl;
    }
    canvasChi2PIDvsP->Print("output_mt/chi2pid_vs_p_pion.pdf");

    // Clean up
    delete canvasChi2PIDvsP;
    inputFile->Close();
    delete inputFile;

    // Disable multi-threading (good practice)
    ROOT::DisableImplicitMT();

    cout << "Program completed successfully." << endl;
    return 0;
}