#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TCanvas.h>
#include <ROOT/TTreeProcessorMT.hxx> // Include for multi-threading
#include <TStopwatch.h>
#include <cmath>
#include <ROOT/TThreadedObject.hxx> // For thread-safe histogram handling

void betaForMomentum() {
    TStopwatch timer;
    timer.Start();

    // Open the ROOT file
    TFile* file = TFile::Open("../charged_particles.root", "READ");
    if (!file || !file->IsOpen()) {
        std::cerr << "Error: Unable to open the ROOT file!" << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("charged_particle");
    if (!tree) {
        std::cerr << "Error: Unable to access the TTree!" << std::endl;
        return;
    }

    // Use TThreadedObject for thread-local histograms
    ROOT::TThreadedObject<TH2F> threaded_histo(
        "threaded_histo",
        "Beta vs Momentum;Momentum (GeV);Beta",
        200, 0, 8, 200, 0.88, 1.02);

    // Multi-threaded tree processing
    ROOT::TTreeProcessorMT processor(*tree);
    processor.Process([&](TTreeReader& reader) {
        TTreeReaderValue<int> charge(reader, "charge");
        TTreeReaderValue<short> status(reader, "status");
        TTreeReaderValue<float> px(reader, "px");
        TTreeReaderValue<float> py(reader, "py");
        TTreeReaderValue<float> pz(reader, "pz");
        TTreeReaderValue<float> beta(reader, "beta");

        // Access the thread-local histogram (shared_ptr)
        auto local_histo = threaded_histo.Get();

        while (reader.Next()) {
            // Check for forward detector and positive charge
            if ((abs(*status) / 2000 == 1) && (*charge > 0)) {
                float momentum = sqrt((*px) * (*px) + (*py) * (*py) + (*pz) * (*pz));
                local_histo->Fill(momentum, *beta); // Use -> for shared_ptr
            }
        }
    });

    // Merge thread-local histograms into the global histogram
    TH2F* h_beta_vs_momentum = new TH2F(
        "h_beta_vs_momentum",
        "Beta vs Momentum;Momentum (GeV);Beta",
        200, 0, 8, 200, 0.88, 1.02);

    threaded_histo.Merge([&](std::shared_ptr<TH2F> main, std::vector<std::shared_ptr<TH2F>>& histos) {
        for (auto& histo : histos) {
            main->Add(histo.get());
        }
        *h_beta_vs_momentum = *main; // Copy the merged histogram into the final one
    });

    // Draw the histogram to canvas
    TCanvas* canvas = new TCanvas("canvas", "Beta vs. Momentum", 800, 600);
    canvas->SetLogz();
    h_beta_vs_momentum->SetMarkerStyle(20);
    h_beta_vs_momentum->Draw("COLZ");
    h_beta_vs_momentum->SetStats(0); 
    canvas->SaveAs("output.pdf");

    // Close the ROOT file
    //file->Close();

    timer.Stop();
    std::cout << "Execution completed in " << timer.RealTime() << " seconds (real time), "
              << timer.CpuTime() << " seconds (CPU time)." << std::endl;
}
