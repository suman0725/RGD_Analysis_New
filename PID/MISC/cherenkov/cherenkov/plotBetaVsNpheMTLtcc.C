#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TCanvas.h>
#include <ROOT/TTreeProcessorMT.hxx> // Include for multi-threading
#include <TStopwatch.h>
#include <cmath>

void plotBetaVsNpheMTLtcc() {
    TStopwatch timer;
    timer.Start();

    // Open the ROOT file
    TFile* file = TFile::Open("../LTCC/charged_particles_LTCC_HTCC.root", "READ");
    if (!file || !file->IsOpen()) {
        std::cerr << "Error: Unable to open the ROOT file!" << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("charged_particle");
    if (!tree) {
        std::cerr << "Error: Unable to access the TTree!" << std::endl;
        return;
    }

    // Create 2D histogram for nphe vs. momentum
    TH2F* h_nphe_vs_beta = new TH2F(
        "h_nphe_vs_Beta",
        "nphe vs. Beta for Particles with Forward Detector and Charge > 0;Momentum (GeV/c);nphe",
        200, 0, 1, 200, 0, 80);
    TH1F* h_nphe = new TH1F ("h_nphe", "nphe distribution for positive particles ; nphe", 200, 0, 150);

    // Multi-threaded tree processing
    ROOT::TTreeProcessorMT processor(*tree);

    processor.Process([&](TTreeReader& reader) {
        TTreeReaderValue<int> charge(reader, "charge");
        TTreeReaderValue<int> detector(reader, "detector");
        TTreeReaderValue<short> status(reader, "status");
        TTreeReaderValue<float> px(reader, "px");
        TTreeReaderValue<float> py(reader, "py");
        TTreeReaderValue<float> pz(reader, "pz");
        TTreeReaderValue<float> nphe(reader, "nphe");
        TTreeReaderValue<float> beta(reader, "beta");

        while (reader.Next()) {
            // Check for forward detector and positive charge
            if ((abs(*status) / 2000 == 1) && (*charge > 0) && (*detector == 16 )) {
             /* float momentum = sqrt((*px) * (*px) + (*py) * (*py) + (*pz) * (*pz));
                if (2 < momentum < 2.25) {
                    h_nphe_vs_beta->Fill(*beta, *nphe);
            } */
           h_nphe->Fill(*nphe);
            }
        }
    });

    // Create a canvas to draw the histogram
    /* TCanvas* canvas = new TCanvas("canvas", "nphe vs. beta", 800, 600);
    canvas->SetLogz(); // Set log scale for the z-axis
    h_nphe_vs_beta->SetMarkerStyle(20);
    h_nphe_vs_beta->Draw("COLZ");

    // Save the histogram to a file
    canvas->SaveAs("nphe_vs_beta_mt_ltcc.png"); */
    TCanvas* canvas1 = new TCanvas ("nphe", "nphe distribution", 800, 600);
    canvas1->SetLogz(); 
    h_nphe->Draw(); 
    canvas1->SaveAs("nphe for +ve particles");
    // Close the ROOT file
    //file->Close();

    timer.Stop();
    std::cout << "Execution completed in " << timer.RealTime() << " seconds (real time), "
              << timer.CpuTime() << " seconds (CPU time)." << std::endl;
}
