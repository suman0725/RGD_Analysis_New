#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <ROOT/TTreeProcessorMT.hxx> // For multi-threading
#include <TStopwatch.h>
#include <cmath>
#include <vector>

void betaForMomentumSlicesAllNormalize() {
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
        file->Close();
        return;
    }

    // Define momentum slice ranges and beta minimum values
    const float sliceStep = 0.1; // 100 MeV slices
    const float sliceStart = 0.6;
    const float sliceEnd = 3.6;
    const std::vector<float> betaMin = {
        0.50, 0.56, 0.62, 0.66, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80,
        0.82, 0.84, 0.86, 0.86, 0.88, 0.88, 0.90, 0.90, 0.90, 0.90,
        0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.94, 0.94, 0.94, 0.94,
        0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
        0.94, 0.94, 0.94, 0.94, 0.95, 0.94, 0.96, 0.96, 0.96, 0.96,
        0.96, 0.96, 0.96, 0.96, 0.96, 0.96
    };

    // Create histograms for each slice
    std::vector<TH1F*> histograms;

    for (size_t i = 0; i < betaMin.size(); ++i) {
        float sliceMin = sliceStart + i * sliceStep;
        float sliceMax = sliceMin + sliceStep;

        TH1F* h_beta_slice = new TH1F(
            Form("h_beta_slice_%zu", i),
            Form("Beta for Momentum Slice [%.1f, %.1f] GeV", sliceMin, sliceMax),
            200, betaMin[i], 1.02
        );
        histograms.push_back(h_beta_slice);
    }

    // Process the tree
    ROOT::TTreeProcessorMT processor(*tree);
    processor.Process([&](TTreeReader& reader) {
        TTreeReaderValue<int> charge(reader, "charge");
        TTreeReaderValue<short> status(reader, "status");
        TTreeReaderValue<float> px(reader, "px");
        TTreeReaderValue<float> py(reader, "py");
        TTreeReaderValue<float> pz(reader, "pz");
        TTreeReaderValue<float> beta(reader, "beta");

        while (reader.Next()) {
            // Check for forward detector and positive charge
            if ((abs(*status) / 2000 == 1) && (*charge > 0)) {
                float momentum = sqrt((*px) * (*px) + (*py) * (*py) + (*pz) * (*pz));

                for (size_t i = 0; i < betaMin.size(); ++i) {
                    float sliceMin = sliceStart + i * sliceStep;
                    float sliceMax = sliceMin + sliceStep;

                    if (momentum >= sliceMin && momentum < sliceMax) {
                        histograms[i]->Fill(*beta);
                        break;
                    }
                }
            }
        }
    });

    // Create a canvas and save normalized histograms
    TCanvas* canvas = new TCanvas("canvas", "Beta Histograms", 800, 600);
    canvas->Print("beta_histograms_all.pdf[");

    for (size_t i = 0; i < histograms.size(); ++i) {
        // Normalize the histogram
        if (histograms[i]->Integral() > 0) {
            histograms[i]->Scale(1.0 / histograms[i]->Integral());
        }

        // Set axis titles and y-axis maximum
        histograms[i]->GetXaxis()->SetTitle("Beta");
        histograms[i]->GetYaxis()->SetTitle("Normalized Counts");
        histograms[i]->SetMaximum(0.03);

        // Draw and save to PDF
        histograms[i]->Draw();
        canvas->Update();
        canvas->Print("beta_histograms_all.pdf");
    }

    canvas->Print("beta_histograms_all.pdf]");

    // Clean up
    for (auto hist : histograms) {
        delete hist;
    }

    file->Close();
    delete canvas;

    timer.Stop();
    std::cout << "Processing completed in " << timer.RealTime() << " seconds." << std::endl;
}
