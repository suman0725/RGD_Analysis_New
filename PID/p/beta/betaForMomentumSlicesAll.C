#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <ROOT/TTreeProcessorMT.hxx> // For multi-threading
#include <TStopwatch.h>
#include <cmath>
#include <vector>

void betaForMomentumSlicesAll() {
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

    // Define momentum slice range
    const float sliceStart = 0.2;
    const float sliceEnd = 6.5;
    const float sliceStep = 0.1;

    // Define beta minimums for each slice
    std::vector<float> betaMin = {0.15, 0.25, 0.35, 0.45,  0.5, 0.56, 0.62, 0.66, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80,
                                   0.82, 0.84, 0.86, 0.86, 0.88, 0.88, 0.90, 0.90, 0.90, 0.90,
                                   0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.94, 0.94, 0.94, 0.94,
                                   0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94,
                                   0.94, 0.94, 0.94, 0.95, 0.94, 0.96, 0.96, 0.96, 0.96, 0.96,
                                   0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96};

    // Prepare histograms for each slice
    std::vector<TH1F*> histograms;
    for (size_t i = 0; i < betaMin.size(); ++i) {
        float sliceMin = sliceStart + i * sliceStep;
        float sliceMax = sliceMin + sliceStep;

        // Create histogram for this slice
        std::string histName = "h_beta_slice_" + std::to_string(int(sliceMin * 100)) + "_" + std::to_string(int(sliceMax * 100));
        std::string histTitle = "Beta for Momentum Slice [" + std::to_string(sliceMin) + ", " + std::to_string(sliceMax) + "] GeV";
        histograms.push_back(new TH1F(histName.c_str(), histTitle.c_str(), 200, betaMin[i], 1.02));
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

                // Fill the appropriate histogram
                for (size_t i = 0; i < histograms.size(); ++i) {
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

    // Save histograms to a PDF
    TCanvas* canvas = new TCanvas("canvas", "Beta for Momentum Slices", 800, 600);
    canvas->Print("output_beta_slices_all.pdf["); // Open PDF file

    for (size_t i = 0; i < histograms.size(); ++i) {
        histograms[i]->Draw();
        canvas->Print("output_beta_slices_all.pdf");
    }

    canvas->Print("output_beta_slices_all.pdf]"); // Close PDF file

    // Clean up
    for (auto hist : histograms) {
        delete hist;
    }
    delete canvas;
    file->Close();

    timer.Stop();
    std::cout << "Processing completed in " << timer.RealTime() << " seconds." << std::endl;
}
