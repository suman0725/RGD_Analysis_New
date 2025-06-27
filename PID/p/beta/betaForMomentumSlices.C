#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <ROOT/TTreeProcessorMT.hxx> // For multi-threading
#include <TStopwatch.h>
#include <cmath>
#include <vector>

void betaForMomentumSlices() {
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

    // Define momentum slice ranges
    const float minMomentum = 0.6;
    const float maxMomentum = 6.5;
    const float sliceWidth = 0.1;
    const int numSlices = static_cast<int>((maxMomentum - minMomentum) / sliceWidth);

    // Create a vector to store histograms for each slice
    std::vector<TH1F*> histograms;

    // Create canvas for saving multiple plots to a single PDF
    TCanvas* canvas = new TCanvas("canvas", "Beta in Momentum Slices", 800, 600);
    canvas->Print("output_beta_slices.pdf["); // Open PDF file

    // Loop over momentum slices
    for (int i = 0; i < numSlices; i++) {
        float sliceMin = minMomentum + i * sliceWidth;
        float sliceMax = sliceMin + sliceWidth;

        // Create histogram for this slice
        std::string histName = "h_beta_slice_" + std::to_string(i);
        std::string histTitle = "Beta for Momentum Slice: [" + std::to_string(sliceMin) + ", " + std::to_string(sliceMax) + "] GeV";
        TH1F* h_beta_slice = new TH1F(histName.c_str(), histTitle.c_str(), 200, 0.5, 1.1); // Beta range starts from 0.5
        histograms.push_back(h_beta_slice);

        // Read and process the tree
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
                    if (momentum >= sliceMin && momentum < sliceMax) {
                        h_beta_slice->Fill(*beta);
                    }
                }
            }
        });

        // Draw the histogram to canvas
        canvas->cd();
        h_beta_slice->GetXaxis()->SetTitle("Beta");
        h_beta_slice->GetYaxis()->SetTitle("Counts");

         // Set y-axis to logarithmic scale
      /*  h_beta_slice->SetMaximum(h_beta_slice->GetMaximum() * 10); // Ensure the histogram has space for the log scale
        canvas->SetLogy(); 
 */
        h_beta_slice->Draw();
        canvas->Print("output_beta_slices.pdf"); // Add this plot to the PDF
    }

    canvas->Print("output_beta_slices.pdf]"); // Close PDF file

    // Clean up
    for (auto hist : histograms) {
        delete hist;
    }
    delete canvas;
    file->Close();

    timer.Stop();
    std::cout << "Execution completed in " << timer.RealTime() << " seconds (real time), "
              << timer.CpuTime() << " seconds (CPU time)." << std::endl;
}
