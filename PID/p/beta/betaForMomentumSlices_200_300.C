#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <ROOT/TTreeProcessorMT.hxx> // For multi-threading
#include <TStopwatch.h>
#include <cmath>
#include <vector>
#include <TF1.h>

void betaForMomentumSlices_200_300() {
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

    // Define the specific momentum slice range: 1.4 to 1.5 GeV
    const float sliceMin = 0.2;
    const float sliceMax = 0.3;

    // Create a vector to store histograms for this slice
    std::vector<TH1F*> histograms;

    // Create histogram for this slice
    TH1F* h_beta_slice = new TH1F("h_beta_slice_0.2_0.3", "Beta for Momentum Slice [0.2, 0.3] GeV", 200, 0.5, 1.5); // Beta range from 0.75 to 1.05
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

    // Normalize the histogram (divide by its integral)
    TH1F* h_beta_slice_normalized = (TH1F*)h_beta_slice->Clone("h_beta_slice_normalized");
    h_beta_slice_normalized->Scale(1.0 / h_beta_slice_normalized->Integral());  // Normalize the histogram

    // Create Canvas 2 (Normalized)
    TCanvas* canvas_normalized = new TCanvas("canvas_normalized", "Beta for Momentum Slice [0.2, 0.3] GeV (Normalized)", 800, 600);
    canvas_normalized->Print("output_beta_0.2_0.3_normalized.pdf["); // Open PDF file
    canvas_normalized->cd();

    h_beta_slice_normalized->GetXaxis()->SetTitle("Beta");
    h_beta_slice_normalized->GetYaxis()->SetTitle("Normalized Counts");
    h_beta_slice_normalized->SetMaximum(0.03);
    h_beta_slice_normalized->Draw();

    // Fit Gaussian functions to the three regions
    // 1st peak: 0.80 to 0.88
     TF1* gauss1 = new TF1("gauss1", "gaus(0)", 0.7, 0.98); // Define a Gaussian between 0.80 and 0.88
    gauss1->SetParameters(100, 0.84, 0.02);  // Initial guess for amplitude, mean, and sigma
    int fitStatus1 = h_beta_slice_normalized->Fit(gauss1, "R");
    if (fitStatus1 != 0) {
        std::cerr << "Fit for gauss1 failed!" << std::endl;
    } else {
        std::cout << "Fit for gauss1 successful!" << std::endl;
    }
    std::cout << "Gaussian 1: Mean = " << gauss1->GetParameter(1) 
              << ", Sigma = " << gauss1->GetParameter(2) << std::endl;
    /*
    // 2nd peak: 0.93 to 0.96
    TF1* gauss2 = new TF1("gauss2", "gaus(0)", 0.93, 0.96); // Define a Gaussian between 0.93 and 0.96
    gauss2->SetParameters(100, 0.94, 0.02);  // Initial guess for amplitude, mean, and sigma
    int fitStatus2 = h_beta_slice_normalized->Fit(gauss2, "R");
    if (fitStatus2 != 0) {
        std::cerr << "Fit for gauss2 failed!" << std::endl;
    } else {
        std::cout << "Fit for gauss2 successful!" << std::endl;
    }
    std::cout << "Gaussian 2: Mean = " << gauss2->GetParameter(1) 
              << ", Sigma = " << gauss2->GetParameter(2) << std::endl;


    // Overlay the fits
    gauss1->Draw("same");
    gauss2->Draw("same");
     */

    // Add the plot to the PDF
    canvas_normalized->Print("output_beta_0.2_0.3_normalized.pdf");
    canvas_normalized->Print("output_beta_0.2_0.3_normalized.pdf]"); // Close PDF file for normalized

    // Clean up
    for (auto hist : histograms) {
        delete hist;
    }
    delete canvas_normalized;
    file->Close();

    timer.Stop();
    std::cout << "Execution completed in " << timer.RealTime() << " seconds (real time), "
              << timer.CpuTime() << " seconds (CPU time)." << std::endl;
}
