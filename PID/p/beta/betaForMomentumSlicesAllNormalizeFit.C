#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <ROOT/TTreeProcessorMT.hxx> // For multi-threading
#include <TStopwatch.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <TF1.h>

void betaForMomentumSlicesAllNormalizeFit() {
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
    const size_t numSlices = static_cast<size_t>((sliceEnd - sliceStart) / sliceStep);

    const std::vector<float> betaMin = {
        0.50, 0.56, 0.62, 0.66, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80,
        0.82, 0.84, 0.86, 0.86, 0.88, 0.88, 0.90, 0.90, 0.90, 0.90,
        0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.94, 0.94, 0.94, 0.94,
        0.94
    };
   const std::vector<std::pair<float, float>> fitRangesLeft = {
        {0.50, 0.64}, {0.56, 0.69}, {0.63, 0.72}, {0.67, 0.78}, {0.70, 0.79},
        {0.73, 0.82}, {0.76, 0.83}, {0.78, 0.86}, {0.80, 0.87}, {0.82, 0.89},
        {0.84, 0.90}, {0.84, 0.93}, {0.86, 0.91}, {0.88, 0.92}, {0.88, 0.93},
        {0.89, 0.94}, {0.90, 0.94}, {0.91, 0.95}, {0.91, 0.95}, {0.92, 0.96},
        {0.93, 0.96}, {0.93, 0.96}, {0.935, 0.962}, {0.935, 0.964}, {0.94, 0.97},
        {0.94, 0.97}, {0.95, 0.97}, {0.95, 0.972}, {0.95, 0.974}, {0.955, 0.975},
        {0.955, 0.98}
    };
    const std::vector<std::pair<float, float>> fitRangesRight = {
        {0.92, 1.02}, {0.94, 1.02}, {0.94, 1.02}, {0.95, 1.02}, {0.95, 1.02},
        {0.95, 1.02}, {0.96, 1.02}, {0.96, 1.02}, {0.97, 1.02}, {0.97, 1.02},
        {0.97, 1.02}, {0.98, 1.02}, {0.98, 1.02}, {0.98, 1.02}, {0.98, 1.02},
        {0.985, 1.02}, {0.985, 1.02}, {0.985, 1.02}, {0.985, 1.02}, {0.9875, 1.02},
        {0.988, 1.02}, {0.99, 1.02}, {0.99, 1.02}, {0.99, 1.02}, {0.992, 1.02},
        {0.992, 1.02}, {0.992, 1.02}, {0.992, 1.02}, {0.993, 1.02}, {0.992, 1.02},
        {0.992, 1.02}
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
    canvas->Print("beta_histograms_all_fit.pdf[");

    // Open a text file to save fit results
    std::ofstream fitResultsFile("fit_results_two_gauss.txt");
    for (size_t i = 0; i < histograms.size(); ++i) {
        // Normalize the histogram
        if (histograms[i]->Integral() > 0) {
            histograms[i]->Scale(1.0 / histograms[i]->Integral());
        }

        // Fit two Gaussians to the histogram
        auto fitRangeLeft = fitRangesLeft[i];
        auto fitRangeRight = fitRangesRight[i];

        TF1* fitLeft = new TF1(Form("fitLeft_%zu", i), "gaus", fitRangeLeft.first, fitRangeLeft.second);
        TF1* fitRight = new TF1(Form("fitRight_%zu", i), "gaus", fitRangeRight.first, fitRangeRight.second);

        histograms[i]->Fit(fitLeft, "RQ");
        histograms[i]->Fit(fitRight, "RQ+");

        // Extract the fit parameters for both fits
        double meanLeft = fitLeft->GetParameter(1);
        double sigmaLeft = fitLeft->GetParameter(2);

        double meanRight = fitRight->GetParameter(1);
        double sigmaRight = fitRight->GetParameter(2);

        // Save the results in the text file
        
        float sliceMin = sliceStart + i * sliceStep;
        float sliceMax = sliceMin + sliceStep;

        // Assuming meanLeft, sigmaLeft, meanRight, sigmaRight are already computed for this slice

        // Save the results in the text file with direct values
        fitResultsFile << sliceMin << " " << sliceMax << " "
                    << meanLeft << " " << sigmaLeft << " "
                    << meanRight << " " << sigmaRight << "\n";
        


        // Set axis titles and y-axis maximum
        histograms[i]->GetXaxis()->SetTitle("Beta");
        histograms[i]->GetYaxis()->SetTitle("Normalized Counts");
        histograms[i]->SetMaximum(0.03);

        // Draw and save to PDF
        histograms[i]->Draw();
        fitLeft->SetLineColor(kBlue);
        fitLeft->Draw("same");
        fitRight->SetLineColor(kRed);
        fitRight->Draw("same");
        canvas->Print("beta_histograms_all_fit.pdf");
    }

    canvas->Print("beta_histograms_all_fit.pdf]");
    fitResultsFile.close();

    // Clean up
    for (auto hist : histograms) {
        delete hist;
    }

    file->Close();
    delete canvas;

    timer.Stop();
    std::cout << "Processing completed in " << timer.RealTime() << " seconds." << std::endl;
}
