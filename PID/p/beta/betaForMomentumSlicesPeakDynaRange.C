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
#include <TSpectrum.h>
#include <algorithm>

void betaForMomentumSlicesPeakDynaRange() {
    const int numSlices = 59;
    const float sliceMin = 0.5; // Starting momentum
    const float sliceMax = 3.0; // Ending momentum
    const float sliceWidth = (sliceMax - sliceMin) / numSlices;

    TFile* file = TFile::Open("../charged_particles.root", "READ");
    TTree* tree = (TTree*)file->Get("charged_particle");

    ROOT::TTreeProcessorMT processor(*tree);

    for (int i = 0; i < numSlices; ++i) {
        float currentMin = sliceMin + i * sliceWidth;
        float currentMax = currentMin + sliceWidth;

        // Step 1: Determine beta range dynamically for the current slice
        float minBeta = 0.5;  // Starting point for the slice range
        float maxBeta = 1.0;  // Temporary maxBeta to find the peaks

        processor.Process([&](TTreeReader& reader) {
            TTreeReaderValue<int> charge(reader, "charge");
            TTreeReaderValue<short> status(reader, "status");
            TTreeReaderValue<float> px(reader, "px");
            TTreeReaderValue<float> py(reader, "py");
            TTreeReaderValue<float> pz(reader, "pz");
            TTreeReaderValue<float> beta(reader, "beta");

            while (reader.Next()) {
                if (abs(*status) / 2000 == 1 && *charge > 0) {
                    float momentum = sqrt((*px) * (*px) + (*py) * (*py) + (*pz) * (*pz));
                    if (momentum >= currentMin && momentum < currentMax) {
                        minBeta = std::min(minBeta, *beta); // Update the dynamic minBeta
                        maxBeta = std::max(maxBeta, *beta); // Update the dynamic maxBeta
                    }
                }
            }
        });

        // Cap maxBeta at 1.02, ensuring it never exceeds this value
        maxBeta = std::min(maxBeta, 1.02f);

        // Step 2: Create the histogram with the dynamic beta range
        TH1F* h_beta = new TH1F(Form("h_beta_slice_%d", i),
                                Form("Beta for Momentum Slice [%.2f, %.2f] GeV", currentMin, currentMax),
                                200, minBeta, maxBeta);

        // Step 3: Fill the histogram with the updated beta range
        processor.Process([&](TTreeReader& reader) {
            TTreeReaderValue<int> charge(reader, "charge");
            TTreeReaderValue<short> status(reader, "status");
            TTreeReaderValue<float> px(reader, "px");
            TTreeReaderValue<float> py(reader, "py");
            TTreeReaderValue<float> pz(reader, "pz");
            TTreeReaderValue<float> beta(reader, "beta");

            while (reader.Next()) {
                if (abs(*status) / 2000 == 1 && *charge > 0) {
                    float momentum = sqrt((*px) * (*px) + (*py) * (*py) + (*pz) * (*pz));
                    if (momentum >= currentMin && momentum < currentMax) {
                        h_beta->Fill(*beta); // Fill with the dynamic beta range
                    }
                }
            }
        });

        // Normalize the histogram
        h_beta->Scale(1.0 / h_beta->Integral());

        // Use TSpectrum for peak detection
        TSpectrum spectrum(3); // Looking for up to 3 peaks
        int nPeaks = spectrum.Search(h_beta, 2, "", 0.1); // Adjust sensitivity if needed
        double* peakPositions = spectrum.GetPositionX();

        // Sort peaks in ascending order for consistency
        std::vector<double> sortedPeaks(peakPositions, peakPositions + nPeaks);
        std::sort(sortedPeaks.begin(), sortedPeaks.end());

        // Step 4: Adjust the lower limit of the histogram dynamically based on the first peak
        float lowerLimit = minBeta; // Default is the starting point of the range
        if (sortedPeaks.size() > 0) {
            lowerLimit = sortedPeaks[0] - 0.02; // Start just before the first peak
        }

        // Cap the upper limit at 1.02, no matter what
        float upperLimit = 1.02f;

        // Step 5: Create and fit the histogram
        TH1F* h_beta_final = new TH1F(Form("h_beta_final_slice_%d", i),
                                      Form("Beta for Momentum Slice [%.2f, %.2f] GeV", currentMin, currentMax),
                                      200, lowerLimit, upperLimit);

        // Step 6: Fill the histogram with the updated beta range
        processor.Process([&](TTreeReader& reader) {
            TTreeReaderValue<int> charge(reader, "charge");
            TTreeReaderValue<short> status(reader, "status");
            TTreeReaderValue<float> px(reader, "px");
            TTreeReaderValue<float> py(reader, "py");
            TTreeReaderValue<float> pz(reader, "pz");
            TTreeReaderValue<float> beta(reader, "beta");

            while (reader.Next()) {
                if (abs(*status) / 2000 == 1 && *charge > 0) {
                    float momentum = sqrt((*px) * (*px) + (*py) * (*py) + (*pz) * (*pz));
                    if (momentum >= currentMin && momentum < currentMax) {
                        h_beta_final->Fill(*beta); // Fill the final histogram
                    }
                }
            }
        });

        // Normalize the final histogram
        h_beta_final->Scale(1.0 / h_beta_final->Integral());

        // Use TSpectrum for peak detection again
        nPeaks = spectrum.Search(h_beta_final, 2, "", 0.1);
        peakPositions = spectrum.GetPositionX();

        // Sort peaks in ascending order for consistency
        sortedPeaks.clear();
        for (int j = 0; j < nPeaks; ++j) {
            sortedPeaks.push_back(peakPositions[j]);
        }
        std::sort(sortedPeaks.begin(), sortedPeaks.end());

        // Fit Gaussians for each detected peak
        for (int j = 0; j < 3 && j < nPeaks; ++j) {
            double peak = sortedPeaks[j];
            TF1* gauss = new TF1(Form("gauss_peak_%d_slice_%d", j, i), "gaus(0)", peak - 0.02, peak + 0.02);
            gauss->SetParameters(h_beta_final->GetBinContent(h_beta_final->FindBin(peak)), peak, 0.01);
            h_beta_final->Fit(gauss, "R");
        }

        // Save the histogram to a file or canvas
        TCanvas* canvas = new TCanvas(Form("canvas_slice_%d", i), "Momentum Slice", 800, 600);
        h_beta_final->Draw();

        // Open a multi-page PDF
        if (i == 0) {
            canvas->SaveAs("output_all_slices_dynaRange.pdf("); // Open the PDF on the first slice
        } else if (i == numSlices - 1) {
            canvas->SaveAs("output_all_slices_dynaRange.pdf)"); // Close the PDF on the last slice
        } else {
            canvas->SaveAs("output_all_slices.pdf"); // Add pages in between
        }

        delete h_beta_final;
        delete canvas;
    }

    file->Close();
}

