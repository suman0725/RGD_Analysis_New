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
#include <vector>
#include <cmath>

void betaForMomentumSlicesPeakD() {
    const int numSlices = 59;
    const float sliceMin = 0.2; // Starting momentum
    const float sliceMax = 3.0; // Ending momentum
    const float sliceWidth = (sliceMax - sliceMin) / numSlices;

    TFile* file = TFile::Open("../charged_particles.root", "READ");
    TTree* tree = (TTree*)file->Get("charged_particle");

    ROOT::TTreeProcessorMT processor(*tree);

    for (int i = 0; i < numSlices; ++i) {
        float currentMin = sliceMin + i * sliceWidth;
        float currentMax = currentMin + sliceWidth;

        TH1F* h_beta = new TH1F(Form("h_beta_slice_%d", i),
                                Form("Beta for Momentum Slice [%.2f, %.2f] GeV", currentMin, currentMax),
                                200, 0.75, 1.05);

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
                        h_beta->Fill(*beta);
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

        // Fit Gaussians for each detected peak
        for (int j = 0; j < 3; ++j) {
            double peak = sortedPeaks[j];
            TF1* gauss = new TF1(Form("gauss_peak_%d_slice_%d", j, i), "gaus(0)", peak - 0.02, peak + 0.02);
            gauss->SetParameters(h_beta->GetBinContent(h_beta->FindBin(peak)), peak, 0.01);
            h_beta->Fit(gauss, "R");
        }

        // Save the histogram to a file or canvas
        TCanvas* canvas = new TCanvas(Form("canvas_slice_%d", i), "Momentum Slice", 800, 600);
        h_beta->Draw();

        // Open a multi-page PDF
        if (i == 0) {
            canvas->SaveAs("output_all_slices.pdf("); // Open the PDF on the first slice
        } else if (i == numSlices - 1) {
            canvas->SaveAs("output_all_slices.pdf)"); // Close the PDF on the last slice
        } else {
            canvas->SaveAs("output_all_slices.pdf"); // Add pages in between
        }

        delete h_beta;
        delete canvas;
            }

    file->Close();
}

