#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <ROOT/TTreeProcessorMT.hxx> // Include for multi-threading
#include <TStopwatch.h>
#include <iostream>
#include <cmath>

using namespace std;

void plotBetaMTNormalize() {
    TStopwatch timer;
    timer.Start();

    // Open the ROOT file
    TFile* file = TFile::Open("charged_particles.root", "READ");
    TTree* tree = (TTree*)file->Get("charged_particle");

    // Create a PDF to store the histograms
    TCanvas* canvas = new TCanvas("canvas", "Beta Distributions", 1200, 800);
    canvas->Print("beta_histograms_mt_normalize.pdf[");

    // Momentum ranges
    const float pStep = 0.25;
    const float pMax = 8.5;

    // Multi-threaded tree processing
    ROOT::TTreeProcessorMT processor(*tree);
    for (float pLow = 2.0; pLow < pMax; pLow += pStep) {
        float pHigh = pLow + pStep;

        // Histogram for this momentum range
        TH1F* hist = new TH1F(Form("hist_beta_%.2f_%.2f", pLow, pHigh),
                              Form("Beta Distribution for %.2f < p < %.2f", pLow, pHigh),
                              200, 0.88, 1.02);

        processor.Process([&](TTreeReader& reader) {
            TTreeReaderValue<int> charge(reader, "charge");
            TTreeReaderValue<short> status(reader, "status");
            TTreeReaderValue<float> px(reader, "px");
            TTreeReaderValue<float> py(reader, "py");
            TTreeReaderValue<float> pz(reader, "pz");
            TTreeReaderValue<float> beta(reader, "beta");

            while (reader.Next()) {
                if (abs(*status) / 2000 == 1 && *charge > 0) {
                    float p = sqrt((*px) * (*px) + (*py) * (*py) + (*pz) * (*pz));
                    if (p >= pLow && p < pHigh) {
                        hist->Fill(*beta);
                    }
                }
            }
        });

        // Normalize the histogram by scaling it to unit integral
        hist->Scale(1.0 / hist->Integral());

        // Plot the normalized histogram
        hist->Draw();
        canvas->Print("beta_histograms_mt_normalize.pdf");

        delete hist; // Cleanup
    }

    // Close the PDF and ROOT file
    canvas->Print("beta_histograms_mt_normalize.pdf]");
    file->Close();

    timer.Stop();
    cout << "Histograms saved to beta_histogram_mt_normalize.pdf" << endl;
    cout << "Total time: " << timer.RealTime() << " seconds (real time), "
         << timer.CpuTime() << " seconds (CPU time)" << endl;
}
