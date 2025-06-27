#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

using namespace std;
using namespace ROOT;

// Function to compute the chi2pid cut value based on momentum
double getChi2Cut(double p) {
    const double C = 0.9795; // Scaling factor
    if (p <= 2.75) {
        return 3.0 * C; // 3.0 * 0.9795 ≈ 2.9385
    } else {
        return C * (-0.20 + 3.00 * exp(-p / 4.51) + 37.26 * exp(-p / 0.87));
    }
}

int main() {
    // Start total program timer
    auto start = chrono::high_resolution_clock::now();

    // Step 1: Enable multi-threading
    ROOT::EnableImplicitMT(4); // Adjust based on your CPU (e.g., 4, 8)
    cout << "Multi-threading enabled for RDataFrame processing." << endl;

    // Step 2: Open the ROOT file and create RDataFrames
    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/Skim/pkptreeCxC_9_test_modified.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file pkptreeCxC_9_test_modified.root\n";
        ROOT::DisableImplicitMT();
        return 1;
    }
    ROOT::RDataFrame df_pions("EB_pid_pions", file); // Pion tree
    ROOT::RDataFrame df_kaons("EB_pid_kaons", file); // Kaon tree (for reference, though not used here)

    // Step 3: Create output directories
    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/beta_plots", kTRUE);

    // Step 4: Set up momentum bins (1.0 to 7.0 GeV, width 0.3)
    vector<double> pBins;
    double pMin = 1.0, pMax = 7.0;
    int nBins = static_cast<int>((pMax - pMin) / 0.3); // 6 bins (1.0-1.3, 1.3-1.6, ..., 6.7-7.0)
    for (int i = 0; i <= nBins; i++) {
        pBins.push_back(pMin + i * 0.3);
    }

    // Step 5: Initialize canvas and print start
    TCanvas* canvas = new TCanvas("canvas", "Beta Distribution per Momentum Bin", 1200, 800);
    canvas->Print("output/beta_plots/beta_distributions_log.pdf[");

    const double C = 0.9795;
    const double chi2NegCut = -3.0 * C; // For pion selection

    // Step 6: Process each momentum bin
    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];

        // Filter for momentum range (uncut)
        auto filtered_df_pions_uncut = df_pions.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});

        // Filter for momentum range with chi2pid cut
        auto filtered_df_pions_cut = df_pions.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"})
                                            .Filter([chi2NegCut](float p, float chi2) {
                                                double chi2Cut = getChi2Cut(p);
                                                return chi2 > chi2NegCut && chi2 < chi2Cut;
                                            }, {"p", "orig_chi2pid"});

        // Create histogram models for beta
        ROOT::RDF::TH1DModel modelBetaPionsUncut(
            TString::Format("beta_pions_uncut_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #beta; Counts", pLow, pHigh),
            100, 0.96, 1.02 // Match your specified range
        );
        ROOT::RDF::TH1DModel modelBetaPionsCut(
            TString::Format("beta_pions_cut_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #beta; Counts", pLow, pHigh),
            100, 0.96, 1.02
        );

        // Fill histograms
        auto histo_pions_uncut = filtered_df_pions_uncut.Histo1D(modelBetaPionsUncut, "beta");
        auto histo_pions_cut = filtered_df_pions_cut.Histo1D(modelBetaPionsCut, "beta");

        // Disable stats box
        histo_pions_uncut->SetStats(0);
        histo_pions_cut->SetStats(0);

        // Plot data points with log scale
        canvas->Clear();
        histo_pions_uncut->SetMarkerStyle(20); // Filled circles
        histo_pions_uncut->SetMarkerSize(1.2);
        histo_pions_uncut->SetMarkerColor(kBlue);
        histo_pions_uncut->Draw("P"); // Plot as points
        histo_pions_cut->SetMarkerStyle(20); // Filled circles
        histo_pions_cut->SetMarkerSize(1.2);
        histo_pions_cut->SetMarkerColor(kRed);
        histo_pions_cut->Draw("P SAME"); // Overlay points

        // Set log scale
        canvas->SetLogy();

        // Add legend with entry counts
        TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg->AddEntry(histo_pions_uncut.GetPtr(), 
                      TString::Format("Pions (Uncut) Entries: %lld", (long long)histo_pions_uncut->GetEntries()), "p");
        leg->AddEntry(histo_pions_cut.GetPtr(), 
                      TString::Format("Pions (Cut) Entries: %lld", (long long)histo_pions_cut->GetEntries()), "p");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);
        leg->Draw();

        // Update and save
        canvas->Update();
        canvas->Print(TString::Format("output/beta_plots/beta_distributions_log.pdf"));
        delete leg;
    }

    // Step 7: Close PDF
    canvas->Print("output/beta_plots/beta_distributions_log.pdf]");

    // Step 8: Clean up
    delete canvas;
    file->Close();
    delete file;

    // Disable multi-threading
    ROOT::DisableImplicitMT();
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Total program time: " << duration.count() / 1000.0 << " seconds" << endl;

    cout << "Program completed successfully." << endl;
    return 0;
}