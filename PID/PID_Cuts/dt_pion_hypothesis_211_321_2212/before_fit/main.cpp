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
#include <TStyle.h>

using namespace std;
using namespace ROOT;

int main() {
    // Start total program timer
    auto start = chrono::high_resolution_clock::now();

    // Step 1: Enable multi-threading
    ROOT::EnableImplicitMT(4); // Adjust based on your CPU
    cout << "Multi-threading enabled for RDataFrame processing." << endl;

    // Step 2: Set style to remove stat box
    gStyle->SetOptStat(0);

    // Step 3: Open the ROOT file and create RDataFrame
    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/Skim/pkptreeCxC_9_test_modified.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file pkptreeCxC_9_test_modified.root\n";
        ROOT::DisableImplicitMT();
        return 1;
    }
    ROOT::RDataFrame df_all("EB_all_pion_assumed", file); // All PID 211, 321, 2212

    // Step 4: Create output directories
    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/dt_plots", kTRUE);

    // Step 5: Set up momentum bins (1.0 to 7.0 GeV, width 0.3)
    vector<double> pBins;
    double pMin = 1.0, pMax = 7.0;
    int nBins = static_cast<int>((pMax - pMin) / 0.3); // 6 bins
    for (int i = 0; i <= nBins; i++) {
        pBins.push_back(pMin + i * 0.3);
    }

    // Step 6: Initialize canvas and print start
    TCanvas* canvas = new TCanvas("canvas", "Delta t Distribution per Momentum Bin", 1200, 800);
    canvas->Print("output/dt_plots/dt_distributions_linear.pdf[");

    // Step 7: Process each momentum bin
    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];

        // Filter for momentum range (no chi2pid cut)
        auto filtered_df = df_all.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});

        // Create histogram model for dt
        ROOT::RDF::TH1DModel modelDt(
            TString::Format("dt_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #Delta t (ns); Counts", pLow, pHigh),
            100, -2.0, 2.0 // Typical dt range, adjust if needed
        );

        // Fill histogram
        auto histo_dt = filtered_df.Histo1D(modelDt, "dt");

        // Plot histogram with linear scale
        canvas->Clear();
        histo_dt->SetMarkerStyle(20); // Filled circles
        histo_dt->SetMarkerSize(1.2);
        histo_dt->SetMarkerColor(kBlack);
        histo_dt->Draw("P"); // Plot as points

        // Add legend with entry counts
        TLegend* leg = new TLegend(0.7, 0.6, 0.9, 0.8); // Shifted down
        leg->AddEntry(histo_dt.GetPtr(), Form("Delta t (%.0f)", histo_dt->GetEntries()), "p");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);
        leg->Draw();

        // Update and save
        canvas->Update();
        canvas->Print(TString::Format("output/dt_plots/dt_distributions_linear.pdf"));
        delete leg;
    }

    // Step 8: Close PDF
    canvas->Print("output/dt_plots/dt_distributions_linear.pdf]");

    // Step 9: Clean up
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