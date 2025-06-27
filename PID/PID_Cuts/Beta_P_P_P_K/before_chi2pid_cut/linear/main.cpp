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
#include <TLine.h>

using namespace std;
using namespace ROOT;

// Calculate theoretical beta
double getTheoreticalBeta(double p, double mass) {
    return p / sqrt(p * p + mass * mass);
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
    ROOT::RDataFrame df_kaons("EB_pid_kaons", file); // Kaon tree

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
    canvas->Print("output/beta_plots/beta_distributions_linear.pdf[");


    // Particle masses (GeV/c^2)
    const double pionMass = 0.1396; // Pion mass
    const double kaonMass = 0.4937; // Kaon mass


    // Step 6: Process each momentum bin
    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];
        double pMid = (pLow + pHigh) / 2.0;

        // Filter for momentum range (no chi2pid cut)
        auto filtered_df_pions = df_pions.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});
        auto filtered_df_kaons = df_kaons.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});

        // Create histogram models for beta (match your specified range)
        ROOT::RDF::TH1DModel modelBetaPions(
            TString::Format("beta_pions_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #beta; Counts", pLow, pHigh),
            100, 0.96, 1.03 // Match your earlier plot range
        );
        ROOT::RDF::TH1DModel modelBetaKaons(
            TString::Format("beta_kaons_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #beta; Counts", pLow, pHigh),
            100, 0.96, 1.03
        );

        // Fill histograms
        auto histo_pions = filtered_df_pions.Histo1D(modelBetaPions, "beta");
        auto histo_kaons = filtered_df_kaons.Histo1D(modelBetaKaons, "beta");

        // Plot data points with linear scale
        canvas->Clear();
        histo_pions->SetMarkerStyle(20); // Filled circles
        histo_pions->SetMarkerSize(1.2);
        histo_pions->SetMarkerColor(kBlue);
        histo_pions->Draw("P"); // Plot as points
        histo_kaons->SetMarkerStyle(20); // Filled circles
        histo_kaons->SetMarkerSize(1.2);
        histo_kaons->SetMarkerColor(kGreen);
        histo_kaons->Draw("P SAME"); // Overlay points

        // Remove log scale to use linear scale (default)
        // canvas->SetLogy(); // Commented out to ensure linear scale

         // Calculate theoretical beta
         double betaPion = getTheoreticalBeta(pMid, pionMass);
         double betaKaon = getTheoreticalBeta(pMid, kaonMass);

        // Get maximum height for vertical lines (use larger of the two histograms)
        double maxHeight = max(histo_pions->GetMaximum(), histo_kaons->GetMaximum()) * 1.1;

        // Add theoretical vertical lines
        TLine* linePion = new TLine(betaPion, 0, betaPion, maxHeight);
        linePion->SetLineColor(kBlue);
        linePion->SetLineStyle(2); // Dashed
        linePion->SetLineWidth(2);
        linePion->SetLineColorAlpha(kBlue, 0.3); // 30% opacity for fade
        linePion->Draw();

        TLine* lineKaon = new TLine(betaKaon, 0, betaKaon, maxHeight);
        lineKaon->SetLineColor(kGreen);
        lineKaon->SetLineStyle(2); // Dashed
        lineKaon->SetLineWidth(2);
        lineKaon->SetLineColorAlpha(kGreen, 0.3); // 30% opacity for fade
        lineKaon->Draw();

        // Add legend
        TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg->AddEntry(histo_pions.GetPtr(), "Pions (Uncut)", "p");
        leg->AddEntry(histo_kaons.GetPtr(), "Kaons (Uncut)", "p");
        leg->AddEntry(linePion, "Pion #beta_{theory}", "l");
        leg->AddEntry(lineKaon, "Kaon #beta_{theory}", "l");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);
        leg->Draw();

        // Update and save
        canvas->Update();
        canvas->Print(TString::Format("output/beta_plots/beta_distributions_linear.pdf"));
        delete leg;
    }

    // Step 7: Close PDF
    canvas->Print("output/beta_plots/beta_distributions_linear.pdf]");

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