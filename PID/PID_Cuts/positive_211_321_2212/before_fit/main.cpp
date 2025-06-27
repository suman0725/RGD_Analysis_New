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
#include <TLine.h>

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

    // Step 3: Open the ROOT file and create RDataFrames
    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/Skim/pkptreeCxC_9_test_modified.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file pkptreeCxC_9_test_modified.root\n";
        ROOT::DisableImplicitMT();
        return 1;
    }
    ROOT::RDataFrame df_pos_pion_assumed("EB_pos_pion_assumed", file); // All positive charge particles
    ROOT::RDataFrame df_all_pion_assumed("EB_all_pion_assumed", file); // All PID 211, 321, 2212
    ROOT::RDataFrame df_pions("EB_pid_pions", file); // Pion tree
    ROOT::RDataFrame df_kaons("EB_pid_kaons", file); // Kaon tree
    ROOT::RDataFrame df_protons("EB_pid_protons", file); // Proton tree

    // Step 4: Create output directories
    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/beta_plots", kTRUE);

    // Step 5: Set up momentum bins (1.0 to 7.0 GeV, width 0.3)
    vector<double> pBins;
    double pMin = 1.0, pMax = 7.0;
    int nBins = static_cast<int>((pMax - pMin) / 0.3); // 6 bins
    for (int i = 0; i <= nBins; i++) {
        pBins.push_back(pMin + i * 0.3);
    }

    // Step 6: Initialize canvas and print start
    TCanvas* canvas = new TCanvas("canvas", "Beta Distribution per Momentum Bin", 1200, 800);
    canvas->Print("output/beta_plots/beta_distributions_linear.pdf[");

    // Step 7: Process each momentum bin
    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];

        // Filter for momentum range (no chi2pid cut)
        auto filtered_df_pos = df_pos_pion_assumed.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});
        auto filtered_df_all = df_all_pion_assumed.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});
        auto filtered_df_pions = df_pions.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});
        auto filtered_df_kaons = df_kaons.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});
        auto filtered_df_protons = df_protons.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});

        // Create histogram models for beta
        ROOT::RDF::TH1DModel modelBetaPos(
            TString::Format("beta_pos_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #beta; Counts", pLow, pHigh),
            100, 0.96, 1.03
        );
        ROOT::RDF::TH1DModel modelBetaAll(
            TString::Format("beta_all_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #beta; Counts", pLow, pHigh),
            100, 0.96, 1.03
        );
        ROOT::RDF::TH1DModel modelBetaPions(
            TString::Format("beta_pions_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #beta; Counts", pLow, pHigh),
            100, 0.96, 1.03
        );
        ROOT::RDF::TH1DModel modelBetaKaons(
            TString::Format("beta_kaons_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #beta; Counts", pLow, pHigh),
            100, 0.96, 1.03
        );
        ROOT::RDF::TH1DModel modelBetaProtons(
            TString::Format("beta_protons_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #beta; Counts", pLow, pHigh),
            100, 0.96, 1.03
        );

        // Fill histograms
        auto histo_pos = filtered_df_pos.Histo1D(modelBetaPos, "beta"); // All positive charge particles
        auto histo_all = filtered_df_all.Histo1D(modelBetaAll, "beta"); // All PID 211, 321, 2212
        auto histo_pions = filtered_df_pions.Histo1D(modelBetaPions, "beta");
        auto histo_kaons = filtered_df_kaons.Histo1D(modelBetaKaons, "beta");
        auto histo_protons = filtered_df_protons.Histo1D(modelBetaProtons, "beta");

        // Plot data points with linear scale
        canvas->Clear();
        histo_pos->SetMarkerStyle(20); // Filled circles
        histo_pos->SetMarkerSize(1.2);
        histo_pos->SetMarkerColor(kOrange);
        histo_pos->Draw("P"); // All positive charge particles (background)
        histo_all->SetMarkerStyle(20);
        histo_all->SetMarkerSize(1.2);
        histo_all->SetMarkerColor(kBlack);
        histo_all->Draw("P SAME"); // Combined PID 211, 321, 2212
        histo_pions->SetMarkerStyle(20);
        histo_pions->SetMarkerSize(1.2);
        histo_pions->SetMarkerColor(kBlue);
        histo_pions->Draw("P SAME"); // Pions
        histo_kaons->SetMarkerStyle(20);
        histo_kaons->SetMarkerSize(1.2);
        histo_kaons->SetMarkerColor(kGreen);
        histo_kaons->Draw("P SAME"); // Kaons
        histo_protons->SetMarkerStyle(20);
        histo_protons->SetMarkerSize(1.2);
        histo_protons->SetMarkerColor(kRed); // Changed to red for protons
        histo_protons->Draw("P SAME"); // Protons

        // Add theoretical beta lines (approximate for midpoint p = (pLow + pHigh) / 2)
        double p_mid = (pLow + pHigh) / 2;
        double beta_pi = p_mid / sqrt(p_mid * p_mid + 0.1396 * 0.1396);
        double beta_K = p_mid / sqrt(p_mid * p_mid + 0.4937 * 0.4937);
        double beta_p = p_mid / sqrt(p_mid * p_mid + 0.9383 * 0.9383);
        TLine* line_pi = new TLine(beta_pi, 0, beta_pi, histo_pos->GetMaximum());
        line_pi->SetLineColor(kBlue); line_pi->SetLineStyle(2); line_pi->Draw(); // Blue for pion
        TLine* line_K = new TLine(beta_K, 0, beta_K, histo_pos->GetMaximum());
        line_K->SetLineColor(kGreen); line_K->SetLineStyle(2); line_K->Draw(); // Green for kaon
        TLine* line_p = new TLine(beta_p, 0, beta_p, histo_pos->GetMaximum());
        line_p->SetLineColor(kRed); line_p->SetLineStyle(2); line_p->Draw(); // Red for proton

        // Add legend with entry counts
        TLegend* leg = new TLegend(0.7, 0.6, 0.9, 0.8); // Shifted down
        leg->AddEntry(histo_pos.GetPtr(), Form("All +ve Charge : %.0f", histo_pos->GetEntries()), "p");
        leg->AddEntry(histo_all.GetPtr(), Form("All PID : %.0f", histo_all->GetEntries()), "p");
        leg->AddEntry(histo_pions.GetPtr(), Form("Pion (211) : %.0f", histo_pions->GetEntries()), "p");
        leg->AddEntry(histo_kaons.GetPtr(), Form("Kaons (321) %.0f", histo_kaons->GetEntries()), "p");
        leg->AddEntry(histo_protons.GetPtr(), Form("Protons (2212) : %.0f", histo_protons->GetEntries()), "p");
        leg->AddEntry(line_pi, "Theory #pi", "l");
        leg->AddEntry(line_K, "Theory K", "l");
        leg->AddEntry(line_p, "Theory p", "l");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);
        leg->Draw();

        // Update and save
        canvas->Update();
        canvas->Print(TString::Format("output/beta_plots/beta_distributions_linear.pdf"));
        delete leg;
        delete line_pi;
        delete line_K;
        delete line_p;
    }

    // Step 8: Close PDF
    canvas->Print("output/beta_plots/beta_distributions_linear.pdf]");

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