#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <chrono>

using namespace std;
using namespace ROOT;

int main() {
    // Start total program timer
    auto start = chrono::high_resolution_clock::now();

    // Step 1: Enable multi-threading
    ROOT::EnableImplicitMT(4); // Adjust based on your CPU (e.g., 4, 8)
    cout << "Multi-threading enabled for RDataFrame processing." << endl;

    // Step 2: Open the ROOT file and create RDataFrame
    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/Skim/pkptreeCxC_9_test_modified.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file pkptreeCxC_9_test_modified.root\n";
        ROOT::DisableImplicitMT();
        return 1;
    }
    ROOT::RDataFrame df("EB_all_pion_assumed", file);
    ROOT::RDataFrame df_pion("EB_pid_pions", file);

    // Step 3: Create output directory
    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/2d_plots", kTRUE);

    // Step 4: Define manual ranges
    const double pMin = 0.0;
    const double pMax = 10.0;
    const double betaMin = 0.8;
    const double betaMax = 1.05;
    const double chi2Min = -15.0;
    const double chi2Max = 15.0;
    const int nBinsX = 200; // Number of bins for momentum (p)
    const int nBinsY_beta = 200; // Number of bins for beta
    const int nBinsY_chi2 = 200; // Number of bins for chi2pid

    // Step 5: Create 2D histogram models with manual ranges
    ROOT::RDF::TH2DModel modelBetaVsP(
        "hBetaVsP",
        "Beta vs Momentum (Before Cut); Momentum (GeV/c); #beta",
        nBinsX, pMin, pMax,
        nBinsY_beta, betaMin, betaMax
    );
    ROOT::RDF::TH2DModel modelChi2VsP(
        "hChi2VsP",
        "Chi2pid vs Momentum (Before Cut); Momentum (GeV/c); chi2pid",
        nBinsX, pMin, pMax,
        nBinsY_chi2, chi2Min, chi2Max
    );
    ROOT::RDF::TH2DModel modelChi2VsP_pion(
        "hChi2VsP_pion",
        "Chi2pid vs Momentum (Pions); Momentum (GeV/c); chi2pid",
        nBinsX, pMin, pMax,
        nBinsY_chi2, chi2Min, chi2Max
    );

    // Step 6: Fill 2D histograms
    auto hBetaVsP = df.Histo2D(modelBetaVsP, "p", "beta");
    auto hChi2VsP = df.Histo2D(modelChi2VsP, "p", "recomputed_chi2pid");
    //auto hChi2VsP_pion = df_pion.Histo2D(modelChi2VsP_pion, "p", "recomputed_chi2pid");
    auto hChi2VsP_pion = df_pion.Histo2D(modelChi2VsP_pion, "p", "orig_chi2pid");

    // Step 7: Create canvases and plot with log scale for z-axis
    TCanvas* canvasBeta = new TCanvas("canvasBeta", "Beta vs Momentum", 800, 600);
    hBetaVsP->Draw("COLZ"); // COLZ for color map
    canvasBeta->SetLogz(); // Set log scale for z-axis
    canvasBeta->Print("output/2d_plots/beta_vs_p_before_cut_logz.pdf");

    TCanvas* canvasChi2 = new TCanvas("canvasChi2", "Chi2pid vs Momentum", 800, 600);
    hChi2VsP->Draw("COLZ"); // COLZ for color map
    canvasChi2->SetLogz(); // Set log scale for z-axis
    canvasChi2->Print("output/2d_plots/chi2pid_vs_p_before_cut_logz.pdf");

    TCanvas* canvasChi2_pion = new TCanvas("canvasChi2_pion", "Chi2pid vs Momentum (Pions)", 800, 600);
    hChi2VsP_pion->Draw("COLZ"); // COLZ for color map
    canvasChi2_pion->SetLogz(); // Set log scale for z-axis
    canvasChi2_pion->Print("output/2d_plots/chi2pid_vs_p_pion_before_cut_logz.pdf");

    // Step 8: Clean up
    delete canvasBeta;
    delete canvasChi2;
    delete canvasChi2_pion;
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