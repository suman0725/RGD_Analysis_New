#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <chrono>
#include <cmath>

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

    // Step 2: Open the ROOT file and create RDataFrame
    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/Skim/pkptreeCxC_9_test_modified.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file pkptreeCxC_9_test_modified.root\n";
        ROOT::DisableImplicitMT();
        return 1;
    }
    ROOT::RDataFrame df("EB_all_pion_assumed", file);

    // Step 3: Create output directory
    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/2d_plots_after_cut", kTRUE);

    // Step 4: Define manual ranges
    const double pMin = 0.0;
    const double pMax = 10.0;
    const double betaMin = 0.95;
    const double betaMax = 1.05;
    const double chi2Min = -3.5;
    const double chi2Max = 3.5;
    const int nBinsX = 200; // Number of bins for momentum (p)
    const int nBinsY_beta = 200; // Number of bins for beta
    const int nBinsY_chi2 = 200; // Number of bins for chi2pid
    const double C = 0.9795; // Scaling factor
    const double chi2NegCut = -3.0 * C; // -3 * 0.9795 ≈ -2.9385

    // Step 5: Apply chi2pid cut
    auto df_cut = df.Filter([chi2NegCut](float p, float chi2) {
    double chi2Cut = getChi2Cut(p);
    return chi2 > chi2NegCut && chi2 < chi2Cut;
    }, {"p", "recomputed_chi2pid"});

    // Step 6: Create 2D histogram models with manual ranges
    ROOT::RDF::TH2DModel modelBetaVsP(
        "hBetaVsP",
        "Beta vs Momentum (After Cut); Momentum (GeV/c); #beta",
        nBinsX, pMin, pMax,
        nBinsY_beta, betaMin, betaMax
    );
    ROOT::RDF::TH2DModel modelChi2VsP(
        "hChi2VsP",
        "Chi2pid vs Momentum (After Cut); Momentum (GeV/c); chi2pid",
        nBinsX, pMin, pMax,
        nBinsY_chi2, chi2Min, chi2Max
    );

    // Step 7: Fill 2D histograms
    auto hBetaVsP = df_cut.Histo2D(modelBetaVsP, "p", "beta");
    auto hChi2VsP = df_cut.Histo2D(modelChi2VsP, "p", "recomputed_chi2pid");

    // Step 8: Create canvases and plot with log scale for z-axis
    TCanvas* canvasBeta = new TCanvas("canvasBeta", "Beta vs Momentum (After Cut)", 800, 600);
    hBetaVsP->Draw("COLZ"); // COLZ for color map
    canvasBeta->SetLogz(); // Set log scale for z-axis
    canvasBeta->Print("output/2d_plots_after_cut/beta_vs_p_after_cut_logz.pdf");

    TCanvas* canvasChi2 = new TCanvas("canvasChi2", "Chi2pid vs Momentum (After Cut)", 800, 600);
    hChi2VsP->Draw("COLZ"); // COLZ for color map
    canvasChi2->SetLogz(); // Set log scale for z-axis
    canvasChi2->Print("output/2d_plots_after_cut/chi2pid_vs_p_after_cut_logz.pdf");

    // Step 9: Clean up
    delete canvasBeta;
    delete canvasChi2;
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