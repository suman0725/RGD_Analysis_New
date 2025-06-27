#include <TFile.h>
#include <TH1F.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
using namespace ROOT;

// chi2 PID cut functions
double getdtCutNeg(double p) {
    return -0.232 + (-6.405) * exp(-8.123 * p); 
}

double getdtCutPos(double p) {
    return 0.225 + 5.245 * exp(-7.298 * p); 
}

// Theoretical beta calculator
double getTheoreticalBeta(double p, double mass) {
    return p / sqrt(p * p + mass * mass);
}

int main() {
    ROOT::EnableImplicitMT(16);

    // Open input file
    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/Skim/pkptreeCxC_9_test_modified.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file.\n";
        ROOT::DisableImplicitMT();
        return 1;
    }

    // Load tree
    ROOT::RDataFrame df_pions("EB_pos_pion_assumed", file);

    // Momentum bin definition
    vector<double> pBins;
    double pMin = 1.0, pMax = 7.0;
    double binWidth = 0.3;
    int nBins = static_cast<int>((pMax - pMin) / binWidth);
    for (int i = 0; i <= nBins; i++) {
        pBins.push_back(pMin + i * binWidth);
    }

    // Use only first bin
    double pLow = pBins[0];
    double pHigh = pBins[1];
    double pCenter = 0.5 * (pLow + pHigh);

    // Filter data for that bin (before cut)
    auto df_pions_before = df_pions.Filter([pLow, pHigh](float p) {
        return p >= pLow && p < pHigh;
    }, {"p"});

    // Create beta histogram before cut
    ROOT::RDF::TH1DModel modelBetaBefore(
        "beta_before_0",
        TString::Format("p: [%.1f, %.1f) GeV/c; #beta; Counts", pLow, pHigh),
        100, 0.96, 1.03
    );
    auto histo_before = df_pions_before.Histo1D(modelBetaBefore, "beta");

    // Apply dt cut (after cut)
    auto df_pions_after = df_pions_before.Filter([](float p, float dt) {
        return dt > getdtCutNeg(p) && dt < getdtCutPos(p);
    }, {"p", "dt"});

    // Create beta histogram after cut
    ROOT::RDF::TH1DModel modelBetaAfter(
        "beta_after_0",
        TString::Format("p: [%.1f, %.1f) GeV/c; #beta; Counts", pLow, pHigh),
        100, 0.96, 1.03
    );
    auto histo_after = df_pions_after.Histo1D(modelBetaAfter, "beta");

    // Theoretical beta values at pCenter
    const double massPion = 0.13957;
    const double massKaon = 0.49367;
    const double massProton = 0.93827;
    double beta_pion = getTheoreticalBeta(pCenter, massPion);
    double beta_kaon = getTheoreticalBeta(pCenter, massKaon);
    double beta_proton = getTheoreticalBeta(pCenter, massProton);

    // Create canvas
    TCanvas* canvas = new TCanvas("canvas", "Beta Histograms with Theory Lines", 800, 600);
    gStyle->SetOptStat(0);

    // Draw histograms
    histo_before->SetLineColor(kBlack);
    histo_before->SetFillColorAlpha(kBlack, 0.3);
    histo_before->Draw("HIST");

    histo_after->SetLineColor(kBlue + 1);
    histo_after->SetFillColorAlpha(kBlue + 1, 0.3);
    histo_after->Draw("HIST SAME");

    // Get max Y for lines
    double maxY = histo_before->GetMaximum() * 1.1;

    // Add vertical lines
    TLine* line_pion = new TLine(beta_pion, 0, beta_pion, maxY);
    line_pion->SetLineColor(kBlue + 1);
    line_pion->SetLineStyle(2);  // dashed
    line_pion->SetLineWidth(2);
    line_pion->Draw();

    TLine* line_kaon = new TLine(beta_kaon, 0, beta_kaon, maxY);
    line_kaon->SetLineColor(kGreen + 2);
    line_kaon->SetLineStyle(2);  // dashed
    line_kaon->SetLineWidth(2);
    line_kaon->Draw();

    TLine* line_proton = new TLine(beta_proton, 0, beta_proton, maxY);
    line_proton->SetLineColor(kRed);
    line_proton->SetLineStyle(2);  // dashed
    line_proton->SetLineWidth(2);
    line_proton->Draw();

    // Add legend
    TLegend* legend = new TLegend(0.62, 0.65, 0.88, 0.88);
    legend->AddEntry(histo_before.GetPtr(), "Before Cut", "f");
    legend->AddEntry(histo_after.GetPtr(), "After Cut", "f");
    legend->AddEntry(line_pion, "Theory #beta (Pion)", "l");
    legend->AddEntry(line_kaon, "Theory #beta (Kaon)", "l");
    legend->AddEntry(line_proton, "Theory #beta (Proton)", "l");
    legend->Draw();

    // Save and clean up
    canvas->SaveAs("beta_histogram_with_theory_lines.png");

    file->Close();
    delete file;
    delete canvas;
    ROOT::DisableImplicitMT();

    cout << "Histogram saved to beta_histogram_with_theory_lines.png" << endl;
    return 0;
}
