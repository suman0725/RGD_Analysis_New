#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TLine.h>
#include <TF1.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cmath>

using namespace std;
using namespace ROOT;

// chi2 PIC cut functions (momentum-dependent)
double getdtCutNeg(double p) {
    return -0.232 + (-6.405) * exp(-8.123 * p); 
}

double getdtCutPos(double p) {
    return 0.225 + 5.245 * exp(-7.298 * p); 
}

// Calculate theoretical beta
double getTheoreticalBeta(double p, double mass) {
    return p / sqrt(p * p + mass * mass);
}

// Double Crystal Ball function
Double_t doubleCrystalBall(Double_t* x, Double_t* par) {
    Double_t t = (x[0] - par[1]) / par[2];
    Double_t absAlphaLeft = TMath::Abs(par[3]);
    Double_t absAlphaRight = TMath::Abs(par[6]);
    Double_t left = (t < -absAlphaLeft) ? 
        par[0] * TMath::Power(par[4] / absAlphaLeft, par[4]) * TMath::Exp(-0.5 * absAlphaLeft * absAlphaLeft) /
        TMath::Power(par[4] / absAlphaLeft - absAlphaLeft - t, par[4]) : 
        par[0] * TMath::Exp(-0.5 * t * t);
    Double_t right = (t > absAlphaRight) ? 
        par[5] * TMath::Power(par[7] / absAlphaRight, par[7]) * TMath::Exp(-0.5 * absAlphaRight * absAlphaRight) /
        TMath::Power(par[7] / absAlphaRight - absAlphaRight + t, par[7]) : 
        par[5] * TMath::Exp(-0.5 * t * t);
    return left + right;
}

int main() {
    auto start = chrono::high_resolution_clock::now();

    ROOT::EnableImplicitMT(16);
    cout << "Multi-threading enabled for RDataFrame processing." << endl;

    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/Skim/pkptreeCxC_9_test_modified.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file pkptreeCxC_9_test_modified.root\n";
        ROOT::DisableImplicitMT();
        return 1;
    }
    ROOT::RDataFrame df_pions("EB_pos_pion_assumed", file);
    gSystem->mkdir("output", kTRUE);

    vector<double> pBins;
    double pMin = 1.0, pMax = 7.0;
    int nBins = static_cast<int>((pMax - pMin) / 0.3);
    for (int i = 0; i <= nBins; i++) {
        pBins.push_back(pMin + i * 0.3);
    }

    TCanvas* canvas = new TCanvas("canvas", "Beta Distribution per Momentum Bin (1.0 to 7.0 GeV/c)", 1200, 800);
    canvas->Print("output/beta_distributions_1.0_to_7.0.pdf[");

    const double pionMass = 0.1396;
    const double kaonMass = 0.4937;
    const double protonMass = 0.9383;

    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];
        double pMid = (pLow + pHigh) / 2.0;

        // Unfiltered (before cut) data
        auto df_pions_before = df_pions.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});
        ROOT::RDF::TH1DModel modelBetaBefore(
            TString::Format("beta_before_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #beta; Counts", pLow, pHigh),
            100, 0.96, 1.03
        );
        auto histo_before = df_pions_before.Histo1D(modelBetaBefore, "beta");
        histo_before->SetStats(0);
        histo_before->Sumw2();
        histo_before->SetMarkerColor(kBlack);
        histo_before->SetMarkerStyle(20);
        histo_before->SetMarkerSize(1.2);
        histo_before->SetMarkerColorAlpha(kBlack, 0.5); // Fade before-cut data

        // Filtered (after cut) data
        auto filtered_df_pions = df_pions.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"})
                                            .Filter([](float p, float dt) {
                                                double dtMin = getdtCutNeg(p);
                                                double dtMax = getdtCutPos(p);
                                                return dt > dtMin && dt < dtMax;
                                            }, {"p", "dt"});
        ROOT::RDF::TH1DModel modelBetaAfter(
            TString::Format("beta_after_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #beta; Counts", pLow, pHigh),
            100, 0.96, 1.03
        );
        auto histo_after = filtered_df_pions.Histo1D(modelBetaAfter, "beta");
        histo_after->SetStats(0);
        histo_after->Sumw2();
        histo_after->SetMarkerColor(kBlue);
        histo_after->SetMarkerStyle(20);
        histo_after->SetMarkerSize(1.2);

        double betaPion = getTheoreticalBeta(pMid, pionMass);
        double betaKaon = getTheoreticalBeta(pMid, kaonMass);
        double betaProton = getTheoreticalBeta(pMid, protonMass);

        canvas->Clear();
        histo_before->Draw("P"); // Draw before-cut data first (faded black)
        histo_after->Draw("P SAME"); // Draw after-cut data on top (blue)

        double maxHeight = max(histo_before->GetMaximum(), histo_after->GetMaximum());
        TLine* linePion = new TLine(betaPion, 0, betaPion, maxHeight);
        linePion->SetLineColor(kBlue);
        linePion->SetLineStyle(2);
        linePion->SetLineWidth(2);
        linePion->SetLineColorAlpha(kBlue, 0.3);
        linePion->Draw();
        TLine* lineKaon = new TLine(betaKaon, 0, betaKaon, maxHeight);
        lineKaon->SetLineColor(kGreen);
        lineKaon->SetLineStyle(2);
        lineKaon->SetLineWidth(2);
        lineKaon->SetLineColorAlpha(kGreen, 0.3);
        lineKaon->Draw();
        TLine* lineProton = new TLine(betaProton, 0, betaProton, maxHeight);
        lineProton->SetLineColor(kRed);
        lineProton->SetLineStyle(2);
        lineProton->SetLineWidth(2);
        lineProton->SetLineColorAlpha(kRed, 0.3);
        lineProton->Draw();

        TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg->AddEntry(histo_before.GetPtr(), "Before Cut (+ve particles)", "p");
        leg->AddEntry(histo_after.GetPtr(), "After Cut (+ve particles)", "p");
        leg->AddEntry(linePion, "Pion #beta_{theory}", "l");
        leg->AddEntry(lineKaon, "Kaon #beta_{theory}", "l");
        leg->AddEntry(lineProton, "Proton #beta_{theory}", "l");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);
        leg->Draw();

        canvas->Update();
        canvas->Print(TString::Format("output/beta_distributions_1.0_to_7.0.pdf"));
        delete leg;
        delete linePion;
        delete lineKaon;
        delete lineProton;
    }

    canvas->Print("output/beta_distributions_1.0_to_7.0.pdf]");
    
    delete canvas;
    file->Close();
    delete file;

    ROOT::DisableImplicitMT();
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Total program time: " << duration.count() / 1000.0 << " seconds" << endl;

    cout << "Program completed successfully." << endl;
    return 0;
}