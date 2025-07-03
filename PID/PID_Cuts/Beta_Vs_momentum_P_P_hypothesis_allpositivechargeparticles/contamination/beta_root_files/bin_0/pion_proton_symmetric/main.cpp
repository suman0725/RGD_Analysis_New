#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>
#include <TLine.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFitResult.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <TStyle.h>
#include <TMath.h>
#include <fstream>

using namespace std;
using namespace ROOT;

// Crystal Ball function
Double_t crystalBallFunc(Double_t* x, Double_t* par) {
    Double_t t = (x[0] - par[1]) / par[2];
    Double_t absAlpha = fabs(par[3]);
    Double_t n = par[4];
    Double_t result = 0.0;
    if (t >= -absAlpha) {
        result = exp(-0.5 * t * t); 
    } else {
        Double_t A = pow(n / absAlpha, n) * exp(-0.5 * absAlpha * absAlpha);
        Double_t B = n / absAlpha - absAlpha;
        result = A / pow(B - t, n);
    }
    return par[0] * result;
}

// Double-sided Crystal Ball function (left tail mirrored from right)
Double_t doubleSidedCrystalBallFunc(Double_t* x, Double_t* par) {
    Double_t t = (x[0] - par[1]) / par[2];
    Double_t absAlphaRight = fabs(par[3]);
    Double_t absAlphaLeft = fabs(par[3]); // Left tail mirrors right
    Double_t nRight = par[4];
    Double_t nLeft = par[4];              // Left tail mirrors right
    Double_t result = 0.0;

    if (t >= -absAlphaLeft && t <= absAlphaRight) {
        result = exp(-0.5 * t * t);
    } else if (t < -absAlphaLeft) {
        Double_t A = pow(nLeft / absAlphaLeft, nLeft) * exp(-0.5 * absAlphaLeft * absAlphaLeft);
        Double_t B = nLeft / absAlphaLeft - absAlphaLeft;
        result = A / pow(B - t, nLeft);
    } else {
        Double_t A = pow(nRight / absAlphaRight, nRight) * exp(-0.5 * absAlphaRight * absAlphaRight);
        Double_t B = nRight / absAlphaRight - absAlphaRight;
        result = A / pow(t - B, nRight);
    }
    return par[0] * result;
}

double getdtCutNeg(double p) {
    return -0.232 + (-6.405) * exp(-8.123 * p);
}
double getdtCutPos(double p) {
    return 0.225 + 5.245 * exp(-7.298 * p);
}

// Total function (π+ only, no kaons, no background)
Double_t totalCrystalBallFunc(Double_t* x, Double_t* par) {
    Double_t cb1 = doubleSidedCrystalBallFunc(x, &par[0]);  // Pion (left tail mirrored from right)
    return cb1;
}

int main() {
    auto start = chrono::high_resolution_clock::now();
    gStyle->SetOptStat(0);

    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/new_RGD_Analysis/PID/charge_particles_custompid/Misidentification/Skim/pkptreeCxC_9_test_modified.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file" << endl;
        return 1;
    }

   

 ROOT::RDataFrame df_all("EB_all_pion_assumed", file);

    double pLow = 1.0, pHigh = 1.3;
    auto df_pions_before = df_all.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});

    auto df_after_cut = df_pions_before.Filter([](float p, float dt) {
        return dt > getdtCutNeg(p) && dt < getdtCutPos(p);
    }, {"p", "dt"});

    ROOT::RDF::TH1DModel modelBeta("beta_fit_cb", "p: [1.0-1.3) GeV/c; #beta; Counts", 200, 0.95, 1.03);
    auto histo_beta_before = df_pions_before.Histo1D(modelBeta, "beta");
    auto histo_beta_after = df_after_cut.Histo1D(modelBeta, "beta");

    // Calculate the true average momentum in the bin
    auto p_mean_action = df_pions_before.Mean("p");
    double p_mid = *p_mean_action;
    
    double beta_pi_theory = p_mid / sqrt(p_mid * p_mid + 0.1396 * 0.1396);

    double binWidth = (1.03 - 0.95) / 200;
    int binPi = static_cast<int>((beta_pi_theory - 0.95) / binWidth) + 1;

    double amp_pi = histo_beta_before->GetBinContent(binPi) * 1.1;

    TCanvas* canvas = new TCanvas("canvas_cb", "Crystal Ball Fit", 1400, 1000);
    canvas->SetGridx(1);
    canvas->SetGridy(1);

    histo_beta_before->SetFillColorAlpha(kBlack, 0.3);
    histo_beta_before->SetLineColor(kBlack);
    histo_beta_before->Draw("HIST");

    histo_beta_after->SetFillColorAlpha(kBlue, 0.3);
    histo_beta_after->SetLineColor(kBlue+1);
    histo_beta_after->SetLineWidth(2);
    histo_beta_after->Draw("HIST SAME");

    TF1* totalCB = new TF1("totalCB", totalCrystalBallFunc, 0.95, 1.03, 5); // 5 parameters (pion only)
    totalCB->SetNpx(1000);
    totalCB->SetRange(0.95, 1.03);
    Double_t initParams[5] = {
        amp_pi, beta_pi_theory, 0.003, 1.5, 3.0          // Pion (amp, mean, sigma, alphaRight, nRight)
    };
    totalCB->SetParameters(initParams);

    for (int i = 0; i < 5; ++i)
        totalCB->SetParName(i, Form("p%d", i));
    // Pion
    totalCB->SetParLimits(0, 0, 210000); totalCB->SetParameter(1, beta_pi_theory); totalCB->SetParLimits(2, 0.002, 0.004);
    totalCB->SetParLimits(3, 0.1, 10.0); totalCB->SetParLimits(4, 1.0, 20.0);

    TFitResultPtr fitResult = histo_beta_before->Fit(totalCB, "RME+S");
    totalCB->SetLineColor(kBlack);
    totalCB->SetLineWidth(4);
    totalCB->Draw("SAME");

    TLine* line_pi = new TLine(beta_pi_theory, 0, beta_pi_theory, histo_beta_before->GetMaximum());
    line_pi->SetLineColor(kBlue+1);
    line_pi->SetLineStyle(2);
    line_pi->SetLineWidth(3);
    line_pi->Draw();

    TLegend* leg = new TLegend(0.65, 0.6, 0.9, 0.9);
    leg->SetTextSize(0.02);
    leg->AddEntry(histo_beta_before.GetPtr(), "+ve particles (Before 3#sigma dt Cut)", "f");
    leg->AddEntry(histo_beta_after.GetPtr(), "+ve pions (After 3#sigma dt Cut)", "f");
    leg->AddEntry(totalCB, "Total Fit", "l");
    leg->AddEntry(line_pi, "#beta_{#pi_theory}", "l");
    leg->Draw(); 

    // Pion peak contribution
    TF1* cb_pion = new TF1("cb_pion", totalCrystalBallFunc, 0.95, 1.03, 5);
    cb_pion->SetNpx(1000);
    Double_t pionParams[5] = {0};
    for (int i = 0; i < 5; ++i) pionParams[i] = totalCB->GetParameter(i);
    cb_pion->SetParameters(pionParams);
    for (int i = 0; i < 5; ++i) cb_pion->FixParameter(i, pionParams[i]);
    cb_pion->SetLineColor(kBlue+2);
    cb_pion->SetLineWidth(2);
    cb_pion->Draw("SAME");

    int minBin = histo_beta_after->FindFirstBinAbove(0);
    int maxBin = histo_beta_after->FindLastBinAbove(0);
    double minBeta = histo_beta_after->GetXaxis()->GetBinLowEdge(minBin);
    double maxBeta = histo_beta_after->GetXaxis()->GetBinUpEdge(maxBin);
    cout << "Min Beta after -3sigma dt cut: " << minBeta << endl;
    cout << "Max Beta after +3sigma dt cut: " << maxBeta << endl;

    // Output to file
    ofstream outFile("fit_parameters_bin9.txt");
    if (outFile.is_open()) {
        outFile << "Fit Parameters for Crystal Ball Fit (Bin 5: p = 3.7 - 4.0 GeV/c)\n";
        outFile << "------------------------------------\n";
        for (int i = 0; i < totalCB->GetNpar(); ++i) {
            outFile << totalCB->GetParName(i) << ": " << totalCB->GetParameter(i)
                    << " +/- " << totalCB->GetParError(i) << "\n";
        }
        outFile << "------------------------------------\n";
        outFile << "NDF: " << fitResult->Ndf() << "\n";
        outFile << "Chi2: " << fitResult->Chi2() << "\n";
        outFile << "Chi2/NDF: " << (fitResult->Ndf() > 0 ? fitResult->Chi2() / fitResult->Ndf() : 0.0) << "\n";
        outFile << "EDM: " << fitResult->Edm() << "\n";
        outFile << "------------------------------------\n";
        outFile.close();
    } else {
        cerr << "Error: Cannot open fit_parameters_bin9.txt for writing" << endl;
    }

    leg->AddEntry(cb_pion, "Pion fit", "l");

    canvas->Print("bin_9.png");

    delete leg;
    delete line_pi;
    delete totalCB;
    delete cb_pion;
    delete canvas;
    file->Close();
    delete file;

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Total program time: " << duration.count() / 1000.0 << " seconds" << endl;

    return 0;
}