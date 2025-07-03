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
#include <iomanip>

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

// Double-sided Crystal Ball function for π+ (right tail mirrored from left)
Double_t doubleSidedCrystalBallFunc(Double_t* x, Double_t* par) {
    Double_t t = (x[0] - par[1]) / par[2];
    Double_t absAlphaLeft = fabs(par[3]);
    Double_t absAlphaRight = fabs(par[3]);
    Double_t nLeft = par[4];
    Double_t nRight = par[4];
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

// Total function (π+ + kaons + linear background)
Double_t totalCrystalBallFunc(Double_t* x, Double_t* par) {
    Double_t cb1 = doubleSidedCrystalBallFunc(x, &par[0]);  // Pion
    Double_t cb2 = crystalBallFunc(x, &par[7]);             // Kaon
    Double_t bg = par[12] + par[13] * x[0];                 // Linear background
    return cb1 + cb2 + bg;
}

int main() {
    auto start = chrono::high_resolution_clock::now();
    gStyle->SetOptStat(0);

    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/Skim/pkptreeCxC_9_test_modified.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file" << endl;
        return 1;
    }

    ROOT::RDataFrame df_all("EB_all_pion_assumed", file);

    double pLow = 1.6, pHigh = 1.9;
    auto df_pions_before = df_all.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});

    auto df_after_cut = df_pions_before.Filter([](float p, float dt) {
        return dt > getdtCutNeg(p) && dt < getdtCutPos(p);
    }, {"p", "dt"});

    ROOT::RDF::TH1DModel modelBeta("beta_fit_cb", "p: [1.6-1.9) GeV/c; #beta; Counts", 200, 0.95, 1.03);
    auto histo_beta_before = df_pions_before.Histo1D(modelBeta, "beta");
    auto histo_beta_after = df_after_cut.Histo1D(modelBeta, "beta");

    double p_mid = (pLow + pHigh) / 2;
    double beta_pi_theory = p_mid / sqrt(p_mid * p_mid + 0.1396 * 0.1396);
    double beta_K_theory = p_mid / sqrt(p_mid * p_mid + 0.4937 * 0.4937);
    double beta_p_theory = p_mid / sqrt(p_mid * p_mid + 0.95383 * 0.95383);

    double binWidth = (1.03 - 0.95) / 200;
    int binPi = static_cast<int>((	beta_pi_theory - 0.95) / binWidth) + 1;
    int binK = static_cast<int>((beta_K_theory - 0.95) / binWidth) + 1;

    double amp_pi = histo_beta_before->GetBinContent(binPi) * 1.1;
    double amp_K = histo_beta_before->GetBinContent(binK) * 1.1;

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

    TF1* totalCB = new TF1("totalCB", totalCrystalBallFunc, 0.95, 1.03, 14);
    totalCB->SetNpx(1000);
    totalCB->SetRange(0.95, 1.03);
    Double_t initParams[14] = {
        amp_pi, beta_pi_theory, 0.006, 2.0, 3.0, 2.0, 3.0, // Pion
        amp_K, beta_K_theory, 0.006, 1.0, 2.0,             // Kaon
        10.0, 5.0                                          // Background
    };
    totalCB->SetParameters(initParams);

    for (int i = 0; i < 14; ++i)
        totalCB->SetParName(i, Form("p%d", i));

    // Pion
    totalCB->SetParLimits(0, 0, 250000); totalCB->SetParameter(1, beta_pi_theory); totalCB->SetParLimits(2, 0.002, 0.008);
    totalCB->SetParLimits(3, 0.1, 20.0); totalCB->SetParLimits(4, 1.0, 20.0);
    totalCB->FixParameter(5, totalCB->GetParameter(3));
    totalCB->FixParameter(6, totalCB->GetParameter(4));

    // Kaon
    totalCB->SetParLimits(7, 0, 10000); totalCB->FixParameter(8, beta_K_theory); totalCB->SetParLimits(9, 0.0025, 0.008);
    totalCB->SetParLimits(10, 0.1, 20.0); totalCB->SetParLimits(11, 1.0, 100.0);

    // Background
    totalCB->SetParLimits(12, 0, 1000); totalCB->SetParLimits(13, -500, 500);

    TFitResultPtr fitResult = histo_beta_before->Fit(totalCB, "RME+S");
    totalCB->SetLineColor(kBlack);
    totalCB->SetLineWidth(2);
    totalCB->Draw("SAME");

    TLine* line_pi = new TLine(beta_pi_theory, 0, beta_pi_theory, histo_beta_before->GetMaximum());
    line_pi->SetLineColor(kBlue+1);
    line_pi->SetLineStyle(2);
    line_pi->SetLineWidth(3);
    line_pi->Draw();

    TLine* line_K = new TLine(beta_K_theory, 0, beta_K_theory, histo_beta_before->GetMaximum());
    line_K->SetLineColor(kGreen);
    line_K->SetLineStyle(2);
    line_K->Draw();

    TLine* line_p = new TLine(beta_p_theory, 0, beta_p_theory, histo_beta_before->GetMaximum());
    line_p->SetLineColor(kRed);
    line_p->SetLineStyle(2);
    line_p->SetLineWidth(3);
    //line_p->Draw();

    TLegend* leg = new TLegend(0.65, 0.6, 0.95, 0.95);
    leg->SetTextSize(0.02);
    leg->AddEntry(histo_beta_before.GetPtr(), "+ve particles (Before 3#sigma dt Cut)", "f");
    leg->AddEntry(histo_beta_after.GetPtr(), "+ve pions (After 3#sigma dt Cut)", "f");
    leg->AddEntry(totalCB, "Total Fit", "l");
    leg->AddEntry(line_pi, "#beta_{#pi_theory}", "l");
    leg->AddEntry(line_K, "#beta_{K_theory}", "l");
    leg->AddEntry(line_p, "#beta_{p_theory}", "l");

    // Pion peak contribution
    TF1* cb_pion = new TF1("cb_pion", totalCrystalBallFunc, 0.95, 1.03, 14);
    cb_pion->SetNpx(1000);
    Double_t pionParams[14] = {0};
    for (int i = 0; i < 7; ++i) pionParams[i] = totalCB->GetParameter(i);
    pionParams[5] = pionParams[3];
    pionParams[6] = pionParams[4];
    for (int i = 7; i < 12; ++i) pionParams[i] = 0;
    pionParams[12] = totalCB->GetParameter(12);
    pionParams[13] = totalCB->GetParameter(13);
    cb_pion->SetParameters(pionParams);
    for (int i = 0; i < 14; ++i) cb_pion->FixParameter(i, pionParams[i]);
    cb_pion->SetLineColor(kBlue+2);
    cb_pion->SetLineWidth(2);
    cb_pion->Draw("SAME");

    // Kaon peak contribution
    TF1* cb_kaon = new TF1("cb_kaon", totalCrystalBallFunc, 0.95, 1.03, 14);
    cb_kaon->SetNpx(1000);
    Double_t kaonParams[14] = {0};
    for (int i = 7; i < 12; ++i) kaonParams[i] = totalCB->GetParameter(i);
    for (int i = 0; i < 7; ++i) kaonParams[i] = 0;
    kaonParams[12] = totalCB->GetParameter(12);
    kaonParams[13] = totalCB->GetParameter(13);
    cb_kaon->SetParameters(kaonParams);
    for (int i = 0; i < 14; ++i) cb_kaon->FixParameter(i, kaonParams[i]);
    cb_kaon->SetLineColor(kGreen);
    cb_kaon->SetLineWidth(2);
    cb_kaon->Draw("SAME");

    // Background component
    TF1* bg = new TF1("bg", "[0] + [1]*x", 0.95, 1.03);
    bg->SetNpx(1000);
    bg->SetParameter(0, totalCB->GetParameter(12));
    bg->SetParameter(1, totalCB->GetParameter(13));
    bg->SetLineColor(kOrange);
    bg->SetLineWidth(2);
    bg->Draw("SAME");

    // Calculate kaon contamination
    double kaonIntegral = cb_kaon->Integral(0.9594, 1.03);
    double pionIntegral = cb_pion->Integral(0.95, 1.03);
    double kaonContamination = (pionIntegral > 0) ? (kaonIntegral / pionIntegral) * 100.0 : 0.0;

    // Output to file
    ofstream outFile("fit_parameters_bin2.txt");
    if (outFile.is_open()) {
        outFile << "Fit Parameters for Crystal Ball Fit (Bin: p = 2.2 - 2.5 GeV/c)\n";
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
        outFile << "Kaon Contamination (Kaon Integral [0.9594, 1.03] / Pion Integral [0.95, 1.03]): " << fixed << setprecision(3) << kaonContamination << "%\n";
        outFile.close();
    } else {
        cerr << "Error: Cannot open fit_parameters_bin2.txt for writing" << endl;
    }

    leg->AddEntry(cb_pion, "Pion fit + background", "l");
    leg->AddEntry(cb_kaon, "Kaon fit + background", "l");
    leg->AddEntry(bg, "Background", "l");
    leg->Draw();

    canvas->Print("bin_2.png");

    delete leg;
    delete line_pi;
    delete line_K;
    delete line_p;
    delete totalCB;
    delete cb_pion;
    delete cb_kaon;
    delete bg;
    delete canvas;
    file->Close();
    delete file;

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Total program time: " << duration.count() / 1000.0 << " seconds" << endl;

    return 0;
}