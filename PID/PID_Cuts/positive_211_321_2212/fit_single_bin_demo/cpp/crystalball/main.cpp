/* 
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>
#include <TLine.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <vector>
#include <chrono>
#include <TStyle.h>
#include <TMath.h>

using namespace std;
using namespace ROOT;


Double_t crystalBallFunc(Double_t* x, Double_t* par) {
    // par[0] = amplitude, par[1] = mean, par[2] = sigma, par[3] = alpha, par[4] = n
    Double_t t = (x[0] - par[1]) / par[2];
    if (par[3] < 0) t = -t;

    Double_t absAlpha = fabs(par[3]);
    Double_t result = 0.0;

    if (t >= -absAlpha) {
        result = par[0] * exp(-0.5 * t * t);
    } else {
        Double_t a = pow(par[4] / absAlpha, par[4]) * exp(-0.5 * absAlpha * absAlpha);
        Double_t b = par[4] / absAlpha - absAlpha;
        result = par[0] * a / pow(b - t, par[4]);
    }

    return result;
}
// Crystal Ball function definition
Double_t totalCrystalBallFunc(Double_t* x, Double_t* par) {
    // Each CB function has 5 parameters
    Double_t cb1 = crystalBallFunc(x, &par[0]);   // params 0–4
    Double_t cb2 = crystalBallFunc(x, &par[5]);   // params 5–9
    Double_t cb3 = crystalBallFunc(x, &par[10]);  // params 10–14
    Double_t poly = par[15] + par[16]*x[0];       // background

    return cb1 + cb2 + cb3 + poly;
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

    gSystem->mkdir("output/fit_plots_cb", kTRUE);

    double pLow = 4.0, pHigh = 4.3;
    auto filtered_df = df_all.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});

    ROOT::RDF::TH1DModel modelBeta("beta_fit_cb", "p: [4.0-4.3) GeV/c; #beta; Counts", 100, 0.95, 1.03);
    auto histo_beta = filtered_df.Histo1D(modelBeta, "beta");

    double p_mid = (pLow + pHigh) / 2;
    double beta_pi = p_mid / sqrt(p_mid * p_mid + 0.1396 * 0.1396);
    double beta_K = p_mid / sqrt(p_mid * p_mid + 0.4937 * 0.4937);
    double beta_p = p_mid / sqrt(p_mid * p_mid + 0.9383 * 0.9383);

    int bin_pi = histo_beta->FindBin(beta_pi);
    int bin_K = histo_beta->FindBin(beta_K);
    int bin_p = histo_beta->FindBin(beta_p);
    double amp_pi = histo_beta->GetBinContent(bin_pi);
    double amp_K = histo_beta->GetBinContent(bin_K);
    double amp_p = histo_beta->GetBinContent(bin_p);

    TCanvas* canvas = new TCanvas("canvas_cb", "Crystal Ball Fit", 1400, 1000);
    canvas->SetGridx(1);
    canvas->SetGridy(1);
    histo_beta->SetMarkerStyle(20);
    histo_beta->SetMarkerSize(1.2);
    histo_beta->SetMarkerColor(kBlack);
    histo_beta->Draw("P");

   TF1* totalCB = new TF1("totalCB", totalCrystalBallFunc, 0.95, 1.03, 17);

    totalCB->SetNpx(1000);

    
    Double_t initParams[17] = {
    amp_pi, beta_pi, 0.003, 1.0, 2.0,   // π+
    amp_K, beta_K,  0.003, 1.0, 2.0,   // K+
    amp_p, beta_p,  0.003, 1.0, 2.0,   // p
    1.0, 0.0                          // background: constant + linear term
};
totalCB->SetParameters(initParams);


    for (int i = 0; i < 17; ++i)
        totalCB->SetParName(i, Form("p%d", i));

    histo_beta->Fit(totalCB, "RME+");
    totalCB->SetLineColor(kBlack);
    totalCB->Draw("SAME");

    TLine* line_pi = new TLine(beta_pi, 0, beta_pi, histo_beta->GetMaximum());
    TLine* line_K = new TLine(beta_K, 0, beta_K, histo_beta->GetMaximum());
    TLine* line_p = new TLine(beta_p, 0, beta_p, histo_beta->GetMaximum());
    for (auto* line : {line_pi, line_K, line_p}) {
        line->SetLineColor(kMagenta);
        line->SetLineStyle(2);
        line->Draw();
    }

    TLegend* leg = new TLegend(0.65, 0.55, 0.85, 0.75);
    leg->AddEntry(histo_beta.GetPtr(), "Data", "p");
    leg->AddEntry(totalCB, "Crystal Ball Total Fit", "l");
    leg->AddEntry(line_pi, "#beta_{theory}", "l");
    leg->Draw();

    canvas->Print("output/fit_plots_cb/beta_cbfit_4.0-4.3.pdf");
    canvas->Print("output/fit_plots_cb/beta_cbfit_4.0-4.3.png");

    delete leg;
    delete line_pi;
    delete line_K;
    delete line_p;
    delete totalCB;
    delete canvas;
    file->Close();
    delete file;

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Total program time: " << duration.count() / 1000.0 << " seconds" << endl;

    return 0;
}
 */
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>
#include <TLine.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <vector>
#include <chrono>
#include <TStyle.h>
#include <TMath.h>

using namespace std;
using namespace ROOT;

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

Double_t totalCrystalBallFunc(Double_t* x, Double_t* par) {
    Double_t cb1 = crystalBallFunc(x, &par[0]);   // params 0–4 for π+
    Double_t cb2 = crystalBallFunc(x, &par[5]);   // params 5–9 for K+
    Double_t poly = par[10] + par[11] * x[0];     // linear background
    return cb1 + cb2 + poly;
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

    gSystem->mkdir("output/fit_plots_cb", kTRUE);

    double pLow = 4.0, pHigh = 4.3;
    auto filtered_df = df_all.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});

    ROOT::RDF::TH1DModel modelBeta("beta_fit_cb", "p: [4.0-4.3) GeV/c; #beta; Counts", 100, 0.95, 1.03);
    auto histo_beta = filtered_df.Histo1D(modelBeta, "beta");

    // Theoretical beta values
    double p_mid = (pLow + pHigh) / 2;
    double beta_pi_theory = p_mid / sqrt(p_mid * p_mid + 0.1396 * 0.1396);
    double beta_K_theory = p_mid / sqrt(p_mid * p_mid + 0.4937 * 0.4937);

    // Calculate bin numbers for theoretical beta values
    int binPi = static_cast<int>((beta_pi_theory - 0.95) / 0.0008) + 1;
    int binK = static_cast<int>((beta_K_theory - 0.95) / 0.0008) + 1;

    // Get amplitudes from histogram at theoretical beta values
    double amp_pi = histo_beta->GetBinContent(binPi);
    double amp_K = histo_beta->GetBinContent(binK);

    // Optional scaling to approximate peak height
    amp_pi *= 1.1;
    amp_K *= 1.1;

    TCanvas* canvas = new TCanvas("canvas_cb", "Crystal Ball Fit", 1400, 1000);
    canvas->SetGridx(1);
    canvas->SetGridy(1);
    histo_beta->SetMarkerStyle(20);
    histo_beta->SetMarkerSize(1.2);
    histo_beta->SetMarkerColor(kBlack);
    histo_beta->Draw("P");

    TF1* totalCB = new TF1("totalCB", totalCrystalBallFunc, 0.95, 1.03, 12); // 10 params for 2 CBs + 2 for background

    totalCB->SetNpx(1000);

    Double_t initParams[12] = {
        amp_pi, beta_pi_theory, 0.005, 1.0, 2.0, // π+ (amp, mean, sigma, alpha, n)
        amp_K, beta_K_theory, 0.005, 1.0, 2.0,   // K+ (amp, mean, sigma, alpha, n)
        1000.0, 0.0                              // background: constant, linear term
    };
    totalCB->SetParameters(initParams);

    for (int i = 0; i < 12; ++i)
        totalCB->SetParName(i, Form("p%d", i));

    // Set parameter limits
    totalCB->SetParLimits(0, 0, 250000);  // Norm1 (π+)
    totalCB->SetParLimits(1, 0.95, 1.03); // Mean1 (π+)
    totalCB->SetParLimits(2, 0.001, 0.1); // Sigma1 (π+)
    totalCB->SetParLimits(3, 0.1, 10.0);  // Alpha1 (π+)
    totalCB->SetParLimits(4, 1.0, 50.0);  // n1 (π+)

    totalCB->SetParLimits(5, 0, 250000);  // Norm2 (K+)
    totalCB->SetParLimits(6, 0.95, 1.03); // Mean2 (K+)
    totalCB->SetParLimits(7, 0.001, 0.1); // Sigma2 (K+)
    totalCB->SetParLimits(8, 0.1, 10.0);  // Alpha2 (K+)
    totalCB->SetParLimits(9, 1.0, 50.0);  // n2 (K+)

    totalCB->SetParLimits(10, 0, 10000);  // Background constant
    totalCB->SetParLimits(11, -1000, 1000); // Background linear term

    histo_beta->Fit(totalCB, "RME+");
    totalCB->SetLineColor(kBlack);
    totalCB->Draw("SAME");

    TLine* line_pi = new TLine(beta_pi_theory, 0, beta_pi_theory, histo_beta->GetMaximum());
    TLine* line_K = new TLine(beta_K_theory, 0, beta_K_theory, histo_beta->GetMaximum());
    for (auto* line : {line_pi, line_K}) {
        line->SetLineColor(kMagenta);
        line->SetLineStyle(2);
        line->Draw();
    }

    TLegend* leg = new TLegend(0.65, 0.55, 0.85, 0.75);
    leg->AddEntry(histo_beta.GetPtr(), "Data", "p");
    leg->AddEntry(totalCB, "Crystal Ball Total Fit", "l");
    leg->AddEntry(line_pi, "#beta_{theory}", "l");
    leg->Draw();

    canvas->Print("output/fit_plots_cb/beta_cbfit_4.0-4.3.pdf");
    canvas->Print("output/fit_plots_cb/beta_cbfit_4.0-4.3.png");

    delete leg;
    delete line_pi;
    delete line_K;
    delete totalCB;
    delete canvas;
    file->Close();
    delete file;

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Total program time: " << duration.count() / 1000.0 << " seconds" << endl;

    return 0;
}

