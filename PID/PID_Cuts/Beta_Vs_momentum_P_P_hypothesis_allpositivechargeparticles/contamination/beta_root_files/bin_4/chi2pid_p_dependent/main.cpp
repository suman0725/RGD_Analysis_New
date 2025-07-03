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
#include <TVirtualFitter.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <TStyle.h>
#include <TMath.h>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace ROOT;

// Crystal Ball function (single-sided)
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

// Double-sided Crystal Ball function for π+ (symmetric tails)
Double_t doubleSidedCrystalBallFunc(Double_t* x, Double_t* par) {
    Double_t t = (x[0] - par[1]) / par[2];
    Double_t absAlpha = fabs(par[3]);
    Double_t n = par[4];
    Double_t result = 0.0;

    if (t >= -absAlpha && t <= absAlpha) {
        result = exp(-0.5 * t * t);
    } else if (t < -absAlpha) {
        Double_t A = pow(n / absAlpha, n) * exp(-0.5 * absAlpha * absAlpha);
        Double_t B = n / absAlpha - absAlpha;
        result = A / pow(B - t, n);
    } else {
        Double_t A = pow(n / absAlpha, n) * exp(-0.5 * absAlpha * absAlpha);
        Double_t B = n / absAlpha - absAlpha;
        result = A / pow(t - B, n);
    }
    return par[0] * result;
}

// Timing cut functions
double getdtCutNeg(double p) {
    return -0.232 + (-6.405) * exp(-8.123 * p);
}
double getdtCutPos(double p) {
    return 0.225 + 5.245 * exp(-7.298 * p);
}

// Chi2PID cut function
double getChi2PIDCutPos(double p) {
    const double sigma = 0.9795;
    if (p <= 2.15) {
        return 5.0 * sigma; // 5σ for p ≤ 2.15 GeV/c
    } else {
        return (-0.21 + 3.07 * exp(-p / 4.51) + 36.51 * exp(-p / 0.87)) * sigma; // For p > 2.15 GeV/c
    }
}
double getChi2PIDCutNeg(double p) {
    const double sigma = 0.9795;
    return 5.0 * sigma; // 5σ for all p
}

// Total function (π+ + kaons + linear background)
Double_t totalCrystalBallFunc(Double_t* x, Double_t* par) {
    Double_t cb1 = doubleSidedCrystalBallFunc(x, &par[0]); // Pion
    Double_t cb2 = crystalBallFunc(x, &par[7]);            // Kaon
    Double_t bg = par[12] + par[13] * x[0];                // Linear background
    return cb1 + cb2 + bg;
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

    double pLow = 2.2, pHigh = 2.5;
    auto df_pions_before = df_all.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});

    // Calculate the true average momentum in the bin
    auto p_mean_action = df_pions_before.Mean("p");
    double p_mid = *p_mean_action;

    // Timing cut (dt)
    auto df_after_dt_cut = df_pions_before.Filter([](float p, float dt) {
        return dt > getdtCutNeg(p) && dt < getdtCutPos(p);
    }, {"p", "dt"});

    // Chi2PID cut
    auto df_after_chi2pid_cut = df_pions_before.Filter([pLow](float p, float recomputed_chi2pid) {
        return recomputed_chi2pid > -getChi2PIDCutNeg(pLow) && recomputed_chi2pid < getChi2PIDCutPos(pLow);
    }, {"p", "recomputed_chi2pid"});

    ROOT::RDF::TH1DModel modelBeta("beta_fit_cb", "p: [2.2-2.5) GeV/c; #beta; Counts", 200, 0.95, 1.03);
    auto histo_beta_before = df_pions_before.Histo1D(modelBeta, "beta");
    auto histo_beta_after_dt = df_after_dt_cut.Histo1D(modelBeta, "beta");
    auto histo_beta_after_chi2pid = df_after_chi2pid_cut.Histo1D(modelBeta, "beta");

    // Check histogram content
    if (histo_beta_before->GetEntries() == 0) {
        cerr << "Error: Histogram histo_beta_before is empty!" << endl;
        file->Close();
        delete file;
        return 1;
    }

    // Determine the effective β range from the dt cut histogram
    double beta_min_after_dt = histo_beta_after_dt->GetBinCenter(histo_beta_after_dt->FindFirstBinAbove(0.001 * histo_beta_after_dt->GetMaximum(), 1));
    double beta_max_after_dt = histo_beta_after_dt->GetBinCenter(histo_beta_after_dt->FindLastBinAbove(0.001 * histo_beta_after_dt->GetMaximum(), 200));

    // Determine the effective β range for chi2pid cut histogram
    double beta_min_after_chi2pid = histo_beta_after_chi2pid->GetBinCenter(histo_beta_after_chi2pid->FindFirstBinAbove(0.001 * histo_beta_after_chi2pid->GetMaximum(), 1));
    double beta_max_after_chi2pid = histo_beta_after_chi2pid->GetBinCenter(histo_beta_after_chi2pid->FindLastBinAbove(0.001 * histo_beta_after_chi2pid->GetMaximum(), 200));

    double beta_pi_theory = p_mid / sqrt(p_mid * p_mid + 0.1396 * 0.1396);
    double beta_K_theory = p_mid / sqrt(p_mid * p_mid + 0.4937 * 0.4937);
    double beta_p_theory = p_mid / sqrt(p_mid * p_mid + 0.9383 * 0.9383);

    double binWidth = (1.03 - 0.95) / 200;
    int binPi = static_cast<int>((beta_pi_theory - 0.95) / binWidth) + 1;
    int binK = static_cast<int>((beta_K_theory - 0.95) / binWidth) + 1;
    int binP = static_cast<int>((beta_p_theory - 0.95) / binWidth) + 1;

    double amp_pi = histo_beta_before->GetBinContent(binPi) * 1.1;
    double amp_K = histo_beta_before->GetBinContent(binK) * 1.1;

    TCanvas* canvas = new TCanvas("canvas_cb", "Crystal Ball Fit", 1400, 1000);
    canvas->SetGridx(1);
    canvas->SetGridy(1);

    histo_beta_before->SetFillColorAlpha(kBlack, 0.3);
    histo_beta_before->SetLineColor(kBlack);
    histo_beta_before->Draw("HIST");

    histo_beta_after_dt->SetFillColorAlpha(kBlue, 0.3);
    histo_beta_after_dt->SetLineColor(kBlue+1);
    histo_beta_after_dt->SetLineWidth(2);
    histo_beta_after_dt->Draw("HIST SAME");

    histo_beta_after_chi2pid->SetFillColorAlpha(kRed, 0.3);
    histo_beta_after_chi2pid->SetLineColor(kRed+1);
    histo_beta_after_chi2pid->SetLineWidth(2);
    histo_beta_after_chi2pid->Draw("HIST SAME");

    // Fit setup
    TF1* totalCB = new TF1("totalCB", totalCrystalBallFunc, 0.95, 1.03, 14);
    totalCB->SetNpx(1000);
    totalCB->SetRange(0.95, 1.03);
    Double_t initParams[14] = {
        amp_pi, beta_pi_theory, 0.003, 2.0, 3.0, 2.0, 3.0, // Pion
        amp_K, beta_K_theory, 0.0025, 1.0, 2.0,             // Kaon
        10.0, 5.0                                          // Background
    };
    totalCB->SetParameters(initParams);

    for (int i = 0; i < 14; ++i)
        totalCB->SetParName(i, Form("p%d", i));

    // Pion parameters
    totalCB->SetParLimits(0, 0, 65000);
    totalCB->FixParameter(1, beta_pi_theory);
    totalCB->SetParLimits(2, 0.002, 0.006);
    totalCB->SetParLimits(3, 0.1, 20.0);
    totalCB->SetParLimits(4, 1.0, 20.0);
    totalCB->FixParameter(5, totalCB->GetParameter(3));
    totalCB->FixParameter(6, totalCB->GetParameter(4));

    // Kaon parameters
    totalCB->SetParLimits(7, 0, 10000);
    totalCB->FixParameter(8, beta_K_theory);
    totalCB->SetParLimits(9, 0.0025, 0.006);
    totalCB->SetParLimits(10, 0.1, 20.0);
    totalCB->SetParLimits(11, 1.0, 100.0);

    // Background parameters
    totalCB->SetParLimits(12, 0, 1000);
    totalCB->SetParLimits(13, -100, 100);

    // Configure minimizer
    TVirtualFitter::SetMaxIterations(10000);
    TVirtualFitter::SetPrecision(0.01);

    // Perform fit
    TFitResultPtr fitResult = histo_beta_before->Fit(totalCB, "RME+S", "", 0.95, 1.03);
    cout << "Fit status: " << fitResult->Status() << endl;
    if (!fitResult->IsValid()) {
        cout << "Fit failed, using approximate errors" << endl;
    } else {
        cout << "Fit converged successfully!" << endl;
    }
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
    line_K->SetLineWidth(3);
    line_K->Draw();

    TLine* line_p = new TLine(beta_p_theory, 0, beta_p_theory, histo_beta_before->GetMaximum());
    line_p->SetLineColor(kRed);
    line_p->SetLineStyle(2);
    line_p->SetLineWidth(3);
    // line_p->Draw();

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

    TLegend* leg = new TLegend(0.65, 0.6, 0.95, 0.95);
    leg->SetTextSize(0.02);
    leg->AddEntry(histo_beta_before.GetPtr(), "+ve particles (Before Cuts)", "f");
    leg->AddEntry(histo_beta_after_dt.GetPtr(), "+ve pions (After 3#sigma dt Cut)", "f");
    leg->AddEntry(histo_beta_after_chi2pid.GetPtr(), "+ve pions (After #chi^{2}_{PID} Cut)", "f");
    leg->AddEntry(totalCB, "Total Fit", "l");
    leg->AddEntry(line_pi, "#beta_{#pi_theory}", "l");
    leg->AddEntry(line_K, "#beta_{K_theory}", "l");
    leg->AddEntry(cb_pion, "Pion fit + background", "l");
    leg->AddEntry(cb_kaon, "Kaon fit + background", "l");
    leg->AddEntry(bg, "Background", "l");
    leg->Draw();

    // Calculate integrals and errors for dt cut
    double pion_integral_before = cb_pion->Integral(0.95, 1.03) / binWidth;
    double pion_integral_error_before = cb_pion->IntegralError(0.95, 1.03, totalCB->GetParameters(), fitResult->GetCovarianceMatrix().GetMatrixArray()) / binWidth;
    if (pion_integral_error_before <= 0 || !fitResult->IsValid()) {
        cout << "Warning: Pion integral error before (dt cut) is zero or negative, using sqrt(N) approximation" << endl;
        pion_integral_error_before = sqrt(pion_integral_before);
    }
    double pion_integral_after_dt = cb_pion->Integral(beta_min_after_dt, beta_max_after_dt) / binWidth;
    double pion_integral_error_after_dt = cb_pion->IntegralError(beta_min_after_dt, beta_max_after_dt, totalCB->GetParameters(), fitResult->GetCovarianceMatrix().GetMatrixArray()) / binWidth;
    if (pion_integral_error_after_dt <= 0 || !fitResult->IsValid()) {
        cout << "Warning: Pion integral error after (dt cut) is zero or negative, using sqrt(N) approximation" << endl;
        pion_integral_error_after_dt = sqrt(pion_integral_after_dt);
    }
    double kaon_integral_after_dt = cb_kaon->Integral(0.9881, 1.0081) / binWidth;
    double kaon_integral_error_after_dt = cb_kaon->IntegralError(0.9881, 1.0081, totalCB->GetParameters(), fitResult->GetCovarianceMatrix().GetMatrixArray()) / binWidth;
    if (kaon_integral_error_after_dt <= 0 || !fitResult->IsValid()) {
        cout << "Warning: Kaon integral error after (dt cut) is zero or negative, using sqrt(N) approximation" << endl;
        kaon_integral_error_after_dt = sqrt(kaon_integral_after_dt);
    }
    double total_after_dt = histo_beta_after_dt->Integral();
    double total_after_error_dt = sqrt(total_after_dt);

    // Calculate efficiency and contamination for dt cut
    double efficiency_dt = pion_integral_after_dt / pion_integral_before;
    double efficiency_error_dt = efficiency_dt * sqrt(pow(pion_integral_error_after_dt / pion_integral_after_dt, 2) + pow(pion_integral_error_before / pion_integral_before, 2));
    double kaon_contamination_dt = kaon_integral_after_dt / total_after_dt;
    double kaon_contamination_error_dt = kaon_contamination_dt * sqrt(pow(kaon_integral_error_after_dt / kaon_integral_after_dt, 2) + pow(total_after_error_dt / total_after_dt, 2));

    // Calculate integrals and errors for chi2pid cut
    double pion_integral_after_chi2pid = cb_pion->Integral(beta_min_after_chi2pid, beta_max_after_chi2pid) / binWidth;
    double pion_integral_error_after_chi2pid = cb_pion->IntegralError(beta_min_after_chi2pid, beta_max_after_chi2pid, totalCB->GetParameters(), fitResult->GetCovarianceMatrix().GetMatrixArray()) / binWidth;
    if (pion_integral_error_after_chi2pid <= 0 || !fitResult->IsValid()) {
        cout << "Warning: Pion integral error after (chi2pid cut) is zero or negative, using sqrt(N) approximation" << endl;
        pion_integral_error_after_chi2pid = sqrt(pion_integral_after_chi2pid);
    }
    double kaon_integral_after_chi2pid = cb_kaon->Integral(0.9845, 1.017) / binWidth;
    double kaon_integral_error_after_chi2pid = cb_kaon->IntegralError(0.9845, 1.017, totalCB->GetParameters(), fitResult->GetCovarianceMatrix().GetMatrixArray()) / binWidth;
    if (kaon_integral_error_after_chi2pid <= 0 || !fitResult->IsValid()) {
        cout << "Warning: Kaon integral error after (chi2pid cut) is zero or negative, using sqrt(N) approximation" << endl;
        kaon_integral_error_after_chi2pid = sqrt(kaon_integral_after_chi2pid);
    }
    double total_after_chi2pid = histo_beta_after_chi2pid->Integral();
    double total_after_error_chi2pid = sqrt(total_after_chi2pid);

    // Calculate efficiency and contamination for chi2pid cut
    double efficiency_chi2pid = pion_integral_after_chi2pid / pion_integral_before;
    double efficiency_error_chi2pid = efficiency_chi2pid * sqrt(pow(pion_integral_error_after_chi2pid / pion_integral_after_chi2pid, 2) + pow(pion_integral_error_before / pion_integral_before, 2));
    double kaon_contamination_chi2pid = kaon_integral_after_chi2pid / total_after_chi2pid;
    double kaon_contamination_error_chi2pid = kaon_contamination_chi2pid * sqrt(pow(kaon_integral_error_after_chi2pid / kaon_integral_after_chi2pid, 2) + pow(total_after_error_chi2pid / total_after_chi2pid, 2));

    // Output to file
    ofstream outFile("fit_parameters_bin4.txt");
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
        outFile << "Fit Status: " << fitResult->Status() << "\n";
        outFile << "IsValid: " << (fitResult->IsValid() ? "Yes" : "No") << "\n";
        outFile << "------------------------------------\n";
        outFile << "Pion Efficiency and Kaon Contamination Results (dt Cut)\n";
        outFile << "--------------------------------\n";
        outFile << "Momentum range: [" << pLow << ", " << pHigh << "] GeV/c\n";
        outFile << "Average momentum: " << p_mid << " GeV/c\n";
        outFile << "β range after dt cut: [" << beta_min_after_dt << ", " << beta_max_after_dt << "]\n";
        outFile << "--------------------------------\n";
        outFile << "True pions before cuts: " << pion_integral_before << " ± " << pion_integral_error_before << "\n";
        outFile << "Pions after dt cut: " << pion_integral_after_dt << " ± " << pion_integral_error_after_dt << "\n";
        outFile << "Kaons after dt cut: " << kaon_integral_after_dt << " ± " << kaon_integral_error_after_dt << "\n";
        outFile << "Total events after dt cut: " << total_after_dt << " ± " << total_after_error_dt << "\n";
        outFile << "--------------------------------\n";
        outFile << "Efficiency (dt cut): " << efficiency_dt * 100 << " ± " << efficiency_error_dt * 100 << " %\n";
        outFile << "Kaon Contamination (dt cut): " << kaon_contamination_dt * 100 << " ± " << kaon_contamination_error_dt * 100 << " %\n";
        outFile << "--------------------------------\n";
        outFile << "Pion Efficiency and Kaon Contamination Results (chi2pid Cut)\n";
        outFile << "--------------------------------\n";
        outFile << "Momentum range: [" << pLow << ", " << pHigh << "] GeV/c\n";
        outFile << "Average momentum: " << p_mid << " GeV/c\n";
        outFile << "β range after chi2pid cut: [" << beta_min_after_chi2pid << ", " << beta_max_after_chi2pid << "]\n";
        outFile << "--------------------------------\n";
        outFile << "True pions before cuts: " << pion_integral_before << " ± " << pion_integral_error_before << "\n";
        outFile << "Pions after chi2pid cut: " << pion_integral_after_chi2pid << " ± " << pion_integral_error_after_chi2pid << "\n";
        outFile << "Kaons after chi2pid cut: " << kaon_integral_after_chi2pid << " ± " << kaon_integral_error_after_chi2pid << "\n";
        outFile << "Total events after chi2pid cut: " << total_after_chi2pid << " ± " << total_after_error_chi2pid << "\n";
        outFile << "--------------------------------\n";
        outFile << "Efficiency (chi2pid cut): " << efficiency_chi2pid * 100 << " ± " << efficiency_error_chi2pid * 100 << " %\n";
        outFile << "Kaon Contamination (chi2pid cut): " << kaon_contamination_chi2pid * 100 << " ± " << kaon_contamination_error_chi2pid * 100 << " %\n";
        outFile << "--------------------------------\n";
        outFile.close();
    } else {
        cerr << "Error: Cannot open fit_parameters_bin4.txt for writing" << endl;
    }

    // Residuals plot
    TCanvas* c_res = new TCanvas("c_res", "Residuals", 800, 600);
    auto hist_res = (TH1D*)histo_beta_before->Clone("hist_res");
    for (int i = 1; i <= hist_res->GetNbinsX(); ++i) {
        double x = hist_res->GetBinCenter(i);
        double data = hist_res->GetBinContent(i);
        double fit_val = totalCB->Eval(x);
        hist_res->SetBinContent(i, data - fit_val);
    }
    hist_res->Draw();
    c_res->Print("residuals_bin4.png");

    canvas->Print("bin_4.png");

    delete leg;
    delete line_pi;
    delete line_K;
    delete line_p;
    delete totalCB;
    delete cb_pion;
    delete cb_kaon;
    delete bg;
    delete canvas;
    delete c_res;
    delete hist_res;
    file->Close();
    delete file;

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Total program time: " << duration.count() / 1000.0 << " seconds" << endl;

    return 0;
}