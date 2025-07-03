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

// Total function (π+ + kaons, no background)
Double_t totalCrystalBallFunc(Double_t* x, Double_t* par) {
    Double_t cb1 = doubleSidedCrystalBallFunc(x, &par[0]); // Pion
    Double_t cb2 = crystalBallFunc(x, &par[7]);            // Kaon
    return cb1 + cb2;
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

    double pLow = 1.3, pHigh = 1.6;
    auto df_pions_before = df_all.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});

    auto df_after_cut = df_pions_before.Filter([](float p, float dt) {
        return dt > getdtCutNeg(p) && dt < getdtCutPos(p);
    }, {"p", "dt"});

    ROOT::RDF::TH1DModel modelBeta("beta_fit_cb", "p: [1.3-1.6) GeV/c; #beta; Counts", 200, 0.95, 1.03);
    auto histo_beta_before = df_pions_before.Histo1D(modelBeta, "beta");
    auto histo_beta_after = df_after_cut.Histo1D(modelBeta, "beta");

    // Check histogram content
    if (histo_beta_before->GetEntries() == 0) {
        cerr << "Error: Histogram histo_beta_before is empty!" << endl;
        file->Close();
        delete file;
        return 1;
    }

    // Determine the effective β range from the blue histogram
    double beta_min_after = histo_beta_after->GetBinCenter(histo_beta_after->FindFirstBinAbove(0.1 * histo_beta_after->GetMaximum()));
    double beta_max_after = histo_beta_after->GetBinCenter(histo_beta_after->FindLastBinAbove(0.1 * histo_beta_after->GetMaximum()));

    // Calculate the true average momentum in the bin
    auto p_mean_action = df_pions_before.Mean("p");
    double p_mid = *p_mean_action;
    double beta_pi_theory = p_mid / sqrt(p_mid * p_mid + 0.1396 * 0.1396);
    double beta_K_theory = p_mid / sqrt(p_mid * p_mid + 0.4937 * 0.4937);
    double beta_p_theory = p_mid / sqrt(p_mid * p_mid + 0.95383 * 0.95383);

    double binWidth = (1.03 - 0.95) / 200;
    int binPi = static_cast<int>((beta_pi_theory - 0.95) / binWidth) + 1;
    int binK = static_cast<int>((beta_K_theory - 0.95) / binWidth) + 1;

    double amp_pi = histo_beta_before.GetPtr()->GetBinContent(binPi) * 1.1;
    double amp_K = histo_beta_before.GetPtr()->GetBinContent(binK) * 1.1;

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

    // Fit setup
    TF1* totalCB = new TF1("totalCB", totalCrystalBallFunc, 0.95, 1.03, 12); // Reduced to 12 parameters
    totalCB->SetNpx(1000);
    totalCB->SetRange(0.95, 1.03);
    Double_t initParams[12] = {
        amp_pi, beta_pi_theory, 0.004, 1.5, 2.5, 1.5, 2.5, // Pion
        amp_K, beta_K_theory, 0.005, 1.0, 2.0              // Kaon
    };
    totalCB->SetParameters(initParams);

    for (int i = 0; i < 12; ++i)
        totalCB->SetParName(i, Form("p%d", i));

    // Pion parameters
    totalCB->SetParLimits(0, 0, 2 * amp_pi);
    totalCB->SetParameter(1, beta_pi_theory);
    totalCB->SetParLimits(1, beta_pi_theory - 0.002, beta_pi_theory + 0.002);
    totalCB->SetParLimits(2, 0.002, 0.006);
    totalCB->SetParLimits(3, 0.1, 50.0);
    totalCB->SetParLimits(4, 1.0, 20.0);
    totalCB->FixParameter(5, totalCB->GetParameter(3));
    totalCB->FixParameter(6, totalCB->GetParameter(4));

    // Kaon parameters
    totalCB->SetParLimits(7, 0, 2 * amp_K);
    totalCB->SetParameter(8, beta_K_theory);
    totalCB->SetParLimits(8, beta_K_theory - 0.01, beta_K_theory + 0.01);
    totalCB->SetParLimits(9, 0.002, 0.008);
    totalCB->SetParLimits(10, 0.1, 50.0);
    totalCB->SetParLimits(11, 1.0, 200.0);

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
    // line_K->Draw();

    TLine* line_p = new TLine(beta_p_theory, 0, beta_p_theory, histo_beta_before->GetMaximum());
    line_p->SetLineColor(kRed);
    line_p->SetLineStyle(2);
    line_p->SetLineWidth(3);
    // line_p->Draw();

    TLegend* leg = new TLegend(0.65, 0.6, 0.95, 0.95);
    leg->SetTextSize(0.02);
    leg->AddEntry(histo_beta_before.GetPtr(), "+ve particles (Before 3#sigma dt Cut)", "f");
    leg->AddEntry(histo_beta_after.GetPtr(), "+ve pions (After 3#sigma dt Cut)", "f");
    leg->AddEntry(totalCB, "Total Fit", "l");
    leg->AddEntry(line_pi, "#beta_{#pi_theory}", "l");

    // Pion peak contribution
    TF1* cb_pion = new TF1("cb_pion", totalCrystalBallFunc, 0.95, 1.03, 12);
    cb_pion->SetNpx(1000);
    Double_t pionParams[12] = {0};
    for (int i = 0; i < 7; ++i) pionParams[i] = totalCB->GetParameter(i);
    pionParams[5] = pionParams[3];
    pionParams[6] = pionParams[4];
    for (int i = 7; i < 12; ++i) pionParams[i] = 0;
    cb_pion->SetParameters(pionParams);
    for (int i = 0; i < 12; ++i) cb_pion->FixParameter(i, pionParams[i]);
    cb_pion->SetLineColor(kBlue+2);
    cb_pion->SetLineWidth(2);
    cb_pion->Draw("SAME");

    // Kaon peak contribution
    TF1* cb_kaon = new TF1("cb_kaon", totalCrystalBallFunc, 0.95, 1.03, 12);
    cb_kaon->SetNpx(1000);
    Double_t kaonParams[12] = {0};
    for (int i = 7; i < 12; ++i) kaonParams[i] = totalCB->GetParameter(i);
    for (int i = 0; i < 7; ++i) kaonParams[i] = 0;
    cb_kaon->SetParameters(kaonParams);
    for (int i = 0; i < 12; ++i) cb_kaon->FixParameter(i, kaonParams[i]);
    cb_kaon->SetLineColor(kGreen);
    cb_kaon->SetLineWidth(2);
    cb_kaon->Draw("SAME");

    TCanvas* c_res = new TCanvas("c_res", "Residuals", 800, 600);
    auto hist_res = (TH1D*)histo_beta_before->Clone("hist_res");
    for (int i = 1; i <= hist_res->GetNbinsX(); ++i) {
        double x = hist_res->GetBinCenter(i);
        double data = hist_res->GetBinContent(i);
        double fit_val = totalCB->Eval(x);
        hist_res->SetBinContent(i, data - fit_val);
    }
    hist_res->Draw();
    c_res->Print("residuals_before.png");

    // Calculate integrals and errors
    double pion_integral_before = cb_pion->Integral(0.95, 1.03) / binWidth;
    double pion_integral_error_before = cb_pion->IntegralError(0.95, 1.03, totalCB->GetParameters(), fitResult->GetCovarianceMatrix().GetMatrixArray()) / binWidth;
    if (pion_integral_error_before <= 0 || !fitResult->IsValid()) {
        cout << "Warning: Pion integral error before is zero or negative, using sqrt(N) approximation" << endl;
        pion_integral_error_before = sqrt(pion_integral_before);
    }
    double pion_integral_after = cb_pion->Integral(beta_min_after, beta_max_after) / binWidth;
    double pion_integral_error_after = cb_pion->IntegralError(beta_min_after, beta_max_after, totalCB->GetParameters(), fitResult->GetCovarianceMatrix().GetMatrixArray()) / binWidth;
    if (pion_integral_error_after <= 0 || !fitResult->IsValid()) {
        cout << "Warning: Pion integral error after is zero or negative, using sqrt(N) approximation" << endl;
        pion_integral_error_after = sqrt(pion_integral_after);
    }
    double kaon_integral_after = cb_kaon->Integral(beta_min_after, beta_max_after) / binWidth;
    double kaon_integral_error_after = cb_kaon->IntegralError(beta_min_after, beta_max_after, totalCB->GetParameters(), fitResult->GetCovarianceMatrix().GetMatrixArray()) / binWidth;
    if (kaon_integral_error_after <= 0 || !fitResult->IsValid()) {
        cout << "Warning: Kaon integral error after is zero or negative, using sqrt(N) approximation" << endl;
        kaon_integral_error_after = sqrt(kaon_integral_after);
    }
    double total_after = histo_beta_after->Integral();
    double total_after_error = sqrt(total_after);

    // Calculate efficiency
    double efficiency = pion_integral_after / pion_integral_before;
    double efficiency_error = efficiency * sqrt(pow(pion_integral_error_after / pion_integral_after, 2) + pow(pion_integral_error_before / pion_integral_before, 2));

    // Calculate kaon contamination
    double kaon_contamination = kaon_integral_after / total_after;
    double kaon_contamination_error = kaon_contamination * sqrt(pow(kaon_integral_error_after / kaon_integral_after, 2) + pow(total_after_error / total_after, 2));

    // Output to file
    ofstream outFile("fit_parameters_bin1.txt");
    if (outFile.is_open()) {
        outFile << "Fit Parameters for Crystal Ball Fit (Bin: p = 1.3 - 1.6 GeV/c)\n";
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
        outFile << "Pion Efficiency and Kaon Contamination Results\n";
        outFile << "--------------------------------\n";
        outFile << "Momentum range: [" << pLow << ", " << pHigh << "] GeV/c\n";
        outFile << "β range after cut: [" << beta_min_after << ", " << beta_max_after << "]\n";
        outFile << "--------------------------------\n";
        outFile << "True pions before cuts: " << pion_integral_before << " ± " << pion_integral_error_before << "\n";
        outFile << "Pions after cuts: " << pion_integral_after << " ± " << pion_integral_error_after << "\n";
        outFile << "Kaons after cuts: " << kaon_integral_after << " ± " << kaon_integral_error_after << "\n";
        outFile << "Total events after cuts: " << total_after << " ± " << total_after_error << "\n";
        outFile << "--------------------------------\n";
        outFile << "Efficiency: " << efficiency * 100 << " ± " << efficiency_error * 100 << " %\n";
        outFile << "Kaon Contamination: " << kaon_contamination * 100 << " ± " << kaon_contamination_error * 100 << " %\n";
        outFile << "--------------------------------\n";
        outFile.close();
    } else {
        cerr << "Error: Cannot open fit_parameters_bin1.txt for writing" << endl;
    }

    leg->AddEntry(cb_pion, "Pion fit", "l");
    leg->AddEntry(cb_kaon, "Kaon fit", "l");
    leg->Draw();

    canvas->Print("bin_1.png");

    delete leg;
    delete line_pi;
    delete line_K;
    delete line_p;
    delete totalCB;
    delete cb_pion;
    delete cb_kaon;
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