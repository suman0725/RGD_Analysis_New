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

using namespace std;
using namespace ROOT;

// Crystal Ball function
Double_t crystalBallFunc(Double_t* x, Double_t* par) {
    Double_t t = (x[0] - par[1]) / par[2]; // (x - mean) / sigma
    Double_t alpha = par[3];               // tail parameter
    Double_t n = par[4];                   // shape parameter
    if (t >= -alpha) {
        return par[0] * exp(-0.5 * t * t); // Gaussian core
    } else {
        Double_t A = pow(n / alpha, n) * exp(-0.5 * alpha * alpha);
        Double_t B = n / alpha - alpha;
        return par[0] * A / pow(B - t, n); // Power-law tail
    }
}

// Combined function for three Crystal Ball fits
Double_t totalFitFunc(Double_t* x, Double_t* par) {
    Double_t pion = crystalBallFunc(x, &par[0]);   // Pion params: [0-4]
    Double_t kaon = crystalBallFunc(x, &par[5]);   // Kaon params: [5-9]
    Double_t proton = crystalBallFunc(x, &par[10]); // Proton params: [10-14]
    return pion + kaon + proton;
}

int main() {
    auto start = chrono::high_resolution_clock::now();
    gStyle->SetOptStat(0);

    // Step 1: Open ROOT file and create RDataFrame
    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/Skim/pkptreeCxC_9_test_modified.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file pkptreeCxC_9_test_modified.root\n";
        return 1;
    }
    ROOT::RDataFrame df_all("EB_all_pion_assumed", file);

    // Step 2: Create output directory
    gSystem->mkdir("output/fit_plots", kTRUE);

    // Step 3: Filter momentum range
    double pLow = 4.0, pHigh = 4.3;
    auto filtered_df = df_all.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});

    // Step 4: Create histogram with adjusted binning
    ROOT::RDF::TH1DModel modelBeta("beta_fit", "p: [4.0-4.3) GeV/c; #beta; Counts", 100, 0.95, 1.03);
    auto histo_beta = filtered_df.Histo1D(modelBeta, "beta");

    // Step 5: Calculate theoretical beta values
    double p_mid = (pLow + pHigh) / 2;
    double beta_pi = p_mid / sqrt(p_mid * p_mid + 0.1396 * 0.1396);
    double beta_K = p_mid / sqrt(p_mid * p_mid + 0.4937 * 0.4937);
    double beta_p = p_mid / sqrt(p_mid * p_mid + 0.9383 * 0.9383);

    // Step 6: Get amplitudes directly at theoretical beta values
    int bin_pi = histo_beta->FindBin(beta_pi);
    int bin_K = histo_beta->FindBin(beta_K);
    int bin_p = histo_beta->FindBin(beta_p);
    double amp_pi = histo_beta->GetBinContent(bin_pi);
    double amp_K = histo_beta->GetBinContent(bin_K);
    double amp_p = histo_beta->GetBinContent(bin_p);
    cout << "Pion amplitude at beta = " << beta_pi << ": " << amp_pi << endl;
    cout << "Kaon amplitude at beta = " << beta_K << ": " << amp_K << endl;
    cout << "Proton amplitude at beta = " << beta_p << ": " << amp_p << endl;

    // Step 7: Initialize canvas
    TCanvas* canvas = new TCanvas("canvas", "Beta Fit with Signals", 1400, 1000);
    canvas->SetGridx(1);
    canvas->SetGridy(1);
    histo_beta->SetMarkerStyle(20);
    histo_beta->SetMarkerSize(1.2);
    histo_beta->SetMarkerColor(kBlack);
    histo_beta->Draw("P");

    // Proton fit (0.96 to 0.98)
    TF1* protonFit = new TF1("protonFit", crystalBallFunc, 0.96, 0.98, 5);
    protonFit->SetParameter(0, amp_p);
    protonFit->FixParameter(1, beta_p);
    protonFit->SetParameter(2, 0.007);
    protonFit->SetParameter(3, 1.0); // alpha
    protonFit->SetParameter(4, 2.0); // n
    histo_beta->Fit(protonFit, "RME+", "", 0.96, 0.98);
    protonFit->SetRange(0.95, 1.03);
    protonFit->SetLineColor(kRed);
    protonFit->Draw("SAME");
    // Print proton fit parameters
    cout << "Proton Fit Parameters:" << endl;
    for (int i = 0; i < 5; i++) {
        cout << "p" << i << " = " << protonFit->GetParameter(i) << " +/- " << protonFit->GetParError(i) << endl;
    }

    // Kaon fit (0.990 to 0.993)
    TF1* kaonFit = new TF1("kaonFit", crystalBallFunc, 0.990, 0.993, 5);
    kaonFit->SetParameter(0, amp_K);
    kaonFit->FixParameter(1, beta_K);
    kaonFit->SetParameter(2, 0.003);
    kaonFit->SetParameter(3, 1.0); // alpha
    kaonFit->SetParameter(4, 2.0); // n
    kaonFit->SetNpx(1000);
    histo_beta->Fit(kaonFit, "RME+", "", 0.990, 0.993);
    kaonFit->SetRange(0.95, 1.03);
    kaonFit->SetLineColor(kGreen);
    kaonFit->Draw("SAME");
    // Print kaon fit parameters
    cout << "Kaon Fit Parameters:" << endl;
    for (int i = 0; i < 5; i++) {
        cout << "p" << i << " = " << kaonFit->GetParameter(i) << " +/- " << kaonFit->GetParError(i) << endl;
    }

    // Pion fit (0.998 to 1.004)
    TF1* pionFit = new TF1("pionFit", crystalBallFunc, 0.998, 1.01, 5);
    pionFit->SetParameter(0, amp_pi);
    pionFit->FixParameter(1, beta_pi);
    pionFit->SetParameter(2, 0.004);
    pionFit->SetParameter(3, 1.0); // alpha
    pionFit->SetParameter(4, 2.0); // n
    pionFit->SetNpx(1000);
    histo_beta->Fit(pionFit, "RMEQ+", "", 0.998, 1.01);
    pionFit->SetRange(0.95, 1.03);
    pionFit->SetLineColor(kBlue);
    pionFit->Draw("SAME");
    // Print pion fit parameters
    cout << "Pion Fit Parameters:" << endl;
    for (int i = 0; i < 5; i++) {
        cout << "p" << i << " = " << pionFit->GetParameter(i) << " +/- " << pionFit->GetParError(i) << endl;
    }

    // Combined fit over full range
    TF1* totalFit = new TF1("totalFit", totalFitFunc, 0.95, 1.03, 15);
    totalFit->SetNpx(1000);
    // Initialize with individual fit parameters
    totalFit->FixParameter(0, pionFit->GetParameter(0));   // Pion amp
    totalFit->FixParameter(1, pionFit->GetParameter(1));   // Pion mean
    totalFit->SetParameter(2, pionFit->GetParameter(2));   // Pion sigma
    totalFit->SetParameter(3, pionFit->GetParameter(3));   // Pion alpha
    totalFit->SetParameter(4, pionFit->GetParameter(4));   // Pion n
    totalFit->SetParameter(5, kaonFit->GetParameter(0));   // Kaon amp
    totalFit->FixParameter(6, kaonFit->GetParameter(1));   // Kaon mean
    totalFit->SetParameter(7, kaonFit->GetParameter(2));   // Kaon sigma
    totalFit->SetParameter(8, kaonFit->GetParameter(3));   // Kaon alpha
    totalFit->SetParameter(9, kaonFit->GetParameter(4));   // Kaon n
    totalFit->SetParameter(10, protonFit->GetParameter(0)); // Proton amp
    totalFit->FixParameter(11, protonFit->GetParameter(1)); // Proton mean
    totalFit->SetParameter(12, protonFit->GetParameter(2)); // Proton sigma
    totalFit->SetParameter(13, protonFit->GetParameter(3)); // Proton alpha
    totalFit->SetParameter(14, protonFit->GetParameter(4)); // Proton n
    // Set wider limits to avoid NaN
    totalFit->SetParLimits(0, 0, amp_pi * 2.0);    // Pion amp
    totalFit->SetParLimits(2, 0.001, 0.008);       // Pion sigma
    totalFit->SetParLimits(3, 0.1, 2.0);           // Pion alpha
    totalFit->SetParLimits(4, 1.0, 5.0);           // Pion n
    totalFit->SetParLimits(5, 0, amp_K * 2.0);     // Kaon amp
    totalFit->SetParLimits(7, 0.001, 0.006);       // Kaon sigma
    totalFit->SetParLimits(8, 0.1, 2.0);           // Kaon alpha
    totalFit->SetParLimits(9, 1.0, 5.0);           // Kaon n
    totalFit->SetParLimits(10, 0, amp_p * 2.0);    // Proton amp
    totalFit->SetParLimits(12, 0.003, 0.012);      // Proton sigma
    totalFit->SetParLimits(13, 0.1, 2.0);          // Proton alpha
    totalFit->SetParLimits(14, 1.0, 5.0);          // Proton n
    histo_beta->Fit(totalFit, "RME+", "", 0.95, 1.03); // Robust fit
    totalFit->SetLineColor(kBlack);
    totalFit->Draw("SAME");
    // Print total fit parameters
    cout << "Total Fit Parameters:" << endl;
    for (int i = 0; i < 15; i++) {
        cout << "p" << i << " = " << totalFit->GetParameter(i) << " +/- " << totalFit->GetParError(i) << endl;
    }

    // Step 10: Draw theoretical beta lines
    TLine* line_pi = new TLine(beta_pi, 0, beta_pi, histo_beta->GetMaximum());
    line_pi->SetLineColor(kBlue);
    line_pi->SetLineStyle(2);
    line_pi->Draw();

    TLine* line_K = new TLine(beta_K, 0, beta_K, histo_beta->GetMaximum());
    line_K->SetLineColor(kGreen);
    line_K->SetLineStyle(2);
    line_K->Draw();

    TLine* line_p = new TLine(beta_p, 0, beta_p, histo_beta->GetMaximum());
    line_p->SetLineColor(kRed);
    line_p->SetLineStyle(2);
    line_p->Draw();

    // Step 11: Add legend
    TLegend* leg = new TLegend(0.65, 0.55, 0.85, 0.75);
    leg->AddEntry(histo_beta.GetPtr(), Form("Data (%.0f)", histo_beta->GetEntries()), "p");
    leg->AddEntry(totalFit, "Total Fit", "l");
    leg->AddEntry(pionFit, "Pion Signal", "l");
    leg->AddEntry(kaonFit, "Kaon Signal", "l");
    leg->AddEntry(protonFit, "Proton Signal", "l");
    leg->AddEntry(line_pi, "Theoretical #beta", "l");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    leg->Draw();

    // Step 12: Save output
    canvas->Print("output/fit_plots/beta_fit_4.0-4.3.pdf");
    canvas->Print("output/fit_plots/beta_fit_4.0-4.3.png");

    // Step 13: Clean up
    delete leg;
    delete line_pi;
    delete line_K;
    delete line_p;
    delete pionFit;
    delete kaonFit;
    delete protonFit;
    delete totalFit;
    delete canvas;
    file->Close();
    delete file;

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Total program time: " << duration.count() / 1000.0 << " seconds" << endl;

    return 0;
}