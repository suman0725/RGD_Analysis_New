

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

    // theoretical beta range for kaon 0.99 - 0.995 , pion 0.995 -1.003, proton 0.97 - 0.98

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
    TF1* protonFit = new TF1("protonFit", "[0]*exp(-0.5*((x-[1])/[2])^2)", 0.96, 0.98);
    protonFit->SetParameter(0, amp_p);
    protonFit->FixParameter(1, beta_p);
    protonFit->SetParameter(2, 0.007);

    histo_beta->Fit(protonFit, "RME+", "", 0.968, 0.982);
    double protonAmp = protonFit->GetParameter(0);
    double protonSigma = protonFit->GetParameter(2);
    cout << "Proton Fit: Amp = " << protonAmp << ", Sigma = " << protonSigma << endl;
    protonFit->SetRange(0.95, 1.03);
    protonFit->SetLineColor(kRed);
    protonFit->Draw("SAME");

    // Kaon fit (0.990 to 0.995)
    TF1* kaonFit = new TF1("kaonFit", "[0]*exp(-0.5*((x-[1])/[2])^2)", 0.99, beta_K);           
    kaonFit->SetParameter(0, amp_K);
    kaonFit->FixParameter(1, beta_K);
    kaonFit->SetParameter(2, 0.003);
   
    kaonFit->SetNpx(1000);
    histo_beta->Fit(kaonFit, "RME+", "", 0.990, 0.993);
    double kaonAmp = kaonFit->GetParameter(0);
    double kaonSigma = kaonFit->GetParameter(2);
    cout << "Kaon Fit: Amp = " << kaonAmp << ", Sigma = " << kaonSigma << endl;
    kaonFit->SetRange(0.95, 1.03);
    kaonFit->SetLineColor(kGreen);
    kaonFit->Draw("SAME");

    // Pion fit (0.9953 to 1.006)
    TF1* pionFit = new TF1("pionFit", "[0]*exp(-0.5*((x-[1])/[2])^2)", 0.9953, 1.006);
    pionFit->SetParameter(0, amp_pi);
    pionFit->FixParameter(1, beta_pi);
    pionFit->SetParameter(2, 0.004);
    // pionFit->SetParLimits(2, 0.0005, 0.02); 
    pionFit->SetNpx(1000);
    histo_beta->Fit(pionFit, "RMEQ+", "", 0.996, 1.004);

    double pionSigma = pionFit->GetParameter(2);
    pionFit->SetLineColor(kBlue);
    pionFit->SetRange(0.95, 1.03);
    
    pionFit->Draw("SAME");

    // Combined fit over full range
    // Combined fit over full range
   
TF1* totalFit = new TF1("totalFit", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[4])/[5])^2) + [6]*exp(-0.5*((x-[7])/[8])^2) + [9] + [10]*x + [11]*exp(-[12]*(x-[13])) + [14]*exp(-[15]*(x-[16]))", 0.95, 1.03);
const char* parNames[17] = {"Pion_Amp", "Pion_Mean", "Pion_Sigma", "Kaon_Amp", "Kaon_Mean", "Kaon_Sigma", "Proton_Amp", "Proton_Mean", "Proton_Sigma", "Bg_Const", "Bg_Linear", "Tail1_Amp", "Tail1_Decay", "Tail1_Center", "Tail2_Amp", "Tail2_Decay", "Tail2_Center"};
for (int i = 0; i < 17; i++) totalFit->SetParName(i, parNames[i]);
Double_t params[17] = {protonAmp, beta_p, protonSigma, kaonAmp, beta_K, kaonSigma, amp_pi, beta_pi, 0.007, 500.0, 0.0, 800.0, 5.0, 0.96, 600.0, 15.0, 1.02};
totalFit->SetParameters(params);
totalFit->SetParLimits(0, protonAmp * 0.8, protonAmp * 1.2);
totalFit->SetParLimits(1, beta_p - 0.005, beta_p + 0.005);
totalFit->SetParLimits(2, protonSigma * 0.8, protonSigma * 1.5);
totalFit->SetParLimits(3, kaonAmp * 0.8, kaonAmp * 1.2);
totalFit->SetParLimits(4, beta_K - 0.0003, beta_K + 0.0003);
totalFit->SetParLimits(5, kaonSigma * 1, kaonSigma * 1.5);

totalFit->SetParLimits(6, amp_pi * 0.8, amp_pi * 1.2);
totalFit->SetParLimits(7, beta_pi - 0.000001, beta_pi + 0.0001);
totalFit->SetParLimits(8, pionSigma * 0.8, pionSigma * 0.94);

totalFit->SetParLimits(9, 0.0, 1000.0);
totalFit->SetParLimits(10, -10, 10.0);

totalFit->SetParLimits(11, 0.0, 1000.0); // Tail1_Amp
totalFit->SetParLimits(12, 0.1, 10.0);   // Tail1_Decay
totalFit->SetParLimits(13, 0.95, 0.98);  // Tail1_Center

totalFit->SetParLimits(14, 0.0, 500.0); // Tail2_Amp (right tail)
totalFit->SetParLimits(15, 0.1, 2.0);    // Tail2_Decay (slower decay for right tail)
totalFit->SetParLimits(16, 1.006, 1.03);  // Tail2_Center

totalFit->SetNpx(1000);
histo_beta->Fit(totalFit, "RME+M", "", 0.95, 1.03);
totalFit->SetLineColor(kBlack);
totalFit->Draw("SAME");

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