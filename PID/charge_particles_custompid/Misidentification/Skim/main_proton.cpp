#include "include/ContaminationUtils.h"
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLine.h>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

int main() {
    // Step 1: Open the ROOT file and load the trees for pions and protons
    TFile* file = new TFile("../Skim/pkptreeCxC_9_test_modified.root");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file ../Skim/pkptreeCxC_7.root\n";
        return 1;
    }
    TTree* treePions = (TTree*)file->Get("EB_pid_pions");
    TTree* treeProtons = (TTree*)file->Get("EB_pid_protons");
    if (!treePions || !treeProtons) {
        cerr << "Error: Cannot load trees\n";
        file->Close();
        delete file;
        return 1;
    }

    // Define maximum entries to process for testing
    Long64_t maxEntries = 1000000;

    // Step 2: Create output directories for saving PDFs and CSVs
    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/pdf/proton", kTRUE);
    gSystem->mkdir("output/csv/proton", kTRUE);

    // Step 3: Set up momentum bins (1 to 7 GeV, 20 bins)
    vector<double> pBins;
    double pMin = 1.0, pMax = 7.0;
    int nBins = 20;
    double binWidth = (pMax - pMin) / nBins;
    for (int i = 0; i <= nBins; i++) pBins.push_back(pMin + i * binWidth);

    // Step 4: Set chi2pid histogram ranges, fit ranges, and fit parameters for each bin
    vector<pair<double, double>> chi2HistRanges(nBins);
    vector<pair<double, double>> chi2FitRanges(nBins);
    vector<tuple<double, double, double>> chi2FitParams(nBins);
    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i];
        if (pLow < 2.0) {
            //chi2HistRanges[i] = {-5.0, 5.0};

            // adjusted range for chi2pid of protons and pions histogram range
            chi2HistRanges[i] = {-10.0, 10.0};
            chi2FitRanges[i] = {-1.3, 1.3};
            chi2FitParams[i] = {1000.0, 0.0, 1.0};
        } else if (pLow < 3.0) {
            //chi2HistRanges[i] = {-4.0, 4.0};
            chi2HistRanges[i] = {-10.0, 10.0};
            chi2FitRanges[i] = {-1.5, 1.5};
            chi2FitParams[i] = {800.0, 0.1, 0.8};
        } else {
            //chi2HistRanges[i] = {-3.0, 3.0};
            chi2HistRanges[i] = {-10.0, 10.0};
            chi2FitRanges[i] = {-1, 1};
            chi2FitParams[i] = {500.0, -0.1, 0.5};
        }
    }

    // Step 5: Set beta histogram ranges and fit ranges for each bin
    vector<pair<double, double>> betaHistRanges(nBins);
    vector<pair<double, double>> betaFitRanges(nBins);
    for (int i = 0; i < nBins; i++) {
        betaHistRanges[i] = {0.95, 1.03}; // Adjusted lower range for protons
        betaFitRanges[i] = {0.95, 1.02};
    }

    // Step 6: Set beta cuts for pions and protons (same for before and after chi2pid cuts)
    vector<double> pionLeftBefore(nBins), pionRightBefore(nBins), protonRightBefore(nBins);
    vector<double> pionLeftAfter(nBins), pionRightAfter(nBins), protonRightAfter(nBins);
    for (int i = 0; i < nBins; i++) {
        
        if (i == 0) { // Bin [1.0-1.3)
            pionLeftBefore[i] = 0.97; pionRightBefore[i] = 1.02; protonRightBefore[i] = 0.85;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  protonRightAfter[i] = 0.85;
        } else if (i == 1) { // Bin [1.3-1.6)
            pionLeftBefore[i] = 0.9755; pionRightBefore[i] = 1.0135; protonRightBefore[i] = 0.88;
            pionLeftAfter[i] = 0.9755;  pionRightAfter[i] = 1.0135;  protonRightAfter[i] = 0.88;
        } else if (i == 2) { // Bin [1.6-1.9)
            pionLeftBefore[i] = 0.98;   pionRightBefore[i] = 1.0140; protonRightBefore[i] = 0.90;
            pionLeftAfter[i] = 0.98;    pionRightAfter[i] = 1.0140;  protonRightAfter[i] = 0.90;
        } else if (i == 3) { // Bin [1.9-2.2)
            pionLeftBefore[i] = 0.9815; pionRightBefore[i] = 1.0145; protonRightBefore[i] = 0.92;
            pionLeftAfter[i] = 0.9815;  pionRightAfter[i] = 1.0145;  protonRightAfter[i] = 0.92;
        } else if (i == 4) { // Bin [2.2-2.5)
            pionLeftBefore[i] = 0.9880; pionRightBefore[i] = 1.0141; protonRightBefore[i] = 0.94;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0141;  protonRightAfter[i] = 0.94;
        } else if (i == 5) { // Bin [2.5-2.8)
            pionLeftBefore[i] = 0.9848; pionRightBefore[i] = 1.013;  protonRightBefore[i] = 0.95;
            pionLeftAfter[i] = 0.9848;  pionRightAfter[i] = 1.013;   protonRightAfter[i] = 0.95;
        } else if (i == 6) { // Bin [2.8-3.1)
            pionLeftBefore[i] = 0.98;  pionRightBefore[i] = 1.03; protonRightBefore[i] = 1.02;
            pionLeftAfter[i] = 0.986;   pionRightAfter[i] = 1.0128;  protonRightAfter[i] = 0.96;
        } else if (i == 7) { // Bin [3.1-3.4)
            pionLeftBefore[i] = 0.98; pionRightBefore[i] = 1.03; protonRightBefore[i] = 1.02;
            pionLeftAfter[i] = 0.9868;  pionRightAfter[i] = 1.0121;  protonRightAfter[i] = 0.97;
        } else if (i == 8) { // Bin [3.4-3.7)
            pionLeftBefore[i] = 0.9879; pionRightBefore[i] = 1.0120; protonRightBefore[i] = 0.975;
            pionLeftAfter[i] = 0.9879;  pionRightAfter[i] = 1.0120;  protonRightAfter[i] = 0.975;
        } else if (i == 9) { // Bin [3.7-4.0)
            pionLeftBefore[i] = 0.9880; pionRightBefore[i] = 1.0120; protonRightBefore[i] = 0.98;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  protonRightAfter[i] = 0.98;
        } else if (i == 10) { // Bin [3.7-4.0)
            pionLeftBefore[i] = 0.98; pionRightBefore[i] = 1.03; protonRightBefore[i] = 1.02;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  protonRightAfter[i] = 0.98;
        } else if (i == 11) { // Bin [3.7-4.0)
            pionLeftBefore[i] = 0.98; pionRightBefore[i] = 1.103; protonRightBefore[i] = 1.02;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  protonRightAfter[i] = 0.98;
        } else if (i == 12) { // Bin [3.7-4.0)
            pionLeftBefore[i] = 0.98; pionRightBefore[i] = 1.03; protonRightBefore[i] = 1.02;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  protonRightAfter[i] = 0.98;
        } else if (i == 13) { // Bin [3.7-4.0)
            pionLeftBefore[i] = 0.98; pionRightBefore[i] = 1.03; protonRightBefore[i] = 1.02;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  protonRightAfter[i] = 0.98;
        } else if (i == 14) { // Bin [3.7-4.0)
            pionLeftBefore[i] = 0.98; pionRightBefore[i] = 1.03; protonRightBefore[i] = 1.02;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  protonRightAfter[i] = 0.98;
        } else if (i == 15) { // Bin [3.7-4.0)
            pionLeftBefore[i] = 0.98; pionRightBefore[i] = 1.03; protonRightBefore[i] = 1.02;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  protonRightAfter[i] = 0.98;
        } else if (i == 16) { // Bin [3.7-4.0)
            pionLeftBefore[i] = 0.98; pionRightBefore[i] = 1.03; protonRightBefore[i] = 1.02;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  protonRightAfter[i] = 0.98;
        } else if (i == 17) { // Bin [3.7-4.0)
            pionLeftBefore[i] = 0.98; pionRightBefore[i] = 1.03; protonRightBefore[i] = 1.02;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  protonRightAfter[i] = 0.98;
        } else if (i == 18) { // Bin [3.7-4.0)
            pionLeftBefore[i] = 0.98; pionRightBefore[i] = 1.03; protonRightBefore[i] = 1.02;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  protonRightAfter[i] = 0.98;
        } else if (i == 19) { // Bin [3.7-4.0)
            pionLeftBefore[i] = 0.98; pionRightBefore[i] = 1.03; protonRightBefore[i] = 1.02;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  protonRightAfter[i] = 0.98;
        }
         else {
            cerr << "Warning: Using default beta cuts for bin [" << pBins[i] << "-" << pBins[i+1] << ")\n";
            pionLeftBefore[i] = 0.9880; pionRightBefore[i] = 1.0120; protonRightBefore[i] = 0.98;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  protonRightAfter[i] = 0.98;
        }
    }

    // Step 7: Create canvases and open PDF files for plotting
    TCanvas* canvasChi2Before = new TCanvas("canvasChi2Before", "Chi2pid Fits (Before)", 1900, 900); //width 1900 was 1200 before 
    canvasChi2Before->SetMargin(0.1, 0.05, 0.15, 0.1);
    TCanvas* canvasChi2After = new TCanvas("canvasChi2After", "Chi2pid Fits (After)", 1900, 900);
    canvasChi2After->SetMargin(0.1, 0.05, 0.15, 0.1);
    TCanvas* canvasBetaBefore = new TCanvas("canvasBetaBefore", "Beta Plots (Before)", 1200, 800);
    canvasBetaBefore->SetMargin(0.1, 0.05, 0.15, 0.1);
    TCanvas* canvasBetaAfter = new TCanvas("canvasBetaAfter", "Beta Plots (After)", 1200, 800);
    canvasBetaAfter->SetMargin(0.1, 0.05, 0.15, 0.1);
    canvasChi2Before->Print("output/pdf/proton/chi2pid_fits_linear_before.pdf[");
    canvasChi2After->Print("output/pdf/proton/chi2pid_fits_linear_after.pdf[");
    canvasBetaBefore->Print("output/pdf/proton/beta_before_cut_linear.pdf[");
    canvasBetaAfter->Print("output/pdf/proton/beta_after_cut_linear.pdf[");

    // Step 8: Open CSV files to save results
    // - proton_to_pion_ratio_before_chi2pid.csv: Stores the proton-to-pion ratio before the chi2pid cut
    // - contamination_after_chi2pid.csv: Stores the contamination (proton-to-pion ratio) after the chi2pid cut
    // - chi2pid_cuts.csv: Stores the pion mean, sigma, and mean ± 3σ ranges from the "before" cut
    ofstream csvFileBefore("output/csv/proton/proton_to_pion_ratio_before_chi2pid.csv");
    csvFileBefore << "Momentum Bin (GeV/c),Proton Left,Contamination (%)\n";
    ofstream csvFileAfter("output/csv/proton/contamination_after_chi2pid.csv");
    csvFileAfter << "Momentum Bin (GeV/c),Proton Left,Proton-to-Pion Ratio (%)\n";
    ofstream csvFileCuts("output/csv/proton/chi2pid_cuts.csv");
    csvFileCuts << "Momentum Bin (GeV/c),Pion Mean,Pion Sigma,Chi2 Min (Mean-3σ),Chi2 Max (Mean+3σ)\n";

    // Step 9: Loop over each momentum bin to process data
    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];

        // --- Before Chi2pid Cut ---
        // Create chi2pid histograms for pions and protons (before cut)
        TH1F* chi2PionsBefore = new TH1F(TString::Format("chi2_pions_before_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                         100, -10, 10);
        TH1F* chi2ProtonsBefore = new TH1F(TString::Format("chi2_protons_before_%d", i),
                                          TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                          100, -10, 10);

        // Fill chi2pid histograms (before cut) with momentum cut
        float p, orig_chi2pid, beta, recomputed_chi2pid;
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("orig_chi2pid", &orig_chi2pid);
        /* for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh) chi2PionsBefore->Fill(orig_chi2pid);
        } */
        //for test
        Long64_t nEntriesPions = min(treePions->GetEntries(), maxEntries);
        for (Long64_t j = 0; j < nEntriesPions; j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh) chi2PionsBefore->Fill(orig_chi2pid);
        }

        treeProtons->SetBranchAddress("p", &p);

        // branch for chi2pid based on the proton hypothesis from EB for protons (2212) 
        //treeProtons->SetBranchAddress("orig_chi2pid", &orig_chi2pid);

        // branch for chi2pid based on the pion hypothesis from EB for protons (2212) 
        treeProtons->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);

        /* for (Long64_t j = 0; j < treeProtons->GetEntries(); j++) {
            treeProtons->GetEntry(j);

            //Filling chi2pid based on the proton hypothesis from EB for protons (2212)
            //if (p >= pLow && p < pHigh) chi2ProtonsBefore->Fill(orig_chi2pid);

            //Filling chi2pid based on the pion hypothesis from EB for protons (2212)
            if (p >= pLow && p < pHigh) chi2ProtonsBefore->Fill(recomputed_chi2pid);
        } */

       //for test
       // Limit the number of entries processed
        Long64_t nEntriesProtons = min(treeProtons->GetEntries(), maxEntries);
        for (Long64_t j = 0; j < nEntriesProtons; j++) {
            treeProtons->GetEntry(j);
            if (p >= pLow && p < pHigh) chi2ProtonsBefore->Fill(recomputed_chi2pid);
        }

        // Set colors for chi2pid histograms
        chi2PionsBefore->SetLineColor(kBlue);
        chi2ProtonsBefore->SetLineColor(kRed); // Changed to red for protons

        // Fit chi2pid histograms with a Gaussian to determine the ±3σ cut
        double pionMeanBefore, pionSigmaBefore, pionChi2Min, pionChi2Max;
        auto [pionAmp, pionMeanInit, pionSigmaInit] = chi2FitParams[i];
        auto [protonAmp, protonMean, protonSigma] = chi2FitParams[i];
        TF1* fitPion = new TF1("gausFit", "gaus", chi2FitRanges[i].first, chi2FitRanges[i].second);
        fitPion->SetParameters(pionAmp, pionMeanInit, pionSigmaInit);
        chi2PionsBefore->Fit(fitPion, "R", "", chi2FitRanges[i].first, chi2FitRanges[i].second);
        pionMeanBefore = fitPion->GetParameter(1);
        pionSigmaBefore = fitPion->GetParameter(2);
        pionChi2Min = pionMeanBefore - 5 * pionSigmaBefore;
        pionChi2Max = pionMeanBefore + 5 * pionSigmaBefore;
        //pionChi2Min = pionMeanBefore - 3 * pionSigmaBefore;
        //pionChi2Max = pionMeanBefore + 3 * pionSigmaBefore; 
        TF1* fitProton = new TF1("gausFit", "gaus", chi2FitRanges[i].first, chi2FitRanges[i].second);
        fitProton->SetParameters(protonAmp, protonMean, protonSigma);
        chi2ProtonsBefore->Fit(fitProton, "R", "", chi2FitRanges[i].first, chi2FitRanges[i].second);
        cout << "Bin [" << pLow << "-" << pHigh << "): Pions chi2 range (before): [" 
             << pionChi2Min << ", " << pionChi2Max << "]\n";

        // Save the pion chi2pid fit parameters (mean, sigma, ±3σ range) to chi2pid_cuts.csv
        csvFileCuts << pLow << "-" << pHigh << "," << pionMeanBefore << "," << pionSigmaBefore << "," 
                    << pionChi2Min << "," << pionChi2Max << "\n";

        // Plot chi2pid histograms (before cut)
        canvasChi2Before->Clear();
        canvasChi2Before->Divide(2, 1);
        canvasChi2Before->cd(1);
        chi2PionsBefore->Draw("HIST");
        if (chi2PionsBefore->GetFunction("gausFit")) {
            chi2PionsBefore->GetFunction("gausFit")->SetLineColor(kBlue);
            chi2PionsBefore->GetFunction("gausFit")->Draw("SAME");
        }
        TLegend* legPions = new TLegend(0.1, 0.75, 0.3, 0.9);
        legPions->AddEntry(chi2PionsBefore, "Pions", "l");
        legPions->SetBorderSize(0);
        legPions->SetFillStyle(0);
        legPions->SetTextSize(0.025);
        legPions->Draw();
        canvasChi2Before->cd(2);
        chi2ProtonsBefore->Draw("HIST");
        if (chi2ProtonsBefore->GetFunction("gausFit")) {
            chi2ProtonsBefore->GetFunction("gausFit")->SetLineColor(kRed);
            chi2ProtonsBefore->GetFunction("gausFit")->Draw("SAME");
        }
        TLegend* legProtons = new TLegend(0.1, 0.75, 0.3, 0.9);
        legProtons->AddEntry(chi2ProtonsBefore, "Protons", "l");
        legProtons->SetBorderSize(0);
        legProtons->SetFillStyle(0);
        legProtons->SetTextSize(0.025);
        legProtons->Draw();
        canvasChi2Before->Update();
        canvasChi2Before->Print("output/pdf/proton/chi2pid_fits_linear_before.pdf");
        delete legPions;
        delete legProtons;

        // Create beta histograms (before cut)
        TH1F* betaPionsBefore = new TH1F(TString::Format("beta_pions_before_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                         70, betaHistRanges[i].first, betaHistRanges[i].second);
        TH1F* betaProtonsBefore = new TH1F(TString::Format("beta_protons_before_%d", i),
                                          TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                          70, betaHistRanges[i].first, betaHistRanges[i].second);
        betaPionsBefore->Sumw2();
        betaProtonsBefore->Sumw2();

        // Fill beta histograms (before cut) with momentum cut
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("beta", &beta);
        /* for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh) {
                betaPionsBefore->Fill(beta);
            }
        } */

       //for test
       for (Long64_t j = 0; j < nEntriesPions; j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh) {
                betaPionsBefore->Fill(beta);
            }
        }
        treeProtons->SetBranchAddress("p", &p);
        treeProtons->SetBranchAddress("beta", &beta);
        /* for (Long64_t j = 0; j < treeProtons->GetEntries(); j++) {
            treeProtons->GetEntry(j);
            if (p >= pLow && p < pHigh) {
                betaProtonsBefore->Fill(beta);
            }
        }
        */
        for (Long64_t j = 0; j < nEntriesProtons; j++) {
            treeProtons->GetEntry(j);
            if (p >= pLow && p < pHigh) {
                betaProtonsBefore->Fill(beta);
            }
        }

        // Fit beta histograms (before cut) to calculate the proton-to-pion ratio
        double pMeanBefore = 0, pSigmaBefore = 0, pConstBefore = 0, protonMeanBefore = 0, protonSigmaBefore = 0, protonConstBefore = 0;
        if (betaPionsBefore->GetEntries() >= 10) {
            TF1* gausPionBefore = new TF1(TString::Format("gaus_%s", betaPionsBefore->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            if (i == 9) { // Bin [3.7-4.0) GeV/c
                gausPionBefore->SetParameters(betaPionsBefore->GetMaximum(), 0.996, 0.004);
                betaPionsBefore->Fit(gausPionBefore, "R", "", 0.995, 1.003);
            } 
            if (i == 5) { // Bin [2.5-2.8) GeV/c
                gausPionBefore->SetParameters(betaPionsBefore->GetMaximum(), 0.9996, 0.0047);
                betaPionsBefore->Fit(gausPionBefore, "R", "", 0.994, 1.006);
            }else {
                gausPionBefore->SetParameters(betaPionsBefore->GetMaximum(), 1.0, 0.005);
                betaPionsBefore->Fit(gausPionBefore, "R", "", betaFitRanges[i].first, betaFitRanges[i].second);
            }
            pMeanBefore = gausPionBefore->GetParameter(1);
            pSigmaBefore = gausPionBefore->GetParameter(2);
            pConstBefore = gausPionBefore->GetParameter(0);
            delete gausPionBefore;
        }
        if (betaProtonsBefore->GetEntries() >= 10) {
            TF1* gausProtonBefore = new TF1(TString::Format("gaus_%s", betaProtonsBefore->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            if (i == 5) { // Lower momentum bins (1.0-2.5 GeV/c)
                gausProtonBefore->SetParameters(betaProtonsBefore->GetMaximum(), 0.85, 0.02);
                betaProtonsBefore->Fit(gausProtonBefore, "R", "", 0.95, 1.0);
            } else if (i == 6) { // Lower momentum bins (1.0-2.5 GeV/c)
                gausProtonBefore->SetParameters(betaProtonsBefore->GetMaximum(), 0.952, 0.02);
                betaProtonsBefore->Fit(gausProtonBefore, "R", "", 0.94, 0.93);
            } else if (i == 7) { // Lower momentum bins (1.0-2.5 GeV/c)
                gausProtonBefore->SetParameters(betaProtonsBefore->GetMaximum(), 0.961, 0.012);
                betaProtonsBefore->Fit(gausProtonBefore, "R", "", 0.95, 0.964);
            } else if (i == 8) { // Lower momentum bins (1.0-2.5 GeV/c)
                gausProtonBefore->SetParameters(betaProtonsBefore->GetMaximum(), 0.966, 0.02);
                betaProtonsBefore->Fit(gausProtonBefore, "R", "", 0.96, 0.967);
            }else if (i == 9) { // Lower momentum bins (1.0-2.5 GeV/c)
                gausProtonBefore->SetParameters(betaProtonsBefore->GetMaximum(), 0.972, 0.02);
                betaProtonsBefore->Fit(gausProtonBefore, "R", "", 0.96, 0.978);
            } else if (i == 10) { // Lower momentum bins (1.0-2.5 GeV/c)
                gausProtonBefore->SetParameters(betaProtonsBefore->GetMaximum(), 0.975, 0.02);
                betaProtonsBefore->Fit(gausProtonBefore, "R", "", 0.968, 0.978);
            }else if (i == 11) { // Lower momentum bins (1.0-2.5 GeV/c)
                gausProtonBefore->SetParameters(betaProtonsBefore->GetMaximum(), 0.978, 0.02);
                betaProtonsBefore->Fit(gausProtonBefore, "R", "", 0.97, 0.984);
            }else if (i == 12) { // Lower momentum bins (1.0-2.5 GeV/c)
                gausProtonBefore->SetParameters(betaProtonsBefore->GetMaximum(), 0.981, 0.02);
                betaProtonsBefore->Fit(gausProtonBefore, "R", "", 0.974, 0.984);
            }else if (i == 13) { // Lower momentum bins (1.0-2.5 GeV/c)
                gausProtonBefore->SetParameters(betaProtonsBefore->GetMaximum(), 0.982, 0.01);
                betaProtonsBefore->Fit(gausProtonBefore, "R", "", 0.974, 0.988);
            }else if (i == 14) { // Lower momentum bins (1.0-2.5 GeV/c)
                gausProtonBefore->SetParameters(betaProtonsBefore->GetMaximum(), 0.985, 0.01);
                betaProtonsBefore->Fit(gausProtonBefore, "R", "", 0.978, 0.988);
            }else if (i == 15) { // Lower momentum bins (1.0-2.5 GeV/c)
                gausProtonBefore->SetParameters(betaProtonsBefore->GetMaximum(), 0.986, 0.01);
                betaProtonsBefore->Fit(gausProtonBefore, "R", "", 0.98, 0.991);
            }else if (i == 16) { // Lower momentum bins (1.0-2.5 GeV/c)
                gausProtonBefore->SetParameters(betaProtonsBefore->GetMaximum(), 0.99, 0.02);
                betaProtonsBefore->Fit(gausProtonBefore, "R", "", 0.978, 0.992);
            }else if (i == 17) { // Lower momentum bins (1.0-2.5 GeV/c)
                gausProtonBefore->SetParameters(betaProtonsBefore->GetMaximum(), 0.99, 0.02);
                betaProtonsBefore->Fit(gausProtonBefore, "R", "", 0.978, 0.993);
            }else if (i == 18) { // Lower momentum bins (1.0-2.5 GeV/c)
                gausProtonBefore->SetParameters(betaProtonsBefore->GetMaximum(), 0.991, 0.02);
                betaProtonsBefore->Fit(gausProtonBefore, "R", "", 0.982, 0.993);
            }else if (i == 19) { // Lower momentum bins (1.0-2.5 GeV/c)
                gausProtonBefore->SetParameters(betaProtonsBefore->GetMaximum(), 0.994, 0.02);
                betaProtonsBefore->Fit(gausProtonBefore, "R", "", 0.98, 0.994);
            }
            
            protonMeanBefore = gausProtonBefore->GetParameter(1);
            protonSigmaBefore = gausProtonBefore->GetParameter(2);
            protonConstBefore = gausProtonBefore->GetParameter(0);
            delete gausProtonBefore;
        }

        // Plot beta histograms (before cut)
        double yMin = 0.0;
        double yMax = -1.0;
        canvasBetaBefore->Clear();
        if (betaPionsBefore->GetEntries() >= 10) {
            betaPionsBefore->SetMarkerStyle(20);
            betaPionsBefore->SetMarkerSize(1);
            betaPionsBefore->SetMarkerColor(kBlue);
            betaPionsBefore->Draw("P");
            yMax = betaPionsBefore->GetMaximum();
            if (auto* fit = betaPionsBefore->GetFunction(TString::Format("gaus_%s", betaPionsBefore->GetName()))) {
                fit->SetLineColor(kBlue);
                fit->SetLineStyle(kDashed);
                fit->Draw("SAME");
                fit->SetRange(0.95, 1.03);
            }
        }
        if (betaProtonsBefore->GetEntries() >= 10) {
            betaProtonsBefore->SetMarkerStyle(20);
            betaProtonsBefore->SetMarkerSize(1);
            betaProtonsBefore->SetMarkerColor(kRed);
            betaProtonsBefore->Draw(betaPionsBefore->GetEntries() >= 10 ? "P SAME" : "P");
            yMax = max(yMax, betaProtonsBefore->GetMaximum());
            if (auto* fit = betaProtonsBefore->GetFunction(TString::Format("gaus_%s", betaProtonsBefore->GetName()))) {
                fit->SetLineColor(kRed);
                fit->SetLineStyle(kDashed);
                fit->Draw("SAME");
                fit->SetRange(0.95, 1.02);
            }
        }
        yMax *= 1.2;
        if (betaPionsBefore->GetEntries() >= 10 || betaProtonsBefore->GetEntries() >= 10) {
            (betaPionsBefore->GetEntries() >= 10 ? betaPionsBefore : betaProtonsBefore)->GetYaxis()->SetRangeUser(yMin, yMax);
        }
        TLegend* leg = new TLegend(0.15, 0.75, 0.35, 0.88);
        if (betaPionsBefore->GetEntries() >= 10) leg->AddEntry(betaPionsBefore, "Pions", "p");
        if (betaProtonsBefore->GetEntries() >= 10) leg->AddEntry(betaProtonsBefore, "Protons", "p");
        leg->SetTextSize(0.03);
        leg->Draw();
        canvasBetaBefore->SetGrid();
        canvasBetaBefore->Update();
        canvasBetaBefore->Print("output/pdf/proton/beta_before_cut_linear.pdf");
        delete leg;

        // Calculate the proton-to-pion ratio before the chi2pid cut
        // - This represents the initial contamination of protons in the pion sample
        // - Uses the beta distributions to compute the overlap between pions and protons
        // - Saved to proton_to_pion_ratio_before_chi2pid.csv with columns: Momentum Bin, Proton Left, Contamination (%)
        double c1Before, c2Before;
        map<string, double> pFitBefore = {{"gaus_mean", pMeanBefore}, {"gaus_sigma", pSigmaBefore}, {"gaus_constant", pConstBefore}};
        map<string, double> protonFitBefore = {{"gaus_mean", protonMeanBefore}, {"gaus_sigma", protonSigmaBefore}, {"gaus_constant", protonConstBefore}};
        double contaminationBefore = calculateOverlapContamination(pFitBefore, protonFitBefore, pionLeftBefore[i], pionRightBefore[i], protonRightBefore[i], c1Before, c2Before);
        double protonLeftBefore = -1.0;
        if (c1Before >= protonMeanBefore && c1Before <= pMeanBefore) {
            protonLeftBefore = c1Before;
        } else if (c2Before >= protonMeanBefore && c2Before <= pMeanBefore) {
            protonLeftBefore = c2Before;
        } else if (contaminationBefore >= 0) {
            protonLeftBefore = (abs(c1Before - protonMeanBefore) < abs(c2Before - protonMeanBefore)) ? c1Before : c2Before;
            cout << "No intersection between protonMeanBefore = " << protonMeanBefore << " and pMeanBefore = " << pMeanBefore 
                 << ", selecting point closest to proton peak: " << protonLeftBefore << endl;
        }

        // Save the proton-to-pion ratio (before cut) to CSV
        if (contaminationBefore < 0 || protonLeftBefore == -1.0) {
            cout << "Bin [" << pLow << "-" << pHigh << "): Proton-to-pion ratio (before cut) calculation failed\n";
            csvFileBefore << pLow << "-" << pHigh << ",N/A,N/A\n";
        } else {
            cout << "Bin [" << pLow << "-" << pHigh << "): Before cut - Pion range: [" << pionLeftBefore[i] << ", " << pionRightBefore[i] 
                 << "], Proton range: [" << protonLeftBefore << ", " << protonRightBefore[i] << "], Proton-to-Pion Ratio: " 
                 << contaminationBefore << "%\n";
            csvFileBefore << pLow << "-" << pHigh << "," << protonLeftBefore << "," << contaminationBefore << "\n";
        }





        // --- After Chi2pid Cut (Mean ± 3 Sigma) ---
        // Create chi2pid histograms (after cut)
        TH1F* chi2PionsAfter = new TH1F(TString::Format("chi2_pions_after_%d", i),
                                        TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                        100, chi2HistRanges[i].first, chi2HistRanges[i].second);
        TH1F* chi2ProtonsAfter = new TH1F(TString::Format("chi2_protons_after_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                         100, chi2HistRanges[i].first, chi2HistRanges[i].second);

        // Fill chi2pid histograms (after cut) with momentum and ±3σ cut based on pion fit
        
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);
        /* for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) {
                chi2PionsAfter->Fill(recomputed_chi2pid);
            }
        } */
       //for test
       for (Long64_t j = 0; j < nEntriesPions; j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) {
                chi2PionsAfter->Fill(recomputed_chi2pid);
            }
        }
        treeProtons->SetBranchAddress("p", &p);
        treeProtons->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);
        /* for (Long64_t j = 0; j < treeProtons->GetEntries(); j++) {
            treeProtons->GetEntry(j);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) {
                chi2ProtonsAfter->Fill(recomputed_chi2pid);
            }
        } */

       //for test
       for (Long64_t j = 0; j < nEntriesProtons; j++) {
            treeProtons->GetEntry(j);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) {
                chi2ProtonsAfter->Fill(recomputed_chi2pid);
            }
        }

        // Set colors for chi2pid histograms (after cut)
        chi2PionsAfter->SetLineColor(kBlue);
        chi2ProtonsAfter->SetLineColor(kRed);

        // Fit chi2pid histograms (after cut) with a Gaussian
        TF1* fitPionAfter = new TF1("gausFit", "gaus", chi2FitRanges[i].first, chi2FitRanges[i].second);
        fitPionAfter->SetParameters(pionAmp, pionMeanInit, pionSigmaInit);
        chi2PionsAfter->Fit(fitPionAfter, "R", "", chi2FitRanges[i].first, chi2FitRanges[i].second);
        TF1* fitProtonAfter = new TF1("gausFit", "gaus", chi2FitRanges[i].first, chi2FitRanges[i].second);
        fitProtonAfter->SetParameters(protonAmp, protonMean, protonSigma);
        chi2ProtonsAfter->Fit(fitProtonAfter, "R", "", chi2FitRanges[i].first, chi2FitRanges[i].second);

        // Plot chi2pid histograms (after cut)
        canvasChi2After->Clear();
        canvasChi2After->Divide(2, 1);
        canvasChi2After->cd(1);
        chi2PionsAfter->Draw("HIST");
        if (chi2PionsAfter->GetFunction("gausFit")) {
            chi2PionsAfter->GetFunction("gausFit")->SetLineColor(kBlue);
            chi2PionsAfter->GetFunction("gausFit")->Draw("SAME");
        }
        TLegend* legPionsAfter = new TLegend(0.1, 0.75, 0.3, 0.9);
        legPionsAfter->AddEntry(chi2PionsAfter, "Pions", "l");
        legPionsAfter->SetBorderSize(0);
        legPionsAfter->SetFillStyle(0);
        legPionsAfter->SetTextSize(0.025);
        legPionsAfter->Draw();
        canvasChi2After->cd(2);
        chi2ProtonsAfter->Draw("HIST");
        if (chi2ProtonsAfter->GetFunction("gausFit")) {
            chi2ProtonsAfter->GetFunction("gausFit")->SetLineColor(kRed);
            chi2ProtonsAfter->GetFunction("gausFit")->Draw("SAME");
        }
        TLegend* legProtonsAfter = new TLegend(0.1, 0.75, 0.3, 0.9);
        legProtonsAfter->AddEntry(chi2ProtonsAfter, "Protons (Pion Hypothesis)", "l");
        legProtonsAfter->SetBorderSize(0);
        legProtonsAfter->SetFillStyle(0);
        legProtonsAfter->SetTextSize(0.025);
        legProtonsAfter->Draw();
        canvasChi2After->Update();
        canvasChi2After->Print("output/pdf/proton/chi2pid_fits_linear_after.pdf");
        delete legPionsAfter;
        delete legProtonsAfter;

        // Create beta histograms (after cut)
        TH1F* betaPionsAfter = new TH1F(TString::Format("beta_pions_after_%d", i),
                                        TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                        70, betaHistRanges[i].first, betaHistRanges[i].second);
        TH1F* betaProtonsAfter = new TH1F(TString::Format("beta_protons_after_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                         70, betaHistRanges[i].first, betaHistRanges[i].second);
        betaPionsAfter->Sumw2();
        betaProtonsAfter->Sumw2();

        // Fill beta histograms (after cut) with momentum and ±3σ cut
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("beta", &beta);
        treePions->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);
        //for test
        /* for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) {
                betaPionsAfter->Fill(beta);
            }
        } */
       //for test
       for (Long64_t j = 0; j < nEntriesPions; j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) {
                betaPionsAfter->Fill(beta);
            }
        }

        treeProtons->SetBranchAddress("p", &p);
        treeProtons->SetBranchAddress("beta", &beta);
        treeProtons->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);
        /* for (Long64_t j = 0; j < treeProtons->GetEntries(); j++) {
            treeProtons->GetEntry(j);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) {
                betaProtonsAfter->Fill(beta);
            }
        } */

       //for test
        for (Long64_t j = 0; j < nEntriesProtons; j++) {
            treeProtons->GetEntry(j);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) {
                betaProtonsAfter->Fill(beta);
            }
        }
        // Fit beta histograms (after cut) for contamination calculation
        double pMeanAfter = 0, pSigmaAfter = 0, pConstAfter = 0, protonMeanAfter = 0, protonSigmaAfter = 0, protonConstAfter = 0;
        if (betaPionsAfter->GetEntries() >= 10) {
            TF1* gausPionAfter = new TF1(TString::Format("gaus_%s", betaPionsAfter->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            if (i <= 3) { 
                gausPionAfter->SetParameters(betaPionsAfter->GetMaximum(), 0.992, 0.02);
                betaPionsAfter->Fit(gausPionAfter, "R", "", 0.984, 1.03);
            } else {
                gausPionAfter->SetParameters(betaPionsAfter->GetMaximum(), 1.0, 0.005);
                betaPionsAfter->Fit(gausPionAfter, "R", "", 0.984, 1.03);
            }
            pMeanAfter = gausPionAfter->GetParameter(1);
            pSigmaAfter = gausPionAfter->GetParameter(2);
            pConstAfter = gausPionAfter->GetParameter(0);
            delete gausPionAfter;
        }
        if (betaProtonsAfter->GetEntries() >= 10) {
            TF1* gausProtonAfter = new TF1(TString::Format("gaus_%s", betaProtonsAfter->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            if (i <= 4) { // Lower momentum bins (1.0-2.5 GeV/c)
                gausProtonAfter->SetParameters(betaProtonsAfter->GetMaximum(), 0.85, 0.002);
                betaProtonsAfter->Fit(gausProtonAfter, "R", "", 0.95, 1.02);
            } else if (i == 10) { // Mid momentum bins (2.5-4.0 GeV/c)
                gausProtonAfter->SetParameters(betaProtonsAfter->GetMaximum(), 0.90, 0.002);
                betaProtonsAfter->Fit(gausProtonAfter, "R", "", 0.95, 1.02);
            }  else if (i == 11) { // Mid momentum bins (2.5-4.0 GeV/c)
                gausProtonAfter->SetParameters(betaProtonsAfter->GetMaximum(), 0.90, 0.002);
                betaProtonsAfter->Fit(gausProtonAfter, "R", "", 0.978, 0.986);
            }else if (i == 15) { // Mid momentum bins (2.5-4.0 GeV/c)
                gausProtonAfter->SetParameters(betaProtonsAfter->GetMaximum(), 0.991, 0.002);
                betaProtonsAfter->Fit(gausProtonAfter, "R", "", 0.988, 0.991);
            }
            else if (i == 16) { // Mid momentum bins (2.5-4.0 GeV/c)
                gausProtonAfter->SetParameters(betaProtonsAfter->GetMaximum(), 0.90, 0.002);
                betaProtonsAfter->Fit(gausProtonAfter, "R", "", 0.988, 0.993);
            }else if (i == 17) { // Mid momentum bins (2.5-4.0 GeV/c)
                gausProtonAfter->SetParameters(betaProtonsAfter->GetMaximum(), 0.992, 0.002);
                betaProtonsAfter->Fit(gausProtonAfter, "R", "", 0.988, 0.993);
            }else if (i == 18) { // Mid momentum bins (2.5-4.0 GeV/c)
                gausProtonAfter->SetParameters(betaProtonsAfter->GetMaximum(), 0.993, 0.002);
                betaProtonsAfter->Fit(gausProtonAfter, "R", "", 0.99, 0.994);
            }else if (i == 19) { // Mid momentum bins (2.5-4.0 GeV/c)
                gausProtonAfter->SetParameters(betaProtonsAfter->GetMaximum(), 0.993, 0.002);
                betaProtonsAfter->Fit(gausProtonAfter, "R", "", 0.99, 0.994);
            }else { // Higher momentum bins (4.0-7.0 GeV/c)
                gausProtonAfter->SetParameters(betaProtonsAfter->GetMaximum(), 0.95, 0.002);
                betaProtonsAfter->Fit(gausProtonAfter, "R", "", 0.95, 1.02);
            }
            protonMeanAfter = gausProtonAfter->GetParameter(1);
            protonSigmaAfter = gausProtonAfter->GetParameter(2);
            protonConstAfter = gausProtonAfter->GetParameter(0);
            delete gausProtonAfter;
        }

        // Plot beta histograms (after cut)
        yMin = 0.0;
        yMax = -1.0;
        canvasBetaAfter->Clear();
        if (betaPionsAfter->GetEntries() >= 10) {
            betaPionsAfter->SetMarkerStyle(20);
            betaPionsAfter->SetMarkerSize(1);
            betaPionsAfter->SetMarkerColor(kBlue);
            betaPionsAfter->Draw("P");
            yMax = betaPionsAfter->GetMaximum();
            if (auto* fit = betaPionsAfter->GetFunction(TString::Format("gaus_%s", betaPionsAfter->GetName()))) {
                fit->SetLineColor(kBlue);
                fit->SetLineStyle(kDashed);
                fit->Draw("SAME");
                fit->SetRange(0.95, 1.03);
            }
        }
        if (betaProtonsAfter->GetEntries() >= 10) {
            betaProtonsAfter->SetMarkerStyle(20);
            betaProtonsAfter->SetMarkerSize(1);
            betaProtonsAfter->SetMarkerColor(kRed);
            betaProtonsAfter->Draw(betaPionsAfter->GetEntries() >= 10 ? "P SAME" : "P");
            yMax = max(yMax, betaProtonsAfter->GetMaximum());
            if (auto* fit = betaProtonsAfter->GetFunction(TString::Format("gaus_%s", betaProtonsAfter->GetName()))) {
                fit->SetLineColor(kRed);
                fit->SetLineStyle(kDashed);
                fit->Draw("SAME");
                fit->SetRange(0.75, 1.03);
            }
        }
        yMax *= 1.2;
        if (betaPionsAfter->GetEntries() >= 10 || betaProtonsAfter->GetEntries() >= 10) {
            (betaPionsAfter->GetEntries() >= 10 ? betaPionsAfter : betaProtonsAfter)->GetYaxis()->SetRangeUser(yMin, yMax);
        }
        TLegend* legAfter = new TLegend(0.15, 0.75, 0.35, 0.88);
        if (betaPionsAfter->GetEntries() >= 10) legAfter->AddEntry(betaPionsAfter, "Pions", "p");
        if (betaProtonsAfter->GetEntries() >= 10) legAfter->AddEntry(betaProtonsAfter, "Protons (Misidentified)", "p");
        legAfter->SetTextSize(0.03);
        legAfter->Draw();
        canvasBetaAfter->SetGrid();
        canvasBetaAfter->Update();
        canvasBetaAfter->Print("output/pdf/proton/beta_after_cut_linear.pdf");
        delete legAfter;

        // Calculate the contamination after the chi2pid cut
        // - This represents the remaining proton-to-pion ratio in the pion sample after applying the ±3σ chi2pid cut
        // - Uses the beta distributions to compute the overlap between pions and protons after the cut
        // - Saved to contamination_after_chi2pid.csv with columns: Momentum Bin, Proton Left, Proton-to-Pion Ratio (%)
        double c1After, c2After;
        map<string, double> pFitAfter = {{"gaus_mean", pMeanAfter}, {"gaus_sigma", pSigmaAfter}, {"gaus_constant", pConstAfter}};
        map<string, double> protonFitAfter = {{"gaus_mean", protonMeanAfter}, {"gaus_sigma", protonSigmaAfter}, {"gaus_constant", protonConstAfter}};
        double protonToPionRatioAfter = calculateOverlapContamination(pFitAfter, protonFitAfter, pionLeftAfter[i], pionRightAfter[i], protonRightAfter[i], c1After, c2After);
        double protonLeftAfter = -1.0;
        if (c1After >= protonMeanAfter && c1After <= pMeanAfter) {
            protonLeftAfter = c1After;
        } else if (c2After >= protonMeanAfter && c2After <= pMeanAfter) {
            protonLeftAfter = c2After;
        } else if (protonToPionRatioAfter >= 0) {
            protonLeftAfter = (abs(c1After - protonMeanAfter) < abs(c2After - protonMeanAfter)) ? c1After : c2After;
            cout << "No intersection between protonMeanAfter = " << protonMeanAfter << " and pMeanAfter = " << pMeanAfter 
                 << ", selecting point closest to proton peak: " << protonLeftAfter << endl;
        }

        // Save the contamination (after cut) to CSV
        if (protonToPionRatioAfter < 0 || protonLeftAfter == -1.0) {
            cout << "Bin [" << pLow << "-" << pHigh << "): Contamination (after cut) calculation failed\n";
            csvFileAfter << pLow << "-" << pHigh << ",N/A,N/A\n";
        } else {
            cout << "Bin [" << pLow << "-" << pHigh << "): After cut - Pion range: [" << pionLeftAfter[i] << ", " << pionRightAfter[i] 
                 << "], Proton range: [" << protonLeftAfter << ", " << protonRightAfter[i] << "], Contamination: " 
                 << protonToPionRatioAfter << "%\n";
            csvFileAfter << pLow << "-" << pHigh << "," << protonLeftAfter << "," << protonToPionRatioAfter << "\n";
        }

        // Clean up histograms and fits for this bin
        delete chi2PionsBefore;
        delete chi2ProtonsBefore;
        delete betaPionsBefore;
        delete betaProtonsBefore;
        delete chi2PionsAfter;
        delete chi2ProtonsAfter;
        delete betaPionsAfter;
        delete betaProtonsAfter;
        delete fitPion;
        delete fitProton;
        delete fitPionAfter;
        delete fitProtonAfter;
    }

    // Step 10: Close PDF and CSV files
    canvasChi2Before->Print("output/pdf/proton/chi2pid_fits_linear_before.pdf]");
    canvasChi2After->Print("output/pdf/proton/chi2pid_fits_linear_after.pdf]");
    canvasBetaBefore->Print("output/pdf/proton/beta_before_cut_linear.pdf]");
    canvasBetaAfter->Print("output/pdf/proton/beta_after_cut_linear.pdf]");
    csvFileBefore.close();
    csvFileAfter.close();
    csvFileCuts.close();

    // Step 11: Clean up
    delete canvasChi2Before;
    delete canvasChi2After;
    delete canvasBetaBefore;
    delete canvasBetaAfter;
    file->Close();
    delete file;

    return 0;
}  