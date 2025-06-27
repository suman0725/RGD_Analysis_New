/* #include "include/ContaminationUtils.h"
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
    // Step 1: Open the ROOT file and load the trees for pions and kaons
    TFile* file = new TFile("../Skim/pkptreeCxC_9_test_modified.root");
    //TFile* file = new TFile("../Skim/pkptreeCxC_9_test.root");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file pkptreeCxC_9_test_modified.root\n";
        return 1;
    }
    TTree* treePions = (TTree*)file->Get("EB_pid_pions");
    TTree* treeKaons = (TTree*)file->Get("EB_pid_kaons");
    TTree* treePionsall = (TTree*)file->Get("EB_all_pion_assumed"); 
    TTree* treeEBPosPionAssumed = (TTree*)file->Get("EB_pos_pion_assumed"); 
    if (!treePions || !treeKaons || !treePionsall || !treeEBPosPionAssumed) {
        cerr << "Error: Cannot load trees\n";
        file->Close();
        delete file;
        return 1;
    }

    // Step 2: Create output directories for saving PDFs and CSVs
    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/pdf/kaon", kTRUE);
    gSystem->mkdir("output/csv/kaon", kTRUE);

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
            chi2HistRanges[i] = {-5.0, 5.0};
            chi2FitRanges[i] = {-1.3, 1.3};
            chi2FitParams[i] = {1000.0, 0.0, 1.0};
        } else if (pLow < 3.0) {
            chi2HistRanges[i] = {-4.0, 4.0};
            chi2FitRanges[i] = {-1.5, 1.5};
            chi2FitParams[i] = {800.0, 0.1, 0.8};
        } else {
            chi2HistRanges[i] = {-3.0, 3.0};
            chi2FitRanges[i] = {-1, 1};
            chi2FitParams[i] = {500.0, -0.1, 0.5};
        }
    }

    // Step 5: Set beta histogram ranges and fit ranges for each bin
    vector<pair<double, double>> betaHistRanges(nBins);
    vector<pair<double, double>> betaFitRanges(nBins);
    for (int i = 0; i < nBins; i++) {
        betaHistRanges[i] = {0.95, 1.03};
        betaFitRanges[i] = {0.95, 1.02};
    }

    // Step 6: Set beta cuts for pions and kaons (same for before and after chi2pid cuts)
    vector<double> pionLeftBefore(nBins), pionRightBefore(nBins), kaonRightBefore(nBins);
    vector<double> pionLeftAfter(nBins), pionRightAfter(nBins), kaonRightAfter(nBins);
    for (int i = 0; i < nBins; i++) {


        if (i == 0) { // Bin [1.0-1.3)
            pionLeftBefore[i] = 0.9880; pionRightBefore[i] = 1.0120; kaonRightBefore[i] = 0.0;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  kaonRightAfter[i] = 0.0;
        } else if (i == 1) { // Bin [1.3-1.6)
            pionLeftBefore[i] = 0.9755; pionRightBefore[i] = 1.0135; kaonRightBefore[i] = 0.9721;
            pionLeftAfter[i] = 0.9755;  pionRightAfter[i] = 1.0135;  kaonRightAfter[i] = 0.9721;
        } else if (i == 2) { // Bin [1.6-1.9)
            pionLeftBefore[i] = 0.98;   pionRightBefore[i] = 1.0140; kaonRightBefore[i] = 0.986;
            pionLeftAfter[i] = 0.98;    pionRightAfter[i] = 1.0140;  kaonRightAfter[i] = 0.986;
        } else if (i == 3) { // Bin [1.9-2.2)
            pionLeftBefore[i] = 0.9815; pionRightBefore[i] = 1.0145; kaonRightBefore[i] = 0.992;
            pionLeftAfter[i] = 0.9815;  pionRightAfter[i] = 1.0145;  kaonRightAfter[i] = 0.992;
        } else if (i == 4) { // Bin [2.2-2.5)
            pionLeftBefore[i] = 0.9880; pionRightBefore[i] = 1.0141; kaonRightBefore[i] = 0.9971;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0141;  kaonRightAfter[i] = 0.9971;
        } else if (i == 5) { // Bin [2.5-2.8)
            pionLeftBefore[i] = 0.9848; pionRightBefore[i] = 1.013;  kaonRightBefore[i] = 1.0;
            pionLeftAfter[i] = 0.9848;  pionRightAfter[i] = 1.013;   kaonRightAfter[i] = 1.0;
        } else if (i == 6) { // Bin [2.8-3.1)
            pionLeftBefore[i] = 0.986;  pionRightBefore[i] = 1.0128; kaonRightBefore[i] = 1.003;
            pionLeftAfter[i] = 0.986;   pionRightAfter[i] = 1.0128;  kaonRightAfter[i] = 1.003;
        } else if (i == 7) { // Bin [3.1-3.4)
            pionLeftBefore[i] = 0.9868; pionRightBefore[i] = 1.0121; kaonRightBefore[i] = 1.0048;
            pionLeftAfter[i] = 0.9868;  pionRightAfter[i] = 1.0121;  kaonRightAfter[i] = 1.0048;
        } else if (i == 8) { // Bin [3.4-3.7)
            pionLeftBefore[i] = 0.9879; pionRightBefore[i] = 1.0120; kaonRightBefore[i] = 1.0049;
            pionLeftAfter[i] = 0.9879;  pionRightAfter[i] = 1.0120;  kaonRightAfter[i] = 1.0049;
        } else if (i == 9) { // Bin [3.7-4.0)
            pionLeftBefore[i] = 0.9880; pionRightBefore[i] = 1.0120; kaonRightBefore[i] = 1.012;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  kaonRightAfter[i] = 1.012;
        } else {
            cerr << "Warning: Using default beta cuts for bin [" << pBins[i] << "-" << pBins[i+1] << ")\n";
            pionLeftBefore[i] = 0.9880; pionRightBefore[i] = 1.0120; kaonRightBefore[i] = 1.0;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  kaonRightAfter[i] = 1.0;
        }
    }

    // Step 7: Create canvases and open PDF files for plotting
    TCanvas* canvasChi2Before = new TCanvas("canvasChi2Before", "Chi2pid Fits (Before)", 1900, 900);
    canvasChi2Before->SetMargin(0.1, 0.05, 0.15, 0.1);
    TCanvas* canvasChi2After = new TCanvas("canvasChi2After", "Chi2pid Fits (After)", 1900, 900);
     canvasChi2After->SetMargin(0.1, 0.05, 0.15, 0.1);
    TCanvas* canvasBetaBefore = new TCanvas("canvasBetaBefore", "Beta Plots (Before)", 1200, 800);
    canvasBetaBefore->SetMargin(0.1, 0.05, 0.15, 0.1);
    TCanvas* canvasBetaAfter = new TCanvas("canvasBetaAfter", "Beta Plots (After)", 800, 600);
    canvasBetaAfter->SetMargin(0.1, 0.05, 0.15, 0.1);
    TCanvas* canvaschi2pidall = new TCanvas("canvaschi2pidall", "chi2pid Plot", 1200, 800); 
    TCanvas* canvaschi2pidallafter = new TCanvas("canvaschi2pidallafter", "chi2pid Plot (After)", 1200, 800); 
    TCanvas* canvaschi2pidposbefore = new TCanvas("canvaschi2pidpos", "chi2pid Plot +ve Charge Particles (Before)", 1200, 800); 
    TCanvas* canvaschi2pidposafter = new TCanvas("canvaschi2pidposafter", "chi2pid Plot +ve Charge Particles (After)", 1200, 800);

    canvaschi2pidall->Print("output/pdf/kaon/chi2pidforall.pdf[");
    canvaschi2pidallafter->Print("output/pdf/kaon/chi2pidforallafter.pdf[");
    canvasChi2Before->Print("output/pdf/kaon/chi2pid_fits_linear_before.pdf[");
    canvasChi2After->Print("output/pdf/kaon/chi2pid_fits_linear_after.pdf[");
    canvasBetaBefore->Print("output/pdf/kaon/beta_before_cut_linear.pdf[");
    canvasBetaAfter->Print("output/pdf/kaon/beta_after_cut_linear.pdf[");
    canvaschi2pidposbefore->Print("output/pdf/kaon/chi2pidpos_before.pdf["); 
    canvaschi2pidposafter->Print("output/pdf/kaon/chi2pidpos_after.pdf["); 

    // Step 8: Open CSV files to save results
    ofstream csvFileBefore("output/csv/kaon/kaon_to_pion_ratio_before_chi2pid.csv");
    csvFileBefore << "Momentum Bin (GeV/c),Kaon Left,Contamination (%)\n";
    ofstream csvFileAfter("output/csv/kaon/contamination_after_chi2pid.csv");
    csvFileAfter << "Momentum Bin (GeV/c),Kaon Left,Kaon-to-Pion Ratio (%)\n";
    ofstream csvFileCuts("output/csv/kaon/chi2pid_cuts.csv");
    csvFileCuts << "Momentum Bin (GeV/c),Pion Mean,Pion Sigma,Chi2 Min (Mean-3σ),Chi2 Max (Mean+3σ)\n";

    // Step 9: Loop over each momentum bin to process data
    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];




        // --- Before Chi2pid Cut ---

        
        // Create chi2pid histograms for pions, kaons, and proton (EB PID: 211, 321, 2212) under pion hypothesis 

        TH1F* chi2pidall = new TH1F(TString::Format("chi2pidall_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                         100, -15, 15);

        // Fiil chi2pidall histograms 
        float p, chi2pid, beta, orig_chi2pid, recomputed_chi2pid;

        treePionsall->SetBranchAddress("p", &p); 
        treePionsall->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid); 
        for ( Long64_t j = 0; j < treePionsall->GetEntries(); j++){
            treePionsall->GetEntry(j); 
            if (p >= pLow && p < pHigh) chi2pidall->Fill(recomputed_chi2pid); 
        }

        // Plot chi2pidall in canvas 
        canvaschi2pidall->Clear(); 
        chi2pidall-> Draw("HIST"); 
        canvaschi2pidall->Update(); 
        canvaschi2pidall->Print("output/pdf/kaon/chi2pidforall.pdf");

        //--------------------------------//--------------------------------------------------------//-----------------------------//

        // Create chi2pid histograms for all +ve charge particles under pion hypothesis 

        TH1F* chi2pidposbefore = new TH1F(TString::Format("chi2pidposbefore_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                         100, -15, 15);

        // Fiil chi2pidall histograms 
       
        treeEBPosPionAssumed->SetBranchAddress("p", &p); 
        treeEBPosPionAssumed->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid); 
        for ( Long64_t j = 0; j < treePionsall->GetEntries(); j++){
            treeEBPosPionAssumed->GetEntry(j); 
            if (p >= pLow && p < pHigh) chi2pidposbefore->Fill(recomputed_chi2pid); 
        }

        // Plot chi2pidall in canvas 
        canvaschi2pidposbefore->Clear(); 
        chi2pidposbefore-> Draw("HIST"); 
        canvaschi2pidposbefore->Update(); 
        canvaschi2pidposbefore->Print("output/pdf/kaon/chi2pidpos_before.pdf");

        //--------------------------------//--------------------------------------------------------//-----------------------------//

        // Create chi2pid histograms for pions and kaons (before cut)
        TH1F* chi2PionsBefore = new TH1F(TString::Format("chi2_pions_before_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                         100, -5, 5);
        TH1F* chi2KaonsBefore = new TH1F(TString::Format("chi2_kaons_before_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                         100, -5, 5);

        // Fill chi2pid histograms (before cut) with momentum cut
        
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("orig_chi2pid", &orig_chi2pid);
        for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh) chi2PionsBefore->Fill(orig_chi2pid);
        }
        treeKaons->SetBranchAddress("p", &p);
        treeKaons->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);
        for (Long64_t j = 0; j < treeKaons->GetEntries(); j++) {
            treeKaons->GetEntry(j);
            if (p >= pLow && p < pHigh) chi2KaonsBefore->Fill(recomputed_chi2pid);
        }

        // Set colors for chi2pid histograms
        chi2PionsBefore->SetLineColor(kBlue);
        chi2KaonsBefore->SetLineColor(kGreen);

        // Fit chi2pid histograms with a Gaussian to determine the ±3σ cut
        double pionMeanBefore, pionSigmaBefore, pionChi2Min, pionChi2Max;
        auto [pionAmp, pionMeanInit, pionSigmaInit] = chi2FitParams[i];
        auto [kaonAmp, kaonMean, kaonSigma] = chi2FitParams[i];
        TF1* fitPion = new TF1("gausFit", "gaus", chi2FitRanges[i].first, chi2FitRanges[i].second);
        fitPion->SetParameters(pionAmp, pionMeanInit, pionSigmaInit);
        chi2PionsBefore->Fit(fitPion, "R", "", chi2FitRanges[i].first, chi2FitRanges[i].second);
        pionMeanBefore = fitPion->GetParameter(1);
        pionSigmaBefore = fitPion->GetParameter(2);
        pionChi2Min = pionMeanBefore - 5 * pionSigmaBefore;
        pionChi2Max = pionMeanBefore + 5 * pionSigmaBefore;
        //pionChi2Min = pionMeanBefore - 3 * pionSigmaBefore;
        //pionChi2Max = pionMeanBefore + 3 * pionSigmaBefore;
        
        cout << "Bin [" << pLow << "-" << pHigh << "): Pions chi2 range (before): [" 
             << pionChi2Min << ", " << pionChi2Max << "]\n";

        // Save the pion chi2pid fit parameters to CSV
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
        chi2KaonsBefore->Draw("HIST");
        if (chi2KaonsBefore->GetFunction("gausFit")) {
            chi2KaonsBefore->GetFunction("gausFit")->SetLineColor(kGreen);
            chi2KaonsBefore->GetFunction("gausFit")->Draw("SAME");
        }
        TLegend* legKaons = new TLegend(0.1, 0.75, 0.3, 0.9);
        legKaons->AddEntry(chi2KaonsBefore, "Kaons", "l");
        legKaons->SetBorderSize(0);
        legKaons->SetFillStyle(0);
        legKaons->SetTextSize(0.025);
        legKaons->Draw();
        canvasChi2Before->Update();
        canvasChi2Before->Print("output/pdf/kaon/chi2pid_fits_linear_before.pdf");
        delete legPions;
        delete legKaons;

        
        //--------------------------------//--------------------------------------------------------//-----------------------------//


        // Create beta histograms (before cut)
        TH1F* betaPionsBefore = new TH1F(TString::Format("beta_pions_before_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                         70, betaHistRanges[i].first, betaHistRanges[i].second);
        TH1F* betaKaonsBefore = new TH1F(TString::Format("beta_kaons_before_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                         70, betaHistRanges[i].first, betaHistRanges[i].second);
        betaPionsBefore->Sumw2();
        betaKaonsBefore->Sumw2();

        // Fill beta histograms (before cut) with momentum cut
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("beta", &beta);
        for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh) {
                betaPionsBefore->Fill(beta);
            }
        }
        treeKaons->SetBranchAddress("p", &p);
        treeKaons->SetBranchAddress("beta", &beta);
        for (Long64_t j = 0; j < treeKaons->GetEntries(); j++) {
            treeKaons->GetEntry(j);
            if (p >= pLow && p < pHigh) {
                betaKaonsBefore->Fill(beta);
            }
        }

        // Fit beta histograms (before cut) to calculate the kaon-to-pion ratio
        double pMeanBefore = 0, pSigmaBefore = 0, pConstBefore = 0, kMeanBefore = 0, kSigmaBefore = 0, kConstBefore = 0;
        if (betaPionsBefore->GetEntries() >= 10) {
            TF1* gausPionBefore = new TF1(TString::Format("gaus_%s", betaPionsBefore->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            if (i == 9) { // Bin [3.7-4.0) GeV/c
                gausPionBefore->SetParameters(betaPionsBefore->GetMaximum(), 0.996, 0.004);
                betaPionsBefore->Fit(gausPionBefore, "R", "", 0.994, 1.006);
            } else {
                gausPionBefore->SetParameters(betaPionsBefore->GetMaximum(), 1.0, 0.005);
                betaPionsBefore->Fit(gausPionBefore, "R", "", betaFitRanges[i].first, betaFitRanges[i].second);
            }
            pMeanBefore = gausPionBefore->GetParameter(1);
            pSigmaBefore = gausPionBefore->GetParameter(2);
            pConstBefore = gausPionBefore->GetParameter(0);
            delete gausPionBefore;
        }
        if (betaKaonsBefore->GetEntries() >= 10) {
            TF1* gausKaonBefore = new TF1(TString::Format("gaus_%s", betaKaonsBefore->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            if (i == 9) { // Bin [3.7-4.0) GeV/c
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.986, 0.9945);
            } else if (i == 8) { // Bin [3.4-3.7) GeV/c
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.988, 0.9945);
            } else if (i == 10) { // Bin [4.0-4.3) GeV/c
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.987, 0.995);
            } else if (i == 11) { // Bin [4.3-4.6) GeV/c
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.996, 0.012);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.992, 0.996);
            } else if (i == 12) { // Bin [4.6-4.9) GeV/c
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.012);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.991, 0.996);
            } else if (i == 13) { // Bin [4.9-5.2) GeV/c
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.012);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.991, 0.995);
            } else if (i == 14) { // Bin [5.2-5.5) GeV/c
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.991, 0.996);
            } else if (i == 15) { // Bin [5.5-5.8) GeV/c
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.991, 0.996);
            } else if (i == 16) { // Bin [5.8-6.1) GeV/c
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.991, 0.997);
            } else if (i == 17) { // Bin [6.1-6.4) GeV/c
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.991, 0.997);
            } else if (i == 18) { // Bin [6.4-6.7) GeV/c
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.982, 0.997);
            } else if (i == 5) { // Bin [2.5-2.8) GeV/c
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.974, 0.992);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.974, 0.992);
            } else {
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.99, 0.005);
                betaKaonsBefore->Fit(gausKaonBefore, "R");
            }
            kMeanBefore = gausKaonBefore->GetParameter(1);
            kSigmaBefore = gausKaonBefore->GetParameter(2);
            kConstBefore = gausKaonBefore->GetParameter(0);
            delete gausKaonBefore;
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
        if (betaKaonsBefore->GetEntries() >= 10) {
            betaKaonsBefore->SetMarkerStyle(20);
            betaKaonsBefore->SetMarkerSize(1);
            betaKaonsBefore->SetMarkerColor(kGreen);
            betaKaonsBefore->Draw(betaPionsBefore->GetEntries() >= 10 ? "P SAME" : "P");
            yMax = max(yMax, betaKaonsBefore->GetMaximum());
            if (auto* fit = betaKaonsBefore->GetFunction(TString::Format("gaus_%s", betaKaonsBefore->GetName()))) {
                fit->SetLineColor(kGreen);
                fit->SetLineStyle(kDashed);
                fit->Draw("SAME");
                fit->SetRange(0.95, 1.03);
            }
        }
        yMax *= 1.2;
        if (betaPionsBefore->GetEntries() >= 10 || betaKaonsBefore->GetEntries() >= 10) {
            (betaPionsBefore->GetEntries() >= 10 ? betaPionsBefore : betaKaonsBefore)->GetYaxis()->SetRangeUser(yMin, yMax);
        }
        TLegend* leg = new TLegend(0.15, 0.75, 0.35, 0.88);
        if (betaPionsBefore->GetEntries() >= 10) leg->AddEntry(betaPionsBefore, "Pions", "p");
        if (betaKaonsBefore->GetEntries() >= 10) leg->AddEntry(betaKaonsBefore, "Kaons", "p");
        leg->SetTextSize(0.03);
        leg->Draw();
        canvasBetaBefore->SetGrid();
        canvasBetaBefore->Update();
        canvasBetaBefore->Print("output/pdf/kaon/beta_before_cut_linear.pdf");
        delete leg;

        //--------------------------------//--------------------------------------------------------//-----------------------------//


        // Calculate the kaon-to-pion ratio before the chi2pid cut
        double c1Before, c2Before;
        map<string, double> pFitBefore = {{"gaus_mean", pMeanBefore}, {"gaus_sigma", pSigmaBefore}, {"gaus_constant", pConstBefore}};
        map<string, double> kFitBefore = {{"gaus_mean", kMeanBefore}, {"gaus_sigma", kSigmaBefore}, {"gaus_constant", kConstBefore}};
        double contaminationBefore = calculateOverlapContamination(pFitBefore, kFitBefore, pionLeftBefore[i], pionRightBefore[i], kaonRightBefore[i], c1Before, c2Before);
        double kaonLeftBefore = -1.0;
        if (c1Before >= kMeanBefore && c1Before <= pMeanBefore) {
            kaonLeftBefore = c1Before;
        } else if (c2Before >= kMeanBefore && c2Before <= pMeanBefore) {
            kaonLeftBefore = c2Before;
        } else if (contaminationBefore >= 0) {
            kaonLeftBefore = (abs(c1Before - kMeanBefore) < abs(c2Before - kMeanBefore)) ? c1Before : c2Before;
            cout << "No intersection between kMeanBefore = " << kMeanBefore << " and pMeanBefore = " << pMeanBefore 
                 << ", selecting point closest to kaon peak: " << kaonLeftBefore << endl;
        }

        // Save the kaon-to-pion ratio (before cut) to CSV
        if (contaminationBefore < 0 || kaonLeftBefore == -1.0) {
            cout << "Bin [" << pLow << "-" << pHigh << "): Kaon-to-pion ratio (before cut) calculation failed\n";
            csvFileBefore << pLow << "-" << pHigh << ",N/A,N/A\n";
        } else {
            cout << "Bin [" << pLow << "-" << pHigh << "): Before cut - Pion range: [" << pionLeftBefore[i] << ", " << pionRightBefore[i] 
                 << "], Kaon range: [" << kaonLeftBefore << ", " << kaonRightBefore[i] << "], Kaon-to-Pion Ratio: " 
                 << contaminationBefore << "%\n";
            csvFileBefore << pLow << "-" << pHigh << "," << kaonLeftBefore << "," << contaminationBefore << "\n";
        }



        //  ################-----------------------################--------------------################------------------------################
        //    ################-----------------------################--------------------################------------------------################
        //   ################-----------------------################--------------------################------------------------################
        //    ################-----------------------################--------------------################------------------------################
         //   ################-----------------------################--------------------################------------------------################  
        

        // --- After Chi2pid Cut (Mean ± 3 Sigma) ---
        // ---Fill Chipid for all after ---

        // Create chi2pid histograms for pions, kaons, and proton (EB PID: 211, 321, 2212) under pion hypothesis  after

        TH1F* chi2pidallafter = new TH1F(TString::Format("chi2pidallafter_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                         100, -15, 15);

        // Fiil chi2pidall histograms 


        treePionsall->SetBranchAddress("p", &p); 
        treePionsall->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid); 
        for ( Long64_t j = 0; j < treePionsall->GetEntries(); j++){
            treePionsall->GetEntry(j); 
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) chi2pidallafter->Fill(recomputed_chi2pid); 
        }

        // Plot chi2pidall in canvas 
        canvaschi2pidallafter->Clear(); 
        chi2pidallafter-> Draw("HIST"); 
        canvaschi2pidallafter->Update(); 
        canvaschi2pidallafter->Print("output/pdf/kaon/chi2pidforallafter.pdf");

        //--------------------------------//--------------------------------------------------------//-----------------------------//


        // Create chi2pid histograms for +ve charge particles under pion hypothesis  after

        TH1F* chi2pidposafter = new TH1F(TString::Format("chi2pidposafter_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                         100, -15, 15);

        // Fiil chi2pidpos histograms 


        treeEBPosPionAssumed->SetBranchAddress("p", &p); 
        treeEBPosPionAssumed->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid); 
        for ( Long64_t j = 0; j < treeEBPosPionAssumed->GetEntries(); j++){
            treeEBPosPionAssumed->GetEntry(j); 
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) chi2pidposafter->Fill(recomputed_chi2pid); 
        }

        // Plot chi2pidall in canvas 
        canvaschi2pidposafter->Clear(); 
        chi2pidposafter-> Draw("HIST"); 
        canvaschi2pidposafter->Update(); 
        canvaschi2pidposafter->Print("output/pdf/kaon/chi2pidpos_after.pdf");


        //--------------------------------//--------------------------------------------------------//-----------------------------//


        // Create chi2pid histograms (after cut)
        TH1F* chi2PionsAfter = new TH1F(TString::Format("chi2_pions_after_%d", i),
                                        TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                        100, chi2HistRanges[i].first, chi2HistRanges[i].second);
        TH1F* chi2KaonsAfter = new TH1F(TString::Format("chi2_kaons_after_%d", i),
                                        TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                        100, chi2HistRanges[i].first, chi2HistRanges[i].second);

        // Fill chi2pid histograms (after cut) with momentum and ±3σ cut
         
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);
        for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) {
                chi2PionsAfter->Fill(orig_chi2pid);
            }
        }

        treeKaons->SetBranchAddress("p", &p);
        treeKaons->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);
        for (Long64_t j = 0; j < treeKaons->GetEntries(); j++) {
            treeKaons->GetEntry(j);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) {
                chi2KaonsAfter->Fill(recomputed_chi2pid);
            }
        }

        // Set colors for chi2pid histograms (after cut)
        chi2PionsAfter->SetLineColor(kBlue);
        chi2KaonsAfter->SetLineColor(kGreen);

        // Fit chi2pid histograms (after cut) with a Gaussian
        TF1* fitPionAfter = new TF1("gausFit", "gaus", chi2FitRanges[i].first, chi2FitRanges[i].second);
        fitPionAfter->SetParameters(pionAmp, pionMeanInit, pionSigmaInit);
        chi2PionsAfter->Fit(fitPionAfter, "R", "", chi2FitRanges[i].first, chi2FitRanges[i].second);
        
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
        chi2KaonsAfter->Draw("HIST");
        if (chi2KaonsAfter->GetFunction("gausFit")) {
            chi2KaonsAfter->GetFunction("gausFit")->SetLineColor(kGreen);
            chi2KaonsAfter->GetFunction("gausFit")->Draw("SAME");
        }
        TLegend* legKaonsAfter = new TLegend(0.1, 0.75, 0.3, 0.9);
        legKaonsAfter->AddEntry(chi2KaonsAfter, "Kaons (Pion Hypothesis)", "l");
        legKaonsAfter->SetBorderSize(0);
        legKaonsAfter->SetFillStyle(0);
        legKaonsAfter->SetTextSize(0.025);
        legKaonsAfter->Draw();
        canvasChi2After->Update();
        canvasChi2After->Print("output/pdf/kaon/chi2pid_fits_linear_after.pdf");
        delete legPionsAfter;
        delete legKaonsAfter;

        // Create beta histograms (after cut)
        TH1F* betaPionsAfter = new TH1F(TString::Format("beta_pions_after_%d", i),
                                        TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                        70, betaHistRanges[i].first, betaHistRanges[i].second);
        TH1F* betaKaonsAfter = new TH1F(TString::Format("beta_kaons_after_%d", i),
                                        TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                        70, betaHistRanges[i].first, betaHistRanges[i].second);
        betaPionsAfter->Sumw2();
        betaKaonsAfter->Sumw2();

        // Fill beta histograms (after cut) with momentum and ±3σ cut
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("beta", &beta);
        treePions->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);
        for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) {
                betaPionsAfter->Fill(beta);
            }
        }
        
        treeKaons->SetBranchAddress("p", &p);
        treeKaons->SetBranchAddress("beta", &beta);
        treeKaons->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);
        for (Long64_t j = 0; j < treeKaons->GetEntries(); j++) {
            treeKaons->GetEntry(j);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) {
                betaKaonsAfter->Fill(beta);
            }
        }

        // Fit beta histograms (after cut) for contamination calculation
        double pMeanAfter = 0, pSigmaAfter = 0, pConstAfter = 0, kMeanAfter = 0, kSigmaAfter = 0, kConstAfter = 0;
        if (betaPionsAfter->GetEntries() >= 10) {
            TF1* gausPionAfter = new TF1(TString::Format("gaus_%s", betaPionsAfter->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            if (i == 9) { // Bin [3.7-4.0) GeV/c
                gausPionAfter->SetParameters(betaPionsAfter->GetMaximum(), 0.996, 0.004);
                betaPionsAfter->Fit(gausPionAfter, "R", "", 0.994, 1.006);
            } else {
                gausPionAfter->SetParameters(betaPionsAfter->GetMaximum(), 1.0, 0.005);
                betaPionsAfter->Fit(gausPionAfter, "R", "", betaFitRanges[i].first, betaFitRanges[i].second);
            }
            pMeanAfter = gausPionAfter->GetParameter(1);
            pSigmaAfter = gausPionAfter->GetParameter(2);
            pConstAfter = gausPionAfter->GetParameter(0);
            delete gausPionAfter;
        }
        if (betaKaonsAfter->GetEntries() >= 10) {
            TF1* gausKaonAfter = new TF1(TString::Format("gaus_%s", betaKaonsAfter->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            if (i == 9) { // Bin [3.7-4.0) GeV/c
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.991, 0.996);
            } else if (i == 8) { // Bin [3.4-3.7) GeV/c
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.991, 0.995);
            } else if (i == 10) { // Bin [4.0-4.3) GeV/c
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.991, 0.996);
            } else if (i == 11) { // Bin [4.3-4.6) GeV/c
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.996, 0.012);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.993, 0.9965);
            } else if (i == 12) { // Bin [4.6-4.9) GeV/c
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.012);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.993, 0.9965);
            } else if (i == 13) { // Bin [4.9-5.2) GeV/c
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.012);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.991, 0.995);
            } else if (i == 14) { // Bin [5.2-5.5) GeV/c
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.991, 0.996);
            } else if (i == 15) { // Bin [5.5-5.8) GeV/c
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.991, 0.996);
            } else if (i == 16) { // Bin [5.8-6.1) GeV/c
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.991, 0.997);
            } else if (i == 17) { // Bin [6.1-6.4) GeV/c
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.991, 0.997);
            } else if (i == 7) { // Bin [6.4-6.7) GeV/c
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.984, 0.992);
            } else if (i == 6) { // Bin [2.5-2.8) GeV/c
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.988, 0.992);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.984, 0.992);
            } else {
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.99, 0.005);
                betaKaonsAfter->Fit(gausKaonAfter, "R");
            }
            kMeanAfter = gausKaonAfter->GetParameter(1);
            kSigmaAfter = gausKaonAfter->GetParameter(2);
            kConstAfter = gausKaonAfter->GetParameter(0);
            delete gausKaonAfter;
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
        if (betaKaonsAfter->GetEntries() >= 10) {
            betaKaonsAfter->SetMarkerStyle(20);
            betaKaonsAfter->SetMarkerSize(1);
            betaKaonsAfter->SetMarkerColor(kGreen);
            betaKaonsAfter->Draw(betaPionsAfter->GetEntries() >= 10 ? "P SAME" : "P");
            yMax = max(yMax, betaKaonsAfter->GetMaximum());
            if (auto* fit = betaKaonsAfter->GetFunction(TString::Format("gaus_%s", betaKaonsAfter->GetName()))) {
                fit->SetLineColor(kGreen);
                fit->SetLineStyle(kDashed);
                fit->Draw("SAME");
                if (i == 6) {fit->SetRange(0.97, 1.01);}
                else {fit->SetRange(0.95, 1.03);}
            }
        }
        yMax *= 1.2;
        if (betaPionsAfter->GetEntries() >= 10 || betaKaonsAfter->GetEntries() >= 10) {
            (betaPionsAfter->GetEntries() >= 10 ? betaPionsAfter : betaKaonsAfter)->GetYaxis()->SetRangeUser(yMin, yMax);
        }
        TLegend* legAfter = new TLegend(0.15, 0.75, 0.35, 0.88);
        if (betaPionsAfter->GetEntries() >= 10) legAfter->AddEntry(betaPionsAfter, "Pions", "p");
        if (betaKaonsAfter->GetEntries() >= 10) legAfter->AddEntry(betaKaonsAfter, "Kaons (Misidentified)", "p");
        legAfter->SetTextSize(0.03);
        legAfter->Draw();
        canvasBetaAfter->SetGrid();
        canvasBetaAfter->Update();
        canvasBetaAfter->Print("output/pdf/kaon/beta_after_cut_linear.pdf");
        delete legAfter;

        // Calculate the contamination after the chi2pid cut
        double c1After, c2After;
        map<string, double> pFitAfter = {{"gaus_mean", pMeanAfter}, {"gaus_sigma", pSigmaAfter}, {"gaus_constant", pConstAfter}};
        map<string, double> kFitAfter = {{"gaus_mean", kMeanAfter}, {"gaus_sigma", kSigmaAfter}, {"gaus_constant", kConstAfter}};
        double kaonToPionRatioAfter = calculateOverlapContamination(pFitAfter, kFitAfter, pionLeftAfter[i], pionRightAfter[i], kaonRightAfter[i], c1After, c2After);
        double kaonLeftAfter = -1.0;
        if (c1After >= kMeanAfter && c1After <= pMeanAfter) {
            kaonLeftAfter = c1After;
        } else if (c2After >= kMeanAfter && c2After <= pMeanAfter) {
            kaonLeftAfter = c2After;
        } else if (kaonToPionRatioAfter >= 0) {
            kaonLeftAfter = (abs(c1After - kMeanAfter) < abs(c2After - kMeanAfter)) ? c1After : c2After;
            cout << "No intersection between kMeanAfter = " << kMeanAfter << " and pMeanAfter = " << pMeanAfter 
                 << ", selecting point closest to kaon peak: " << kaonLeftAfter << endl;
        }

        // Save the contamination (after cut) to CSV
        if (kaonToPionRatioAfter < 0 || kaonLeftAfter == -1.0) {
            cout << "Bin [" << pLow << "-" << pHigh << "): Contamination (after cut) calculation failed\n";
            csvFileAfter << pLow << "-" << pHigh << ",N/A,N/A\n";
        } else {
            cout << "Bin [" << pLow << "-" << pHigh << "): After cut - Pion range: [" << pionLeftAfter[i] << ", " << pionRightAfter[i] 
                 << "], Kaon range: [" << kaonLeftAfter << ", " << kaonRightAfter[i] << "], Contamination: " 
                 << kaonToPionRatioAfter << "%\n";
            csvFileAfter << pLow << "-" << pHigh << "," << kaonLeftAfter << "," << kaonToPionRatioAfter << "\n";
        }

        // Clean up histograms and fits for this bin
        delete chi2PionsBefore;
        delete chi2KaonsBefore;
        delete betaPionsBefore;
        delete betaKaonsBefore;
        delete chi2PionsAfter;
        delete chi2KaonsAfter;
        delete betaPionsAfter;
        delete betaKaonsAfter;
        delete chi2pidall; 
        delete chi2pidallafter; 
        delete chi2pidposbefore; 
        delete chi2pidposafter; 
        delete fitPion;
       
        delete fitPionAfter;
        
    }

    // Step 10: Close PDF and CSV files
    canvaschi2pidall->Print("output/pdf/kaon/chi2pidforall.pdf]"); 
    canvaschi2pidallafter->Print("output/pdf/kaon/chi2pidforallafter.pdf]"); 
    canvaschi2pidposbefore->Print("output/pdf/kaon/chi2pidpos_before.pdf]"); 
    canvaschi2pidposafter->Print("output/pdf/kaon/chi2pidpos_after.pdf]"); 
    canvasChi2Before->Print("output/pdf/kaon/chi2pid_fits_linear_before.pdf]");
    canvasChi2After->Print("output/pdf/kaon/chi2pid_fits_linear_after.pdf]");
    canvasBetaBefore->Print("output/pdf/kaon/beta_before_cut_linear.pdf]");
    canvasBetaAfter->Print("output/pdf/kaon/beta_after_cut_linear.pdf]");
    csvFileBefore.close();
    csvFileAfter.close();
    csvFileCuts.close();

    // Step 11: Clean up
    delete canvasChi2Before;
    delete canvasChi2After;
    delete canvasBetaBefore;
    delete canvasBetaAfter;
    delete canvaschi2pidposbefore; 
    delete canvaschi2pidposafter; 
    file->Close();
    delete file;

    return 0;
}  

 */


#include "include/ContaminationUtils.h"
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLegend.h>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

int main() {
    // Step 1: Open the ROOT file and load the trees for pions and kaons
    TFile* file = new TFile("../Skim/pkptreeCxC_9_test_modified.root");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file pkptreeCxC_9_test_modified.root\n";
        return 1;
    }
    TTree* treePions = (TTree*)file->Get("EB_pid_pions");
    TTree* treeKaons = (TTree*)file->Get("EB_pid_kaons");
    TTree* treePionsall = (TTree*)file->Get("EB_all_pion_assumed"); 
    TTree* treeEBPosPionAssumed = (TTree*)file->Get("EB_pos_pion_assumed"); 
    if (!treePions || !treeKaons || !treePionsall || !treeEBPosPionAssumed) {
        cerr << "Error: Cannot load trees\n";
        file->Close();
        delete file;
        return 1;
    }

    // Step 2: Create output directories for saving PDFs and CSVs
    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/pdf/kaon", kTRUE);
    gSystem->mkdir("output/csv/kaon", kTRUE);

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
            chi2HistRanges[i] = {-5.0, 5.0};
 chi2FitRanges[i] = {-1.3, 1.3};
            chi2FitParams[i] = {1000.0, 0.0, 1.0};
        } else if (pLow < 3.0) {
            chi2HistRanges[i] = {-4.0, 4.0};
            chi2FitRanges[i] = {-1.5, 1.5};
            chi2FitParams[i] = {800.0, 0.1, 0.8};
        } else {
            chi2HistRanges[i] = {-3.0, 3.0};
            chi2FitRanges[i] = {-1, 1};
            chi2FitParams[i] = {500.0, -0.1, 0.5};
        }
    }

    // Step 5: Set beta histogram ranges and fit ranges for each bin
    vector<pair<double, double>> betaHistRanges(nBins);
    vector<pair<double, double>> betaFitRanges(nBins);
    for (int i = 0; i < nBins; i++) {
        betaHistRanges[i] = {0.95, 1.03};
        betaFitRanges[i] = {0.95, 1.02};
    }

    // Step 6: Set beta cuts for pions and kaons
    vector<double> pionLeftBefore(nBins), pionRightBefore(nBins), kaonRightBefore(nBins);
    vector<double> pionLeftAfter(nBins), pionRightAfter(nBins), kaonRightAfter(nBins);
    for (int i = 0; i < nBins; i++) {
        if (i == 0) { // Bin [1.0-1.3)
            pionLeftBefore[i] = 0.9880; pionRightBefore[i] = 1.0120; kaonRightBefore[i] = 0.0;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  kaonRightAfter[i] = 0.0;
        } else if (i == 1) { // Bin [1.3-1.6)
            pionLeftBefore[i] = 0.9755; pionRightBefore[i] = 1.0135; kaonRightBefore[i] = 0.9721;
            pionLeftAfter[i] = 0.9755;  pionRightAfter[i] = 1.0135;  kaonRightAfter[i] = 0.9721;
        } else if (i == 2) { // Bin [1.6-1.9)
            pionLeftBefore[i] = 0.98;   pionRightBefore[i] = 1.0140; kaonRightBefore[i] = 0.986;
            pionLeftAfter[i] = 0.98;    pionRightAfter[i] = 1.0140;  kaonRightAfter[i] = 0.986;
        } else if (i == 3) { // Bin [1.9-2.2)
            pionLeftBefore[i] = 0.9815; pionRightBefore[i] = 1.0145; kaonRightBefore[i] = 0.992;
            pionLeftAfter[i] = 0.9815;  pionRightAfter[i] = 1.0145;  kaonRightAfter[i] = 0.992;
        } else if (i == 4) { // Bin [2.2-2.5)
            pionLeftBefore[i] = 0.9880; pionRightBefore[i] = 1.0141; kaonRightBefore[i] = 0.9971;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0141;  kaonRightAfter[i] = 0.9971;
        } else if (i == 5) { // Bin [2.5-2.8)
            pionLeftBefore[i] = 0.9848; pionRightBefore[i] = 1.013;  kaonRightBefore[i] = 1.0;
            pionLeftAfter[i] = 0.9848;  pionRightAfter[i] = 1.013;   kaonRightAfter[i] = 1.0;
        } else if (i == 6) { // Bin [2.8-3.1)
            pionLeftBefore[i] = 0.986;  pionRightBefore[i] = 1.0128; kaonRightBefore[i] = 1.003;
            pionLeftAfter[i] = 0.986;   pionRightAfter[i] = 1.0128;  kaonRightAfter[i] = 1.003;
        } else if (i == 7) { // Bin [3.1-3.4)
            pionLeftBefore[i] = 0.9868; pionRightBefore[i] = 1.0121; kaonRightBefore[i] = 1.0048;
            pionLeftAfter[i] = 0.9868;  pionRightAfter[i] = 1.0121;  kaonRightAfter[i] = 1.0048;
        } else if (i == 8) { // Bin [3.4-3.7)
            pionLeftBefore[i] = 0.9879; pionRightBefore[i] = 1.0120; kaonRightBefore[i] = 1.0049;
            pionLeftAfter[i] = 0.9879;  pionRightAfter[i] = 1.0120;  kaonRightAfter[i] = 1.0049;
        } else if (i == 9) { // Bin [3.7-4.0)
            pionLeftBefore[i] = 0.9880; pionRightBefore[i] = 1.0120; kaonRightBefore[i] = 1.012;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  kaonRightAfter[i] = 1.012;
        } else {
            cerr << "Warning: Using default beta cuts for bin [" << pBins[i] << "-" << pBins[i+1] << ")\n";
            pionLeftBefore[i] = 0.9880; pionRightBefore[i] = 1.0120; kaonRightBefore[i] = 1.0;
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  kaonRightAfter[i] = 1.0;
        }
    }

    // Step 7: Create canvases for beta plots only (comment out chi2pid canvases)
    
    TCanvas* canvasChi2Before = new TCanvas("canvasChi2Before", "Chi2pid Fits (Before)", 1900, 900);
    canvasChi2Before->SetMargin(0.1, 0.05, 0.15, 0.1);
    TCanvas* canvasChi2After = new TCanvas("canvasChi2After", "Chi2pid Fits (After)", 1900, 900);
    canvasChi2After->SetMargin(0.1, 0.05, 0.15, 0.1);
    TCanvas* canvasBetaBefore = new TCanvas("canvasBetaBefore", "Beta Plots (Before)", 1200, 800);
    canvasBetaBefore->SetMargin(0.1, 0.05, 0.15, 0.1);
    TCanvas* canvasBetaAfter = new TCanvas("canvasBetaAfter", "Beta Plots (After)", 800, 600);
    canvasBetaAfter->SetMargin(0.1, 0.05, 0.15, 0.1);
    /*
    TCanvas* canvaschi2pidall = new TCanvas("canvaschi2pidall", "chi2pid Plot", 1200, 800); 
    TCanvas* canvaschi2pidallafter = new TCanvas("canvaschi2pidallafter", "chi2pid Plot (After)", 1200, 800); 
    TCanvas* canvaschi2pidposbefore = new TCanvas("canvaschi2pidpos", "chi2pid Plot +ve Charge Particles (Before)", 1200, 800); 
    TCanvas* canvaschi2pidposafter = new TCanvas("canvaschi2pidposafter", "chi2pid Plot +ve Charge Particles (After)", 1200, 800);
    */

    // Comment out chi2pid PDF openings
    
    //canvaschi2pidall->Print("output/pdf/kaon/chi2pidforall.pdf[");
    //canvaschi2pidallafter->Print("output/pdf/kaon/chi2pidforallafter.pdf[");
    canvasChi2Before->Print("output/pdf/kaon/chi2pid_fits_linear_before.pdf[");
    canvasChi2After->Print("output/pdf/kaon/chi2pid_fits_linear_after.pdf[");
    //canvaschi2pidposbefore->Print("output/pdf/kaon/chi2pidpos_before.pdf["); 
    //canvaschi2pidposafter->Print("output/pdf/kaon/chi2pidpos_after.pdf["); 
    
    canvasBetaBefore->Print("output/pdf/kaon/beta_before_cut_linear.pdf[");
    canvasBetaAfter->Print("output/pdf/kaon/beta_after_cut_linear.pdf[");

    // Step 8: Open CSV files to save results
    ofstream csvFileBefore("output/csv/kaon/kaon_to_pion_ratio_before_chi2pid.csv");
    csvFileBefore << "Momentum Bin (GeV/c),Kaon Left,Contamination (%)\n";
    ofstream csvFileAfter("output/csv/kaon/contamination_after_chi2pid.csv");
    csvFileAfter << "Momentum Bin (GeV/c),Kaon Left,Kaon-to-Pion Ratio (%)\n";
    ofstream csvFileCuts("output/csv/kaon/chi2pid_cuts.csv");
    csvFileCuts << "Momentum Bin (GeV/c),Pion Mean,Pion Sigma,Chi2 Min (Mean-5σ),Chi2 Max (Mean+5σ)\n";

    // Step 9: Loop over each momentum bin to process data
    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];

        // --- Before Chi2pid Cut ---
        // Create chi2pid histograms for pions and kaons (before cut)
        TH1F* chi2PionsBefore = new TH1F(TString::Format("chi2_pions_before_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                         100, -15, 15);
        TH1F* chi2KaonsBefore = new TH1F(TString::Format("chi2_kaons_before_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                         100, -15, 15);
        /* TH1F* chi2pidall = new TH1F(TString::Format("chi2pidall_%d", i),
                                    TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                    100, -15, 15);
        TH1F* chi2pidposbefore = new TH1F(TString::Format("chi2pidposbefore_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                         100, -15, 15); */
 
        // Fill chi2pid histograms
        float p, chi2pid, beta, orig_chi2pid, recomputed_chi2pid;
        /* treePionsall->SetBranchAddress("p", &p); 
        treePionsall->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid); 
        for (Long64_t j = 0; j < treePionsall->GetEntries(); j++) {
            treePionsall->GetEntry(j); 
            if (p >= pLow && p < pHigh) chi2pidall->Fill(recomputed_chi2pid); 
        }
        treeEBPosPionAssumed->SetBranchAddress("p", &p); 
        treeEBPosPionAssumed->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid); 
        for (Long64_t j = 0; j < treeEBPosPionAssumed->GetEntries(); j++) {
            treeEBPosPionAssumed->GetEntry(j); 
            if (p >= pLow && p < pHigh) chi2pidposbefore->Fill(recomputed_chi2pid); 
        } */
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("orig_chi2pid", &orig_chi2pid);
        for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh) chi2PionsBefore->Fill(orig_chi2pid);
        }
        treeKaons->SetBranchAddress("p", &p);
        treeKaons->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);
        for (Long64_t j = 0; j < treeKaons->GetEntries(); j++) {
            treeKaons->GetEntry(j);
            if (p >= pLow && p < pHigh) chi2KaonsBefore->Fill(recomputed_chi2pid);
        } 

        // Comment out chi2pid plotting
        /*
        canvaschi2pidall->Clear(); 
        chi2pidall->Draw("HIST"); 
        canvaschi2pidall->Update(); 
        canvaschi2pidall->Print("output/pdf/kaon/chi2pidforall.pdf");
        canvaschi2pidposbefore->Clear(); 
        chi2pidposbefore->Draw("HIST"); 
        canvaschi2pidposbefore->Update(); 
        canvaschi2pidposbefore->Print("output/pdf/kaon/chi2pidpos_before.pdf");
        */

        // Fit chi2pid histograms with a Gaussian to determine the ±5σ cut
        double pionMeanBefore, pionSigmaBefore, pionChi2Min, pionChi2Max;
        auto [pionAmp, pionMeanInit, pionSigmaInit] = chi2FitParams[i];
        TF1* fitPion = new TF1("gausFit", "gaus", chi2FitRanges[i].first, chi2FitRanges[i].second);
        fitPion->SetParameters(pionAmp, pionMeanInit, pionSigmaInit);
        chi2PionsBefore->Fit(fitPion, "R", "", chi2FitRanges[i].first, chi2FitRanges[i].second);
        pionMeanBefore = fitPion->GetParameter(1);
        pionSigmaBefore = fitPion->GetParameter(2);
        pionChi2Min = pionMeanBefore - 5 * pionSigmaBefore;
        pionChi2Max = pionMeanBefore + 5 * pionSigmaBefore;
        
        cout << "Bin [" << pLow << "-" << pHigh << "): Pions chi2 range (before): [" 
             << pionChi2Min << ", " << pionChi2Max << "]\n";

        // Save chi2pid fit parameters to CSV
        csvFileCuts << pLow << "-" << pHigh << "," << pionMeanBefore << "," << pionSigmaBefore << "," 
                    << pionChi2Min << "," << pionChi2Max << "\n";

        // Comment out chi2pid plotting (before cut)
        
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
        chi2KaonsBefore->Draw("HIST");
        if (chi2KaonsBefore->GetFunction("gausFit")) {
            chi2KaonsBefore->GetFunction("gausFit")->SetLineColor(kGreen);
            chi2KaonsBefore->GetFunction("gausFit")->Draw("SAME");
        }
        TLegend* legKaons = new TLegend(0.1, 0.75, 0.3, 0.9);
        legKaons->AddEntry(chi2KaonsBefore, "Kaons", "l");
        legKaons->SetBorderSize(0);
        legKaons->SetFillStyle(0);
        legKaons->SetTextSize(0.025);
        legKaons->Draw();
        canvasChi2Before->Update();
        canvasChi2Before->Print("output/pdf/kaon/chi2pid_fits_linear_before.pdf");
        delete legPions;
        delete legKaons;
        

        // Create beta histograms (before cut) with increased bins for large entries
        TH1F* betaPionsBefore = new TH1F(TString::Format("beta_pions_before_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                         100, betaHistRanges[i].first, betaHistRanges[i].second);
        TH1F* betaKaonsBefore = new TH1F(TString::Format("beta_kaons_before_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                         100, betaHistRanges[i].first, betaHistRanges[i].second);
        betaPionsBefore->Sumw2();
        betaKaonsBefore->Sumw2();

        // Fill beta histograms (before cut)
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("beta", &beta);
        for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh) {
                betaPionsBefore->Fill(beta);
            }
        }
        treeKaons->SetBranchAddress("p", &p);
        treeKaons->SetBranchAddress("beta", &beta);
        for (Long64_t j = 0; j < treeKaons->GetEntries(); j++) {
            treeKaons->GetEntry(j);
            if (p >= pLow && p < pHigh) {
                betaKaonsBefore->Fill(beta);
            }
        }

        // Fit beta histograms (before cut)
        double pMeanBefore = 0, pSigmaBefore = 0, pConstBefore = 0, kMeanBefore = 0, kSigmaBefore = 0, kConstBefore = 0;
        if (betaPionsBefore->GetEntries() >= 10) {
            TF1* gausPionBefore = new TF1(TString::Format("gaus_%s", betaPionsBefore->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            if (i == 9) { // Bin [3.7-4.0) GeV/c
                gausPionBefore->SetParameters(betaPionsBefore->GetMaximum(), 0.996, 0.004);
                betaPionsBefore->Fit(gausPionBefore, "R", "", 0.994, 1.006);
            } else {
                gausPionBefore->SetParameters(betaPionsBefore->GetMaximum(), 1.0, 0.005);
                betaPionsBefore->Fit(gausPionBefore, "R", "", betaFitRanges[i].first, betaFitRanges[i].second);
            }
            pMeanBefore = gausPionBefore->GetParameter(1);
            pSigmaBefore = gausPionBefore->GetParameter(2);
            pConstBefore = gausPionBefore->GetParameter(0);
            delete gausPionBefore;
        }
        if (betaKaonsBefore->GetEntries() >= 10) {
            TF1* gausKaonBefore = new TF1(TString::Format("gaus_%s", betaKaonsBefore->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            if (i == 9) {
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.986, 0.9945);
            } else if (i == 8) {
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.988, 0.9945);
            } else if (i == 10) {
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.987, 0.995);
            } else if (i == 11) {
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.996, 0.012);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.992, 0.996);
            } else if (i == 12) {
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.012);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.991, 0.996);
            } else if (i == 13) {
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.012);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.991, 0.995);
            } else if (i == 14) {
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.991, 0.996);
            } else if (i == 15) {
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.991, 0.996);
            } else if (i == 16) {
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.991, 0.997);
            } else if (i == 17) {
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.991, 0.997);
            } else if (i == 18) {
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.982, 0.997);
            } else if (i == 5) {
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.974, 0.992);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.974, 0.992);
            } else if (i == 6) {
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.986, 0.992);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.98, 0.993);
            } else if (i == 7) {
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.989, 0.992);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.982, 0.993);
            }else {
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.99, 0.005);
                betaKaonsBefore->Fit(gausKaonBefore, "R");
            }
            kMeanBefore = gausKaonBefore->GetParameter(1);
            kSigmaBefore = gausKaonBefore->GetParameter(2);
            kConstBefore = gausKaonBefore->GetParameter(0);
            delete gausKaonBefore;
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
        if (betaKaonsBefore->GetEntries() >= 10) {
            betaKaonsBefore->SetMarkerStyle(20);
            betaKaonsBefore->SetMarkerSize(1);
            betaKaonsBefore->SetMarkerColor(kGreen);
            betaKaonsBefore->Draw(betaPionsBefore->GetEntries() >= 10 ? "P SAME" : "P");
            yMax = max(yMax, betaKaonsBefore->GetMaximum());
            if (auto* fit = betaKaonsBefore->GetFunction(TString::Format("gaus_%s", betaKaonsBefore->GetName()))) {
                fit->SetLineColor(kGreen);
                fit->SetLineStyle(kDashed);
                fit->Draw("SAME");
                fit->SetRange(0.95, 1.03);
            }
        }
        yMax *= 1.2;
        if (betaPionsBefore->GetEntries() >= 10 || betaKaonsBefore->GetEntries() >= 10) {
            (betaPionsBefore->GetEntries() >= 10 ? betaPionsBefore : betaKaonsBefore)->GetYaxis()->SetRangeUser(yMin, yMax);
        }
        TLegend* leg = new TLegend(0.15, 0.75, 0.35, 0.88);
        if (betaPionsBefore->GetEntries() >= 10) leg->AddEntry(betaPionsBefore, "Pions", "p");
        if (betaKaonsBefore->GetEntries() >= 10) leg->AddEntry(betaKaonsBefore, "Kaons", "p");
        leg->SetTextSize(0.03);
        leg->Draw();
        canvasBetaBefore->SetGrid();
        canvasBetaBefore->Update();
        canvasBetaBefore->Print("output/pdf/kaon/beta_before_cut_linear.pdf");
        delete leg;

        // Calculate kaon-to-pion ratio before chi2pid cut
        double c1Before, c2Before;
        map<string, double> pFitBefore = {{"gaus_mean", pMeanBefore}, {"gaus_sigma", pSigmaBefore}, {"gaus_constant", pConstBefore}};
        map<string, double> kFitBefore = {{"gaus_mean", kMeanBefore}, {"gaus_sigma", kSigmaBefore}, {"gaus_constant", kConstBefore}};
        double contaminationBefore = calculateOverlapContamination(pFitBefore, kFitBefore, pionLeftBefore[i], pionRightBefore[i], kaonRightBefore[i], c1Before, c2Before);
        double kaonLeftBefore = -1.0;
        if (c1Before >= kMeanBefore && c1Before <= pMeanBefore) {
            kaonLeftBefore = c1Before;
        } else if (c2Before >= kMeanBefore && c2Before <= pMeanBefore) {
            kaonLeftBefore = c2Before;
        } else if (contaminationBefore >= 0) {
            kaonLeftBefore = (abs(c1Before - kMeanBefore) < abs(c2Before - kMeanBefore)) ? c1Before : c2Before;
            cout << "No intersection between kMeanBefore = " << kMeanBefore << " and pMeanBefore = " << pMeanBefore 
                 << ", selecting point closest to kaon peak: " << kaonLeftBefore << endl;
        }

        // Save kaon-to-pion ratio (before cut) to CSV
        if (contaminationBefore < 0 || kaonLeftBefore == -1.0) {
            cout << "Bin [" << pLow << "-" << pHigh << "): Kaon-to-pion ratio (before cut) calculation failed\n";
            csvFileBefore << pLow << "-" << pHigh << ",N/A,N/A\n";
        } else {
            cout << "Bin [" << pLow << "-" << pHigh << "): Before cut - Pion range: [" << pionLeftBefore[i] << ", " << pionRightBefore[i] 
                 << "], Kaon range: [" << kaonLeftBefore << ", " << kaonRightBefore[i] << "], Kaon-to-Pion Ratio: " 
                 << contaminationBefore << "%\n";
            csvFileBefore << pLow << "-" << pHigh << "," << kaonLeftBefore << "," << contaminationBefore << "\n";
        }

        // --- After Chi2pid Cut (±5σ) ---
        /*  TH1F* chi2pidallafter = new TH1F(TString::Format("chi2pidallafter_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                         100, -15, 15);
        TH1F* chi2pidposafter = new TH1F(TString::Format("chi2pidposafter_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                         100, -15, 15); */
        TH1F* chi2PionsAfter = new TH1F(TString::Format("chi2_pions_after_%d", i),
                                        TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                        100, -15, 15);
        TH1F* chi2KaonsAfter = new TH1F(TString::Format("chi2_kaons_after_%d", i),
                                        TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                        100, -15, 15); 

        // Fill chi2pid histograms (after cut)
        /* treePionsall->SetBranchAddress("p", &p); 
        treePionsall->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid); 
        for (Long64_t j = 0; j < treePionsall->GetEntries(); j++) {
            treePionsall->GetEntry(j); 
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) chi2pidallafter->Fill(recomputed_chi2pid); 
        }
        treeEBPosPionAssumed->SetBranchAddress("p", &p); 
        treeEBPosPionAssumed->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid); 
        for (Long64_t j = 0; j < treeEBPosPionAssumed->GetEntries(); j++) {
            treeEBPosPionAssumed->GetEntry(j); 
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) chi2pidposafter->Fill(recomputed_chi2pid); 
        } */
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("orig_chi2pid", &orig_chi2pid);
        treePions->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);
        for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) {
                chi2PionsAfter->Fill(orig_chi2pid);
            }
        }
        treeKaons->SetBranchAddress("p", &p);
        treeKaons->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);
        for (Long64_t j = 0; j < treeKaons->GetEntries(); j++) {
            treeKaons->GetEntry(j);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) {
                chi2KaonsAfter->Fill(recomputed_chi2pid);
            }
        } 

        // Comment out chi2pid plotting (after cut)
        
        /* canvaschi2pidallafter->Clear(); 
        chi2pidallafter->Draw("HIST"); 
        canvaschi2pidallafter->Update(); 
        canvaschi2pidallafter->Print("output/pdf/kaon/chi2pidforallafter.pdf");
        canvaschi2pidposafter->Clear(); 
        chi2pidposafter->Draw("HIST"); 
        canvaschi2pidposafter->Update(); 
        canvaschi2pidposafter->Print("output/pdf/kaon/chi2pidpos_after.pdf"); */
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
        chi2KaonsAfter->Draw("HIST");
        if (chi2KaonsAfter->GetFunction("gausFit")) {
            chi2KaonsAfter->GetFunction("gausFit")->SetLineColor(kGreen);
            chi2KaonsAfter->GetFunction("gausFit")->Draw("SAME");
        }
        TLegend* legKaonsAfter = new TLegend(0.1, 0.75, 0.3, 0.9);
        legKaonsAfter->AddEntry(chi2KaonsAfter, "Kaons (Pion Hypothesis)", "l");
        legKaonsAfter->SetBorderSize(0);
        legKaonsAfter->SetFillStyle(0);
        legKaonsAfter->SetTextSize(0.025);
        legKaonsAfter->Draw();
        canvasChi2After->Update();
        canvasChi2After->Print("output/pdf/kaon/chi2pid_fits_linear_after.pdf");
        delete legPionsAfter;
        delete legKaonsAfter;
        

        // Create beta histograms (after cut) with increased bins
        TH1F* betaPionsAfter = new TH1F(TString::Format("beta_pions_after_%d", i),
                                        TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                        100, betaHistRanges[i].first, betaHistRanges[i].second);
        TH1F* betaKaonsAfter = new TH1F(TString::Format("beta_kaons_after_%d", i),
                                        TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                        100, betaHistRanges[i].first, betaHistRanges[i].second);
        betaPionsAfter->Sumw2();
        betaKaonsAfter->Sumw2();

        // Fill beta histograms (after cut) with ±5σ cut
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("beta", &beta);
        treePions->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);
        for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) {
                betaPionsAfter->Fill(beta);
            }
        }
        treeKaons->SetBranchAddress("p", &p);
        treeKaons->SetBranchAddress("beta", &beta);
        treeKaons->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);
        for (Long64_t j = 0; j < treeKaons->GetEntries(); j++) {
            treeKaons->GetEntry(j);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= pionChi2Min && recomputed_chi2pid <= pionChi2Max) {
                betaKaonsAfter->Fill(beta);
            }
        }

        // Fit beta histograms (after cut)
        double pMeanAfter = 0, pSigmaAfter = 0, pConstAfter = 0, kMeanAfter = 0, kSigmaAfter = 0, kConstAfter = 0;
        if (betaPionsAfter->GetEntries() >= 10) {
            TF1* gausPionAfter = new TF1(TString::Format("gaus_%s", betaPionsAfter->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            if (i == 9) {
                gausPionAfter->SetParameters(betaPionsAfter->GetMaximum(), 0.996, 0.004);
                betaPionsAfter->Fit(gausPionAfter, "R", "", 0.994, 1.006);
            } else {
                gausPionAfter->SetParameters(betaPionsAfter->GetMaximum(), 1.0, 0.005);
                betaPionsAfter->Fit(gausPionAfter, "R", "", betaFitRanges[i].first, betaFitRanges[i].second);
            }
            pMeanAfter = gausPionAfter->GetParameter(1);
            pSigmaAfter = gausPionAfter->GetParameter(2);
            pConstAfter = gausPionAfter->GetParameter(0);
            delete gausPionAfter;
        }
        if (betaKaonsAfter->GetEntries() >= 10) {
            TF1* gausKaonAfter = new TF1(TString::Format("gaus_%s", betaKaonsAfter->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            if (i == 9) {
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.993, 0.01);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.986, 0.995);
            } else if (i == 8) {
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.993, 0.01);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.987, 0.994);
            } else if (i == 10) {
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.988, 0.996);
            } else if (i == 11) {
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.996, 0.012);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.99, 0.996);
            } else if (i == 12) {
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.012);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.992, 0.997);
            } else if (i == 13) {
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.012);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.99, 0.998);
            } else if (i == 14) {
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.99, 0.998);
            } else if (i == 15) {
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.991, 0.999);
            } else if (i == 16) {
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.992, 1);
            } else if (i == 17) {
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.992, 1);
            }
            else if (i == 18) {
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.994, 1);
            }
            else if (i == 19) {
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.994, 1);
            } else if (i == 7) {
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.99, 0.1);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.984, 0.992);
            } else if (i == 6) {
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.988, 0.01);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.983, 0.993);
            } else if (i == 5) {
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.987, 0.01);
                betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.982, 0.993);
            }else {
                gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.99, 0.005);
                betaKaonsAfter->Fit(gausKaonAfter, "R");
            }
            kMeanAfter = gausKaonAfter->GetParameter(1);
            kSigmaAfter = gausKaonAfter->GetParameter(2);
            kConstAfter = gausKaonAfter->GetParameter(0);
            delete gausKaonAfter;
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
        if (betaKaonsAfter->GetEntries() >= 10) {
            betaKaonsAfter->SetMarkerStyle(20);
            betaKaonsAfter->SetMarkerSize(1);
            betaKaonsAfter->SetMarkerColor(kGreen);
            betaKaonsAfter->Draw(betaPionsAfter->GetEntries() >= 10 ? "P SAME" : "P");
            yMax = max(yMax, betaKaonsAfter->GetMaximum());
            if (auto* fit = betaKaonsAfter->GetFunction(TString::Format("gaus_%s", betaKaonsAfter->GetName()))) {
                fit->SetLineColor(kGreen);
                fit->SetLineStyle(kDashed);
                fit->Draw("SAME");
                if (i == 6) {fit->SetRange(0.97, 1.01);}
                else {fit->SetRange(0.95, 1.03);}
            }
        }
        yMax *= 1.2;
        if (betaPionsAfter->GetEntries() >= 10 || betaKaonsAfter->GetEntries() >= 10) {
            (betaPionsAfter->GetEntries() >= 10 ? betaPionsAfter : betaKaonsAfter)->GetYaxis()->SetRangeUser(yMin, yMax);
        }
        TLegend* legAfter = new TLegend(0.15, 0.75, 0.35, 0.88);
        if (betaPionsAfter->GetEntries() >= 10) legAfter->AddEntry(betaPionsAfter, "Pions", "p");
        if (betaKaonsAfter->GetEntries() >= 10) legAfter->AddEntry(betaKaonsAfter, "Kaons (Misidentified)", "p");
        legAfter->SetTextSize(0.03);
        legAfter->Draw();
        canvasBetaAfter->SetGrid();
        canvasBetaAfter->Update();
        canvasBetaAfter->Print("output/pdf/kaon/beta_after_cut_linear.pdf");
        delete legAfter;

        // Calculate contamination after chi2pid cut
        double c1After, c2After;
        map<string, double> pFitAfter = {{"gaus_mean", pMeanAfter}, {"gaus_sigma", pSigmaAfter}, {"gaus_constant", pConstAfter}};
        map<string, double> kFitAfter = {{"gaus_mean", kMeanAfter}, {"gaus_sigma", kSigmaAfter}, {"gaus_constant", kConstAfter}};
        double kaonToPionRatioAfter = calculateOverlapContamination(pFitAfter, kFitAfter, pionLeftAfter[i], pionRightAfter[i], kaonRightAfter[i], c1After, c2After);
        double kaonLeftAfter = -1.0;
        if (c1After >= kMeanAfter && c1After <= pMeanAfter) {
            kaonLeftAfter = c1After;
        } else if (c2After >= kMeanAfter && c2After <= pMeanAfter) {
            kaonLeftAfter = c2After;
        } else if (kaonToPionRatioAfter >= 0) {
            kaonLeftAfter = (abs(c1After - kMeanAfter) < abs(c2After - kMeanAfter)) ? c1After : c2After;
            cout << "No intersection between kMeanAfter = " << kMeanAfter << " and pMeanAfter = " << pMeanAfter 
                 << ", selecting point closest to kaon peak: " << kaonLeftAfter << endl;
        }

        // Save contamination (after cut) to CSV
        if (kaonToPionRatioAfter < 0 || kaonLeftAfter == -1.0) {
            cout << "Bin [" << pLow << "-" << pHigh << "): Contamination (after cut) calculation failed\n";
            csvFileAfter << pLow << "-" << pHigh << ",N/A,N/A\n";
        } else {
            cout << "Bin [" << pLow << "-" << pHigh << "): After cut - Pion range: [" << pionLeftAfter[i] << ", " << pionRightAfter[i] 
                 << "], Kaon range: [" << kaonLeftAfter << ", " << kaonRightAfter[i] << "], Contamination: " 
                 << kaonToPionRatioAfter << "%\n";
            csvFileAfter << pLow << "-" << pHigh << "," << kaonLeftAfter << "," << kaonToPionRatioAfter << "\n";
        }

        // Clean up histograms and fits
        delete chi2PionsBefore;
        delete chi2KaonsBefore;
        delete betaPionsBefore;
        delete betaKaonsBefore;
        delete chi2PionsAfter;
        delete chi2KaonsAfter; 
        delete betaPionsAfter;
        delete betaKaonsAfter;
        //delete chi2pidall; 
       // delete chi2pidallafter; 
       // delete chi2pidposbefore; 
        //delete chi2pidposafter; 
        delete fitPion;
    }

    // Step 10: Close PDF and CSV files
    
    /* canvaschi2pidall->Print("output/pdf/kaon/chi2pidforall.pdf]"); 
    canvaschi2pidallafter->Print("output/pdf/kaon/chi2pidforallafter.pdf]"); 
    canvaschi2pidposbefore->Print("output/pdf/kaon/chi2pidpos_before.pdf]"); 
    canvaschi2pidposafter->Print("output/pdf/kaon/chi2pidpos_after.pdf]");  */
    canvasChi2Before->Print("output/pdf/kaon/chi2pid_fits_linear_before.pdf]");
    canvasChi2After->Print("output/pdf/kaon/chi2pid_fits_linear_after.pdf]");
    
    canvasBetaBefore->Print("output/pdf/kaon/beta_before_cut_linear.pdf]");
    canvasBetaAfter->Print("output/pdf/kaon/beta_after_cut_linear.pdf]");
    csvFileBefore.close();
    csvFileAfter.close();
    csvFileCuts.close();

    // Step 11: Clean up
    
    delete canvasChi2Before;
    delete canvasChi2After;
    delete canvasBetaBefore;
    delete canvasBetaAfter;
    /*
    delete canvaschi2pidall; 
    delete canvaschi2pidallafter; 
    delete canvaschi2pidposbefore; 
    delete canvaschi2pidposafter; 
    */
    file->Close();
    delete file;

    return 0;
}