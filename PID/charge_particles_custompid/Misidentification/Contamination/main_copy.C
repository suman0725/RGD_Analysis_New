/* #include "include/ContaminationUtils.h"
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

int main() {
    // Load ROOT file and trees
    TFile* file = new TFile("../pkptreeCxC_6.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file ../pkptreeCxC_4.root\n";
        return 1;
    }
    TTree* treePions = (TTree*)file->Get("EB_pid_pions");
    TTree* treeKaons = (TTree*)file->Get("EB_pid_kaons");
    if (!treePions || !treeKaons) {
        std::cerr << "Error: Cannot load trees\n";
        file->Close();
        delete file;
        return 1;
    }

    // Setup output directories
    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/pdf", kTRUE);
    gSystem->mkdir("output/csv", kTRUE);

    // Define bins (1 to 4 GeV, 10 uniform bins)
    std::vector<double> pBins;
    double pMin = 1.0, pMax = 7.0;
    int nBins = 20;
    double binWidth = (pMax - pMin) / nBins;
    for (int i = 0; i <= nBins; i++) pBins.push_back(pMin + i * binWidth);

    // Define chi2pid histogram ranges, fit ranges, and fit parameters per bin
    std::vector<std::pair<double, double>> chi2HistRanges(nBins);
    std::vector<std::pair<double, double>> chi2FitRanges(nBins);
    std::vector<std::tuple<double, double, double>> chi2FitParams(nBins);
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

    // Define beta histogram ranges and fit ranges per bin
    std::vector<std::pair<double, double>> betaHistRanges(nBins);
    std::vector<std::pair<double, double>> betaFitRanges(nBins);
    for (int i = 0; i < nBins; i++) {
        betaHistRanges[i] = {0.95, 1.03};
        betaFitRanges[i] = {0.95, 1.02};
    }

    // Initialize canvases and open PDFs
    TCanvas* canvasChi2 = new TCanvas("canvasChi2", "Chi2pid Fits", 1200, 600);
    TCanvas* canvasBeta = new TCanvas("canvasBeta", "Beta Plots", 800, 600);
    TCanvas* canvasBetaPionCut = new TCanvas("canvasBetaPionCut", "Beta Pion Cut Only", 800, 600);
    TCanvas* canvasBetaPionKaonCut = new TCanvas("canvasBetaPionKaonCut", "Beta Pion Kaon Cut", 800, 600);
    canvasChi2->Print("output/pdf/chi2pid_fits_linear.pdf[");
    canvasBeta->Print("output/pdf/beta_before_cut_linear.pdf[");
    canvasBetaPionCut->Print("output/pdf/beta_pion_cut_only_linear.pdf[");
    canvasBetaPionKaonCut->Print("output/pdf/beta_after_pion_kaon_cut_linear.pdf[");
    std::ofstream csvFile("output/csv/contamination_1to4GeV_linear.csv");
    std::ofstream csvFileBefore("output/csv/contamination_before_chi2pid_1to4GeV_linear.csv");
    std::ofstream csvFilePionCut("output/csv/contamination_pion_cut_only_1to4GeV_linear.csv");
    csvFile << "Momentum Bin (GeV/c),c,Contamination (%)\n";
    csvFileBefore << "Momentum Bin (GeV/c),c,Contamination (%)\n";
    csvFilePionCut << "Momentum Bin (GeV/c),c,Contamination (%)\n";

    // Process each bin
    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];

        // Create chi2pid histograms
        TH1F* chi2Pions = new TH1F(TString::Format("chi2_pions_%d", i),
                                   TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                   100, chi2HistRanges[i].first, chi2HistRanges[i].second);
        TH1F* chi2Kaons = new TH1F(TString::Format("chi2_kaons_%d", i),
                                   TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                   100, chi2HistRanges[i].first, chi2HistRanges[i].second);

        // Fill chi2pid histograms with momentum cut
        float p, chi2pid, beta;
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("chi2pid", &chi2pid);
        for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh) chi2Pions->Fill(chi2pid);
        }
        treeKaons->SetBranchAddress("p", &p);
        treeKaons->SetBranchAddress("chi2pid", &chi2pid);
        for (Long64_t j = 0; j < treeKaons->GetEntries(); j++) {
            treeKaons->GetEntry(j);
            if (p >= pLow && p < pHigh) chi2Kaons->Fill(chi2pid);
        }

        // Check for sufficient entries
        if (chi2Pions->GetEntries() < 10 || chi2Kaons->GetEntries() < 10) {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Too few chi2pid entries for "
                      << (chi2Pions->GetEntries() < 10 ? "pions" : "kaons") << ": "
                      << (chi2Pions->GetEntries() < 10 ? chi2Pions->GetEntries() : chi2Kaons->GetEntries()) << "\n";
            csvFile << pLow << "-" << pHigh << ",N/A,N/A\n";
            csvFileBefore << pLow << "-" << pHigh << ",N/A,N/A\n";
            csvFilePionCut << pLow << "-" << pHigh << ",N/A,N/A\n";
            delete chi2Pions;
            delete chi2Kaons;
            continue;
        }

        // Set colors
        chi2Pions->SetLineColor(kBlue);
        chi2Kaons->SetLineColor(kGreen);

        // Fit chi2pid histograms
        double pionChi2Min, pionChi2Max, kaonChi2Min, kaonChi2Max;
        auto [pionAmp, pionMean, pionSigma] = chi2FitParams[i];
        auto [kaonAmp, kaonMean, kaonSigma] = chi2FitParams[i];
        TF1* fitPion = new TF1("gausFit", "gaus", chi2FitRanges[i].first, chi2FitRanges[i].second);
        fitPion->SetParameters(pionAmp, pionMean, pionSigma);
        chi2Pions->Fit(fitPion, "R", "", chi2FitRanges[i].first, chi2FitRanges[i].second);
        pionChi2Min = fitPion->GetParameter(1) - 3 * fitPion->GetParameter(2);
        pionChi2Max = fitPion->GetParameter(1) + 3 * fitPion->GetParameter(2);
        TF1* fitKaon = new TF1("gausFit", "gaus", chi2FitRanges[i].first, chi2FitRanges[i].second);
        fitKaon->SetParameters(kaonAmp, kaonMean, kaonSigma);
        chi2Kaons->Fit(fitKaon, "R", "", chi2FitRanges[i].first, chi2FitRanges[i].second);
        kaonChi2Min = fitKaon->GetParameter(1) - 3 * fitKaon->GetParameter(2);
        kaonChi2Max = fitKaon->GetParameter(1) + 3 * fitKaon->GetParameter(2);
        std::cout << "Bin [" << pLow << "-" << pHigh << "): Pions chi2 range: [" 
                  << pionChi2Min << ", " << pionChi2Max << "], Kaons chi2 range: [" 
                  << kaonChi2Min << ", " << kaonChi2Max << "]\n";

        // Plot chi2pid histograms
        canvasChi2->Clear();
        canvasChi2->Divide(2, 1);
        canvasChi2->cd(1);
        chi2Pions->Draw("HIST");
        if (chi2Pions->GetFunction("gausFit")) {
            chi2Pions->GetFunction("gausFit")->SetLineColor(kBlue);
            chi2Pions->GetFunction("gausFit")->Draw("SAME");
        }
        TLegend* legPions = new TLegend(0.1, 0.75, 0.3, 0.9);
        legPions->AddEntry(chi2Pions, "Pions", "l");
        legPions->SetBorderSize(0);           // Remove border
        legPions->SetFillStyle(0);            // Transparent background
        legPions->SetTextSize(0.025);
        legPions->Draw();
        canvasChi2->cd(2);
        chi2Kaons->Draw("HIST");
        if (chi2Kaons->GetFunction("gausFit")) {
            chi2Kaons->GetFunction("gausFit")->SetLineColor(kGreen);
            chi2Kaons->GetFunction("gausFit")->Draw("SAME");
        }
        TLegend* legKaons = new TLegend(0.1, 0.75, 0.3, 0.9);
        legKaons->AddEntry(chi2Kaons, "Kaons", "l");
        legKaons->SetBorderSize(0);           // Remove border
        legKaons->SetFillStyle(0);            // Transparent background
        legKaons->SetTextSize(0.025);
        legKaons->Draw();
        canvasChi2->Update();
        canvasChi2->Print("output/pdf/chi2pid_fits_linear.pdf");
        delete legPions;
        delete legKaons;
        delete chi2Pions;
        delete chi2Kaons;

        // Define beta cuts for "before chi2pid cuts" section
        double pionLeftBefore, pionRightBefore, kaonRightBefore;
        if (i == 0) { // Bin [1.0-1.3)
            pionLeftBefore = 0.9880;
            pionRightBefore = 1.0120;
            kaonRightBefore = 0.0;
        } else if (i == 1) { // Bin [1.3-1.6)
            pionLeftBefore = 0.9755;
            pionRightBefore = 1.0135;
            kaonRightBefore = 0.9721;
        } else if (i == 2) { // Bin [1.6-1.9)
            pionLeftBefore = 0.98;
            pionRightBefore = 1.0140;
            kaonRightBefore = 0.986;
        } else if (i == 3) { // Bin [1.9-2.2)
            pionLeftBefore = 0.9815;
            pionRightBefore = 1.0145;
            kaonRightBefore = 0.992;
        } else if (i == 4) { // Bin [2.2-2.5)
            pionLeftBefore = 0.9880;
            pionRightBefore = 1.0141;
            kaonRightBefore = 0.9971;
        } else if (i == 5) { // Bin [2.5-2.8)
            pionLeftBefore = 0.9848;
            pionRightBefore = 1.013;
            kaonRightBefore = 1.0;
        } else if (i == 6) { // Bin [2.8-3.1)
            pionLeftBefore = 0.986;
            pionRightBefore = 1.0128;
            kaonRightBefore = 1.003;
        } else if (i == 7) { // Bin [3.1-3.4)
            pionLeftBefore = 0.9868;
            pionRightBefore = 1.0121;
            kaonRightBefore = 1.0048;
        } else if (i == 8) { // Bin [3.4-3.7)
            pionLeftBefore = 0.9879;
            pionRightBefore = 1.0120;
            kaonRightBefore = 1.0049;
        } else if (i == 9) { // Bin [3.7-4.0)
            pionLeftBefore = 0.9880;
            pionRightBefore = 1.0120;
            kaonRightBefore = 1.012;
        } else {
            std::cerr << "Warning: Using default beta cuts for bin [" << pLow << "-" << pHigh << ") in before chi2pid section\n";
            pionLeftBefore = 0.9880;
            pionRightBefore = 1.0120;
            kaonRightBefore = 1.0;
        }
        if (pionLeftBefore >= pionRightBefore) {
            std::cerr << "Error: Invalid pion cuts for bin [" << pLow << "-" << pHigh << ") in before chi2pid section\n";
            csvFileBefore << pLow << "-" << pHigh << ",N/A,N/A\n";
            continue;
        }

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
            if (p >= pLow && p < pHigh && beta >= betaHistRanges[i].first && beta <= betaHistRanges[i].second) {
                betaPionsBefore->Fill(beta);
            }
        }
        treeKaons->SetBranchAddress("p", &p);
        treeKaons->SetBranchAddress("beta", &beta);
        for (Long64_t j = 0; j < treeKaons->GetEntries(); j++) {
            treeKaons->GetEntry(j);
            if (p >= pLow && p < pHigh && beta >= betaHistRanges[i].first && beta <= betaHistRanges[i].second) {
                betaKaonsBefore->Fill(beta);
            }
        }

        // Check for sufficient entries
        if (betaPionsBefore->GetEntries() < 10 || betaKaonsBefore->GetEntries() < 10) {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Too few beta entries for "
                      << (betaPionsBefore->GetEntries() < 10 ? "pions" : "kaons") << ": "
                      << (betaPionsBefore->GetEntries() < 10 ? betaPionsBefore->GetEntries() : betaKaonsBefore->GetEntries()) << "\n";
            csvFile << pLow << "-" << pHigh << ",N/A,N/A\n";
            csvFileBefore << pLow << "-" << pHigh << ",N/A,N/A\n";
            csvFilePionCut << pLow << "-" << pHigh << ",N/A,N/A\n";
            delete betaPionsBefore;
            delete betaKaonsBefore;
            continue;
        }

        // Fit beta histograms (before cut) for contamination
        double pMeanBefore = 0, pSigmaBefore = 0, pConstBefore = 0, kMeanBefore = 0, kSigmaBefore = 0, kConstBefore = 0;
        if (betaPionsBefore->GetEntries() >= 10) {
            TF1* gausPionBefore = new TF1(TString::Format("gaus_%s", betaPionsBefore->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            gausPionBefore->SetParameters(betaPionsBefore->GetMaximum(), 1.0, 0.005);
            betaPionsBefore->Fit(gausPionBefore, "R", "", betaFitRanges[i].first, betaFitRanges[i].second);
            pMeanBefore = gausPionBefore->GetParameter(1);
            pSigmaBefore = gausPionBefore->GetParameter(2);
            pConstBefore = gausPionBefore->GetParameter(0);
            delete gausPionBefore;
        }
        if (betaKaonsBefore->GetEntries() >= 10) {
            TF1* gausKaonBefore = new TF1(TString::Format("gaus_%s", betaKaonsBefore->GetName()), "gaus", 0.95, 1.03);
            if (i == 9) { // Last bin [3.7-4.0) GeV/c
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.992, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.986, 0.9945);
            } else if (i == 8) { // Bin [3.4-3.7) GeV/c
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.987, 0.994);
            } else if (i == 5) { // Bin [2.5-2.8) GeV/c
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.974, 0.992);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.974, 0.992);
            } else {
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.99, 0.005);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.95, 1.0);
            }
            kMeanBefore = gausKaonBefore->GetParameter(1);
            kSigmaBefore = gausKaonBefore->GetParameter(2);
            kConstBefore = gausKaonBefore->GetParameter(0);
            delete gausKaonBefore;
        }

        // Plot beta histograms (before cut)
        canvasBeta->Clear();
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
            }
        }
        if (betaKaonsBefore->GetEntries() >= 10) {
            betaKaonsBefore->SetMarkerStyle(20);
            betaKaonsBefore->SetMarkerSize(1);
            betaKaonsBefore->SetMarkerColor(kGreen);
            betaKaonsBefore->Draw(betaPionsBefore->GetEntries() >= 10 ? "P SAME" : "P");
            yMax = std::max(yMax, betaKaonsBefore->GetMaximum());
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
        canvasBeta->SetGrid();
        canvasBeta->Update();
        canvasBeta->Print("output/pdf/beta_before_cut_linear.pdf");
        delete leg;

        // Calculate contamination (before chi2pid cuts)
        double c1Before, c2Before;
        std::map<std::string, double> pFitBefore = {{"gaus_mean", pMeanBefore}, {"gaus_sigma", pSigmaBefore}, {"gaus_constant", pConstBefore}};
        std::map<std::string, double> kFitBefore = {{"gaus_mean", kMeanBefore}, {"gaus_sigma", kSigmaBefore}, {"gaus_constant", kConstBefore}};
        double contaminationBefore = calculateOverlapContamination(pFitBefore, kFitBefore, pionLeftBefore, pionRightBefore, kaonRightBefore, c1Before, c2Before);
        double kaonLeftBefore = (c1Before >= kMeanBefore && c1Before <= pMeanBefore) ? c1Before : (c2Before >= kMeanBefore && c2Before <= pMeanBefore) ? c2Before : -1.0;

        // Save contamination (before chi2pid cuts)
        if (contaminationBefore < 0 || kaonLeftBefore == -1.0) {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Contamination (before chi2pid) calculation failed\n";
            csvFileBefore << pLow << "-" << pHigh << ",N/A,N/A\n";
        } else {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Before chi2pid - Pion range: [" << pionLeftBefore << ", " << pionRightBefore 
                      << "], Kaon range: [" << kaonLeftBefore << ", " << kaonRightBefore << "], Contamination: " 
                      << contaminationBefore << "%\n";
            csvFileBefore << pLow << "-" << pHigh << "," << kaonLeftBefore << "," << contaminationBefore << "\n";
        }
        delete betaPionsBefore;
        delete betaKaonsBefore;

        // Define beta cuts for "pion chi2pid cut only" section
        double pionLeftPionCut, pionRightPionCut, kaonRightPionCut;
        if (i == 0) { // Bin [1.0-1.3)
            pionLeftPionCut = 0.9880;  // Same as before, editable
            pionRightPionCut = 1.0120;
            kaonRightPionCut = 0.0;
        } else if (i == 1) { // Bin [1.3-1.6)
            pionLeftPionCut = 0.9755;
            pionRightPionCut = 1.0135;
            kaonRightPionCut = 0.9721;
        } else if (i == 2) { // Bin [1.6-1.9)
            pionLeftPionCut = 0.98;
            pionRightPionCut = 1.0140;
            kaonRightPionCut = 0.986;
        } else if (i == 3) { // Bin [1.9-2.2)
            pionLeftPionCut = 0.9815;
            pionRightPionCut = 1.0145;
            kaonRightPionCut = 0.992;
        } else if (i == 4) { // Bin [2.2-2.5)
            pionLeftPionCut = 0.9880;
            pionRightPionCut = 1.0141;
            kaonRightPionCut = 0.9971;
        } else if (i == 5) { // Bin [2.5-2.8)
            pionLeftPionCut = 0.9848;
            pionRightPionCut = 1.013;
            kaonRightPionCut = 1.0;
        } else if (i == 6) { // Bin [2.8-3.1)
            pionLeftPionCut = 0.986;
            pionRightPionCut = 1.0128;
            kaonRightPionCut = 1.003;
        } else if (i == 7) { // Bin [3.1-3.4)
            pionLeftPionCut = 0.9868;
            pionRightPionCut = 1.0121;
            kaonRightPionCut = 1.0048;
        } else if (i == 8) { // Bin [3.4-3.7)
            pionLeftPionCut = 0.9879;
            pionRightPionCut = 1.0120;
            kaonRightPionCut = 1.0049;
        } else if (i == 9) { // Bin [3.7-4.0)
            pionLeftPionCut = 0.9880;
            pionRightPionCut = 1.0120;
            kaonRightPionCut = 1.012;
        } else {
            std::cerr << "Warning: Using default beta cuts for bin [" << pLow << "-" << pHigh << ") in pion chi2pid cut only section\n";
            pionLeftPionCut = 0.9880;
            pionRightPionCut = 1.0120;
            kaonRightPionCut = 1.0;
        }
        if (pionLeftPionCut >= pionRightPionCut) {
            std::cerr << "Error: Invalid pion cuts for bin [" << pLow << "-" << pHigh << ") in pion chi2pid cut only section\n";
            csvFilePionCut << pLow << "-" << pHigh << ",N/A,N/A\n";
            continue;
        }

        // Create beta histograms (pion chi2pid cut only)
        TH1F* betaPionsAfterPionCut = new TH1F(TString::Format("beta_pions_after_pioncut_%d", i),
                                               TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                               70, betaHistRanges[i].first, betaHistRanges[i].second);
        TH1F* betaKaonsAfterPionCut = new TH1F(TString::Format("beta_kaons_after_pioncut_%d", i),
                                               TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                               70, betaHistRanges[i].first, betaHistRanges[i].second);

        betaPionsAfterPionCut->Sumw2();
        betaKaonsAfterPionCut->Sumw2();

        // Fill beta histograms (pion chi2pid cut only)
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("chi2pid", &chi2pid);
        treePions->SetBranchAddress("beta", &beta);
        for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh && chi2pid >= pionChi2Min && chi2pid <= pionChi2Max &&
                beta >= betaHistRanges[i].first && beta <= betaHistRanges[i].second) {
                betaPionsAfterPionCut->Fill(beta);
            }
        }
        treeKaons->SetBranchAddress("p", &p);
        treeKaons->SetBranchAddress("beta", &beta);
        for (Long64_t j = 0; j < treeKaons->GetEntries(); j++) {
            treeKaons->GetEntry(j);
            if (p >= pLow && p < pHigh && beta >= betaHistRanges[i].first && beta <= betaHistRanges[i].second) {
                betaKaonsAfterPionCut->Fill(beta);
            }
        }

        // Check for sufficient entries
        if (betaPionsAfterPionCut->GetEntries() < 10 && betaKaonsAfterPionCut->GetEntries() < 10) {
            csvFilePionCut << pLow << "-" << pHigh << ",N/A,N/A\n";
            delete betaPionsAfterPionCut;
            delete betaKaonsAfterPionCut;
            continue;
        }

        // Fit beta histograms (pion chi2pid cut only)
        double pMeanPionCut = 0, pSigmaPionCut = 0, pConstPionCut = 0, kMeanPionCut = 0, kSigmaPionCut = 0, kConstPionCut = 0;
        if (betaPionsAfterPionCut->GetEntries() >= 10) {
            TF1* gausPionPionCut = new TF1(TString::Format("gaus_%s", betaPionsAfterPionCut->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            gausPionPionCut->SetParameters(betaPionsAfterPionCut->GetMaximum(), 1.0, 0.005);
            betaPionsAfterPionCut->Fit(gausPionPionCut, "R", "", betaFitRanges[i].first, betaFitRanges[i].second);
            pMeanPionCut = gausPionPionCut->GetParameter(1);
            pSigmaPionCut = gausPionPionCut->GetParameter(2);
            pConstPionCut = gausPionPionCut->GetParameter(0);
            delete gausPionPionCut;
        }
        if (betaKaonsAfterPionCut->GetEntries() >= 10) {
            TF1* gausKaonPionCut = new TF1(TString::Format("gaus_%s", betaKaonsAfterPionCut->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            if (i == 9) { // Last bin [3.7-4.0) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.1);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.986, 0.9945);
            } else if (i == 8) { // Bin [3.4-3.7) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.1);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.988, 0.9945);
            } else if (i == 10) { // Bin [3.4-3.7) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.1);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.987, 0.995);
            }
            else if (i == 11) { // Bin [3.4-3.7) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.996, 0.012);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.992, 0.996);
            }
            else if (i == 12) { // Bin [3.4-3.7) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.012);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.991, 0.996);
            }
            else if (i == 13) { // Bin [3.4-3.7) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.012);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.991, 0.995);
            }
            else if (i == 14) { // Bin [3.4-3.7) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.1);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.991, 0.996);
            }
            else if (i == 15) { // Bin [3.4-3.7) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.1);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.991, 0.996);
            }
            else if (i == 16) { // Bin [3.4-3.7) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.1);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.991, 0.997);
            }
            else if (i == 17) { // Bin [3.4-3.7) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.1);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.991, 0.997);
            }
            else if (i == 18) { // Bin [3.4-3.7) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.1);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.982, 0.997);
            }
            
            else {
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.99, 0.005);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R");
            } 
            
            kMeanPionCut = gausKaonPionCut->GetParameter(1);
            kSigmaPionCut = gausKaonPionCut->GetParameter(2);
            kConstPionCut = gausKaonPionCut->GetParameter(0);
            delete gausKaonPionCut;
        }

        // Calculate contamination (pion chi2pid cut only)
        double c1PionCut, c2PionCut;
        std::map<std::string, double> pFitPionCut = {{"gaus_mean", pMeanPionCut}, {"gaus_sigma", pSigmaPionCut}, {"gaus_constant", pConstPionCut}};
        std::map<std::string, double> kFitPionCut = {{"gaus_mean", kMeanPionCut}, {"gaus_sigma", kSigmaPionCut}, {"gaus_constant", kConstPionCut}};
        double contaminationPionCut = calculateOverlapContamination(pFitPionCut, kFitPionCut, pionLeftPionCut, pionRightPionCut, kaonRightPionCut, c1PionCut, c2PionCut);
        double kaonLeftPionCut = -1.0;
        if (c1PionCut >= kMeanPionCut && c1PionCut <= pMeanPionCut) {
            kaonLeftPionCut = c1PionCut;
        } else if (c2PionCut >= kMeanPionCut && c2PionCut <= pMeanPionCut) {
            kaonLeftPionCut = c2PionCut;
        } else if (contaminationPionCut >= 0) {
            // If no intersection lies between kMeanPionCut and pMeanPionCut, choose the point closest to kaon mean
            kaonLeftPionCut = (std::abs(c1PionCut - kMeanPionCut) < std::abs(c2PionCut - kMeanPionCut)) ? c1PionCut : c2PionCut;
            std::cout << "No intersection between kMeanPionCut = " << kMeanPionCut << " and pMeanPionCut = " << pMeanPionCut 
                      << ", selecting point closest to kaon peak: " << kaonLeftPionCut << std::endl;
        }

        // Plot beta histograms (pion chi2pid cut only)
        canvasBetaPionCut->Clear();
        logScale = false;
        yMin = 0.0;
        yMax = -1.0;
        if (betaPionsAfterPionCut->GetEntries() >= 10) {
            betaPionsAfterPionCut->SetMarkerStyle(20);
            betaPionsAfterPionCut->SetMarkerSize(1);
            betaPionsAfterPionCut->SetMarkerColor(kBlue);
            betaPionsAfterPionCut->Draw("P");
            yMax = betaPionsAfterPionCut->GetMaximum();
            if (auto* fit = betaPionsAfterPionCut->GetFunction(TString::Format("gaus_%s", betaPionsAfterPionCut->GetName()))) {
                fit->SetLineColor(kBlue);
                fit->SetLineStyle(kDashed);
                fit->Draw("SAME");
            }
        }
        if (betaKaonsAfterPionCut->GetEntries() >= 10) {
            betaKaonsAfterPionCut->SetMarkerStyle(20);
            betaKaonsAfterPionCut->SetMarkerSize(1);
            betaKaonsAfterPionCut->SetMarkerColor(kGreen);
            betaKaonsAfterPionCut->Draw(betaPionsAfterPionCut->GetEntries() >= 10 ? "P SAME" : "P");
            yMax = std::max(yMax, betaKaonsAfterPionCut->GetMaximum());
            if (auto* fit = betaKaonsAfterPionCut->GetFunction(TString::Format("gaus_%s", betaKaonsAfterPionCut->GetName()))) {
                fit->SetLineColor(kGreen);
                fit->SetLineStyle(kDashed);
                fit->Draw("SAME");
                if ( i >= 11){fit->SetRange(0.98, 1.02);} 
                else {fit->SetRange(0.95, 1.03);}
            }
        }
        yMax *= 1.2;
        if (betaPionsAfterPionCut->GetEntries() >= 10 || betaKaonsAfterPionCut->GetEntries() >= 10) {
            (betaPionsAfterPionCut->GetEntries() >= 10 ? betaPionsAfterPionCut : betaKaonsAfterPionCut)->GetYaxis()->SetRangeUser(yMin, yMax);
        }
        if (kaonLeftPionCut >= 0) {
            TLine* line = new TLine(kaonLeftPionCut, yMin, kaonLeftPionCut, yMax);
            line->SetLineColor(kRed);
            line->SetLineStyle(kSolid);
            line->SetLineWidth(2);
            line->Draw("SAME");
            delete line;
        }
        leg = new TLegend(0.15, 0.75, 0.35, 0.88);
        if (betaPionsAfterPionCut->GetEntries() >= 10) leg->AddEntry(betaPionsAfterPionCut, "Pions", "p");
        if (betaKaonsAfterPionCut->GetEntries() >= 10) leg->AddEntry(betaKaonsAfterPionCut, "Kaons", "p");
        leg->SetTextSize(0.03);
        leg->Draw();
        canvasBetaPionCut->SetGrid();
        canvasBetaPionCut->Update();
        canvasBetaPionCut->Print("output/pdf/beta_pion_cut_only_linear.pdf");
        delete leg;

        // Save contamination (pion chi2pid cut only)
        if (contaminationPionCut < 0 || kaonLeftPionCut == -1.0) {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Contamination (pion cut only) calculation failed\n";
            csvFilePionCut << pLow << "-" << pHigh << ",N/A,N/A\n";
        } else {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Pion cut only - Pion range: [" << pionLeftPionCut << ", " << pionRightPionCut 
                      << "], Kaon range: [" << kaonLeftPionCut << ", " << kaonRightPionCut << "], Contamination: " 
                      << contaminationPionCut << "%\n";
            csvFilePionCut << pLow << "-" << pHigh << "," << kaonLeftPionCut << "," << contaminationPionCut << "\n";
        }
        delete betaPionsAfterPionCut;
        delete betaKaonsAfterPionCut;

        // Define beta cuts for "pion and kaon chi2pid cuts" section
        double pionLeftBoth, pionRightBoth, kaonRightBoth;
        if (i == 0) { // Bin [1.0-1.3)
            pionLeftBoth = 0.9880;  // Same as before, editable
            pionRightBoth = 1.0120;
            kaonRightBoth = 0.0;
        } else if (i == 1) { // Bin [1.3-1.6)
            pionLeftBoth = 0.9755;
            pionRightBoth = 1.0135;
            kaonRightBoth = 0.9721;
        } else if (i == 2) { // Bin [1.6-1.9)
            pionLeftBoth = 0.98;
            pionRightBoth = 1.0140;
            kaonRightBoth = 0.986;
        } else if (i == 3) { // Bin [1.9-2.2)
            pionLeftBoth = 0.9815;
            pionRightBoth = 1.0145;
            kaonRightBoth = 0.992;
        } else if (i == 4) { // Bin [2.2-2.5)
            pionLeftBoth = 0.9880;
            pionRightBoth = 1.0141;
            kaonRightBoth = 0.9971;
        } else if (i == 5) { // Bin [2.5-2.8)
            pionLeftBoth = 0.9848;
            pionRightBoth = 1.013;
            kaonRightBoth = 1.0;
        } else if (i == 6) { // Bin [2.8-3.1)
            pionLeftBoth = 0.986;
            pionRightBoth = 1.0128;
            kaonRightBoth = 1.003;
        } else if (i == 7) { // Bin [3.1-3.4)
            pionLeftBoth = 0.9868;
            pionRightBoth = 1.0121;
            kaonRightBoth = 1.0048;
        } else if (i == 8) { // Bin [3.4-3.7)
            pionLeftBoth = 0.9879;
            pionRightBoth = 1.0120;
            kaonRightBoth = 1.0049;
        } else if (i == 9) { // Bin [3.7-4.0)
            pionLeftBoth = 0.984;
            pionRightBoth = 1.0120;
            kaonRightBoth = 1.012;
        } else {
            std::cerr << "Warning: Using default beta cuts for bin [" << pLow << "-" << pHigh << ") in pion and kaon chi2pid cuts section\n";
            pionLeftBoth = 0.9880;
            pionRightBoth = 1.0120;
            kaonRightBoth = 1.0;
        }
        if (pionLeftBoth >= pionRightBoth) {
            std::cerr << "Error: Invalid pion cuts for bin [" << pLow << "-" << pHigh << ") in pion and kaon chi2pid cuts section\n";
            csvFile << pLow << "-" << pHigh << ",N/A,N/A\n";
            continue;
        }

        // Create beta histograms (after cut, both pions and kaons)
        TH1F* betaPionsAfter = new TH1F(TString::Format("beta_pions_after_%d", i),
                                        TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                        70, betaHistRanges[i].first, betaHistRanges[i].second);
        TH1F* betaKaonsAfter = new TH1F(TString::Format("beta_kaons_after_%d", i),
                                        TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                        70, betaHistRanges[i].first, betaHistRanges[i].second);

        betaPionsAfter->Sumw2();
        betaKaonsAfter->Sumw2();

        // Fill beta histograms (after cut) with momentum and chi2pid cuts
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("chi2pid", &chi2pid);
        treePions->SetBranchAddress("beta", &beta);
        for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh && chi2pid >= pionChi2Min && chi2pid <= pionChi2Max &&
                beta >= betaHistRanges[i].first && beta <= betaHistRanges[i].second) {
                betaPionsAfter->Fill(beta);
            }
        }
        treeKaons->SetBranchAddress("p", &p);
        treeKaons->SetBranchAddress("chi2pid", &chi2pid);
        treeKaons->SetBranchAddress("beta", &beta);
        for (Long64_t j = 0; j < treeKaons->GetEntries(); j++) {
            treeKaons->GetEntry(j);
            if (p >= pLow && p < pHigh && chi2pid >= kaonChi2Min && chi2pid <= kaonChi2Max &&
                beta >= betaHistRanges[i].first && beta <= betaHistRanges[i].second) {
                betaKaonsAfter->Fill(beta);
            }
        }

        // Check for sufficient entries
        if (betaPionsAfter->GetEntries() < 10 && betaKaonsAfter->GetEntries() < 10) {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Too few beta entries for both pion and kaon cuts\n";
            csvFile << pLow << "-" << pHigh << ",N/A,N/A\n";
            delete betaPionsAfter;
            delete betaKaonsAfter;
            continue;
        }

        // Fit beta histograms (after cut, both pions and kaons)
        double pMean = 0, pSigma = 0, pConst = 0, kMean = 0, kSigma = 0, kConst = 0;
        if (betaPionsAfter->GetEntries() >= 10) {
            TF1* gausPion = new TF1(TString::Format("gaus_%s", betaPionsAfter->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            gausPion->SetParameters(betaPionsAfter->GetMaximum(), 1.0, 0.005);
            betaPionsAfter->Fit(gausPion, "R", "", betaFitRanges[i].first, betaFitRanges[i].second);
            pMean = gausPion->GetParameter(1);
            pSigma = gausPion->GetParameter(2);
            pConst = gausPion->GetParameter(0);
            delete gausPion;
        }
        if (betaKaonsAfter->GetEntries() >= 10) {
            TF1* gausKaon = new TF1(TString::Format("gaus_%s", betaKaonsAfter->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            if (i == 9) { // Last bin [3.7-4.0) GeV/c
                gausKaon->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
                betaKaonsAfter->Fit(gausKaon, "R", "", 0.982, 0.996);
            } else if( i == 8 ){
                gausKaon->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
                betaKaonsAfter->Fit(gausKaon, "R", "", 0.982, 0.995);
            }else if( i == 7 ){
                gausKaon->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
                betaKaonsAfter->Fit(gausKaon, "R", "", 0.98, 0.992);
            }
            else if (i == 6) { // Bin [2.5-2.8) GeV/c
                gausKaon->SetParameters(betaKaonsAfter->GetMaximum(), 0.974, 0.992);
                betaKaonsAfter->Fit(gausKaon, "R", "", 0.974, 0.992);
            }
            
            else {
                gausKaon->SetParameters(betaKaonsAfter->GetMaximum(), 0.99, 0.005);
                betaKaonsAfter->Fit(gausKaon, "R");
            }  
            
            

            kMean = gausKaon->GetParameter(1);
            kSigma = gausKaon->GetParameter(2);
            kConst = gausKaon->GetParameter(0);
            delete gausKaon;
        }

        // Calculate contamination (after chi2pid cuts, both pions and kaons)
        double c1, c2;
        std::map<std::string, double> pFit = {{"gaus_mean", pMean}, {"gaus_sigma", pSigma}, {"gaus_constant", pConst}};
        std::map<std::string, double> kFit = {{"gaus_mean", kMean}, {"gaus_sigma", kSigma}, {"gaus_constant", kConst}};
        double contamination = calculateOverlapContamination(pFit, kFit, pionLeftBoth, pionRightBoth, kaonRightBoth, c1, c2);
        double kaonLeft = (c1 >= kMean && c1 <= pMean) ? c1 : (c2 >= kMean && c2 <= pMean) ? c2 : -1.0;

        // Plot beta histograms (pion and kaon chi2pid cut)
        canvasBetaPionKaonCut->Clear();
        logScale = false;
        yMin = 0.0;
        yMax = -1.0;
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
            }
        }
        if (betaKaonsAfter->GetEntries() >= 10) {
            betaKaonsAfter->SetMarkerStyle(20);
            betaKaonsAfter->SetMarkerSize(1);
            betaKaonsAfter->SetMarkerColor(kGreen);
            betaKaonsAfter->Draw(betaPionsAfter->GetEntries() >= 10 ? "P SAME" : "P");
            yMax = std::max(yMax, betaKaonsAfter->GetMaximum());
            if (auto* fit = betaKaonsAfter->GetFunction(TString::Format("gaus_%s", betaKaonsAfter->GetName()))) {
                fit->SetLineColor(kGreen);
                fit->SetLineStyle(kDashed);
                fit->Draw("SAME");
                fit->SetRange(0.95, 1.03);
            }
        }
        yMax *= 1.2;
        if (betaPionsAfter->GetEntries() >= 10 || betaKaonsAfter->GetEntries() >= 10) {
            (betaPionsAfter->GetEntries() >= 10 ? betaPionsAfter : betaKaonsAfter)->GetYaxis()->SetRangeUser(yMin, yMax);
        }
        if (kaonLeft >= 0) {
            TLine* line = new TLine(kaonLeft, yMin, kaonLeft, yMax);
            line->SetLineColor(kRed);
            line->SetLineStyle(kSolid);
            line->SetLineWidth(2);
            line->Draw("SAME");
            delete line;
        }
        leg = new TLegend(0.15, 0.75, 0.35, 0.88);
        if (betaPionsAfter->GetEntries() >= 10) leg->AddEntry(betaPionsAfter, "Pions", "p");
        if (betaKaonsAfter->GetEntries() >= 10) leg->AddEntry(betaKaonsAfter, "Kaons", "p");
        leg->SetTextSize(0.03);
        leg->Draw();
        canvasBetaPionKaonCut->SetGrid();
        canvasBetaPionKaonCut->Update();
        canvasBetaPionKaonCut->Print("output/pdf/beta_after_pion_kaon_cut_linear.pdf");
        delete leg;

        // Save contamination (after chi2pid cuts, both pions and kaons)
        if (contamination < 0 || kaonLeft == -1.0) {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Contamination (both cuts) calculation failed\n";
            csvFile << pLow << "-" << pHigh << ",N/A,N/A\n";
        } else {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Both cuts - Pion range: [" << pionLeftBoth << ", " << pionRightBoth 
                      << "], Kaon range: [" << kaonLeft << ", " << kaonRightBoth << "], Contamination: " 
                      << contamination << "%\n";
            csvFile << pLow << "-" << pHigh << "," << kaonLeft << "," << contamination << "\n";
        }
        delete betaPionsAfter;
        delete betaKaonsAfter;
    }

    // Close PDFs and CSVs
    canvasChi2->Print("output/pdf/chi2pid_fits_linear.pdf]");
    canvasBeta->Print("output/pdf/beta_before_cut_linear.pdf]");
    canvasBetaPionCut->Print("output/pdf/beta_pion_cut_only_linear.pdf]");
    canvasBetaPionKaonCut->Print("output/pdf/beta_after_pion_kaon_cut_linear.pdf]");
    csvFile.close();
    csvFileBefore.close();
    csvFilePionCut.close();

    // Clean up
    delete canvasChi2;
    delete canvasBeta;
    delete canvasBetaPionCut;
    delete canvasBetaPionKaonCut;
    file->Close();
    delete file;

    return 0;
} */

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

int main() {
    // Load ROOT file and trees
    TFile* file = new TFile("../pkptreeCxC_5.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file ../pkptreeCxC_6.root\n";
        return 1;
    }
    TTree* treePions = (TTree*)file->Get("EB_pid_pions");
    TTree* treeKaons = (TTree*)file->Get("EB_pid_kaons");
    if (!treePions || !treeKaons) {
        std::cerr << "Error: Cannot load trees\n";
        file->Close();
        delete file;
        return 1;
    }

    // Setup output directories
    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/pdf", kTRUE);
    gSystem->mkdir("output/csv", kTRUE);

    // Define bins (1 to 7 GeV, 20 uniform bins)
    std::vector<double> pBins;
    double pMin = 1.0, pMax = 7.0;
    int nBins = 20;
    double binWidth = (pMax - pMin) / nBins;
    for (int i = 0; i <= nBins; i++) pBins.push_back(pMin + i * binWidth);

    // Define chi2pid histogram ranges, fit ranges, and fit parameters per bin
    std::vector<std::pair<double, double>> chi2HistRanges(nBins);
    std::vector<std::pair<double, double>> chi2FitRanges(nBins);
    std::vector<std::tuple<double, double, double>> chi2FitParams(nBins);
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

    // Define beta histogram ranges and fit ranges per bin
    std::vector<std::pair<double, double>> betaHistRanges(nBins);
    std::vector<std::pair<double, double>> betaFitRanges(nBins);
    for (int i = 0; i < nBins; i++) {
        betaHistRanges[i] = {0.95, 1.03};
        betaFitRanges[i] = {0.95, 1.02};
    }

    // Initialize canvases and open PDFs
    TCanvas* canvasChi2 = new TCanvas("canvasChi2", "Chi2pid Fits", 1200, 600);
    TCanvas* canvasBeta = new TCanvas("canvasBeta", "Beta Plots", 800, 600);
    TCanvas* canvasBetaPionCut = new TCanvas("canvasBetaPionCut", "Beta Pion Cut Only", 800, 600);
    TCanvas* canvasBetaSigmaBin3 = new TCanvas("canvasBetaSigmaBin3", "Beta Sigma from Bin 3", 800, 600);
    canvasChi2->Print("output/pdf/chi2pid_fits_linear.pdf[");
    canvasBeta->Print("output/pdf/beta_before_cut_linear.pdf[");
    canvasBetaPionCut->Print("output/pdf/beta_pion_cut_only_linear.pdf[");
    canvasBetaSigmaBin3->Print("output/pdf/beta_sigma_from_bin3_linear.pdf[");
    std::ofstream csvFileBefore("output/csv/contamination_before_chi2pid_1to4GeV_linear.csv");
    std::ofstream csvFilePionCut("output/csv/contamination_pion_cut_only_1to4GeV_linear.csv");
    std::ofstream csvFileSigmaBin3("output/csv/contamination_sigma_from_bin3_1to4GeV_linear.csv");
    csvFileBefore << "Momentum Bin (GeV/c),c,Contamination (%)\n";
    csvFilePionCut << "Momentum Bin (GeV/c),c,Contamination (%)\n";
    csvFileSigmaBin3 << "Momentum Bin (GeV/c),c,Contamination (%)\n";

    // Store sigma from bin 3 ([1.6, 1.9) GeV/c)
    double pionSigmaBin3 = 0.0;

    // Process each bin
    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];

        // Create chi2pid histograms
        TH1F* chi2Pions = new TH1F(TString::Format("chi2_pions_%d", i),
                                   TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                   100, chi2HistRanges[i].first, chi2HistRanges[i].second);
        TH1F* chi2Kaons = new TH1F(TString::Format("chi2_kaons_%d", i),
                                   TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                   100, chi2HistRanges[i].first, chi2HistRanges[i].second);

        // Fill chi2pid histograms with momentum cut
        float p, chi2pid, beta;
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("chi2pid", &chi2pid);
        for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh) chi2Pions->Fill(chi2pid);
        }
        treeKaons->SetBranchAddress("p", &p);
        treeKaons->SetBranchAddress("chi2pid", &chi2pid);
        for (Long64_t j = 0; j < treeKaons->GetEntries(); j++) {
            treeKaons->GetEntry(j);
            if (p >= pLow && p < pHigh) chi2Kaons->Fill(chi2pid);
        }

        // Check for sufficient entries
        if (chi2Pions->GetEntries() < 10 || chi2Kaons->GetEntries() < 10) {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Too few chi2pid entries for "
                      << (chi2Pions->GetEntries() < 10 ? "pions" : "kaons") << ": "
                      << (chi2Pions->GetEntries() < 10 ? chi2Pions->GetEntries() : chi2Kaons->GetEntries()) << "\n";
            csvFileBefore << pLow << "-" << pHigh << ",N/A,N/A\n";
            csvFilePionCut << pLow << "-" << pHigh << ",N/A,N/A\n";
            csvFileSigmaBin3 << pLow << "-" << pHigh << ",N/A,N/A\n";
            delete chi2Pions;
            delete chi2Kaons;
            continue;
        }

        // Set colors
        chi2Pions->SetLineColor(kBlue);
        chi2Kaons->SetLineColor(kGreen);

        // Fit chi2pid histograms
        double pionChi2Min, pionChi2Max;
        auto [pionAmp, pionMean, pionSigma] = chi2FitParams[i];
        auto [kaonAmp, kaonMean, kaonSigma] = chi2FitParams[i];
        TF1* fitPion = new TF1("gausFit", "gaus", chi2FitRanges[i].first, chi2FitRanges[i].second);
        fitPion->SetParameters(pionAmp, pionMean, pionSigma);
        chi2Pions->Fit(fitPion, "R", "", chi2FitRanges[i].first, chi2FitRanges[i].second);
        pionChi2Min = fitPion->GetParameter(1) - 3 * fitPion->GetParameter(2);
        pionChi2Max = fitPion->GetParameter(1) + 3 * fitPion->GetParameter(2);
        if (i == 2) { // Bin [1.6, 1.9)
            pionSigmaBin3 = fitPion->GetParameter(2);
            std::cout << "Bin [1.6-1.9): Pion sigma = " << pionSigmaBin3 << "\n";
        }
        TF1* fitKaon = new TF1("gausFit", "gaus", chi2FitRanges[i].first, chi2FitRanges[i].second);
        fitKaon->SetParameters(kaonAmp, kaonMean, kaonSigma);
        chi2Kaons->Fit(fitKaon, "R", "", chi2FitRanges[i].first, chi2FitRanges[i].second);
        std::cout << "Bin [" << pLow << "-" << pHigh << "): Pions chi2 range: [" 
                  << pionChi2Min << ", " << pionChi2Max << "]\n";

        // Plot chi2pid histograms
        canvasChi2->Clear();
        canvasChi2->Divide(2, 1);
        canvasChi2->cd(1);
        chi2Pions->Draw("HIST");
        if (chi2Pions->GetFunction("gausFit")) {
            chi2Pions->GetFunction("gausFit")->SetLineColor(kBlue);
            chi2Pions->GetFunction("gausFit")->Draw("SAME");
        }
        TLegend* legPions = new TLegend(0.1, 0.75, 0.3, 0.9);
        legPions->AddEntry(chi2Pions, "Pions", "l");
        legPions->SetBorderSize(0);
        legPions->SetFillStyle(0);
        legPions->SetTextSize(0.025);
        legPions->Draw();
        canvasChi2->cd(2);
        chi2Kaons->Draw("HIST");
        if (chi2Kaons->GetFunction("gausFit")) {
            chi2Kaons->GetFunction("gausFit")->SetLineColor(kGreen);
            chi2Kaons->GetFunction("gausFit")->Draw("SAME");
        }
        TLegend* legKaons = new TLegend(0.1, 0.75, 0.3, 0.9);
        legKaons->AddEntry(chi2Kaons, "Kaons", "l");
        legKaons->SetBorderSize(0);
        legKaons->SetFillStyle(0);
        legKaons->SetTextSize(0.025);
        legKaons->Draw();
        canvasChi2->Update();
        canvasChi2->Print("output/pdf/chi2pid_fits_linear.pdf");
        delete legPions;
        delete legKaons;
        delete chi2Pions;
        delete chi2Kaons;

        // Define beta cuts for "before chi2pid cuts" section
        double pionLeftBefore, pionRightBefore, kaonRightBefore;
        if (i == 0) { // Bin [1.0-1.3)
            pionLeftBefore = 0.9880;
            pionRightBefore = 1.0120;
            kaonRightBefore = 0.0;
        } else if (i == 1) { // Bin [1.3-1.6)
            pionLeftBefore = 0.9755;
            pionRightBefore = 1.0135;
            kaonRightBefore = 0.9721;
        } else if (i == 2) { // Bin [1.6-1.9)
            pionLeftBefore = 0.98;
            pionRightBefore = 1.0140;
            kaonRightBefore = 0.986;
        } else if (i == 3) { // Bin [1.9-2.2)
            pionLeftBefore = 0.9815;
            pionRightBefore = 1.0145;
            kaonRightBefore = 0.992;
        } else if (i == 4) { // Bin [2.2-2.5)
            pionLeftBefore = 0.9880;
            pionRightBefore = 1.0141;
            kaonRightBefore = 0.9971;
        } else if (i == 5) { // Bin [2.5-2.8)
            pionLeftBefore = 0.9848;
            pionRightBefore = 1.013;
            kaonRightBefore = 1.0;
        } else if (i == 6) { // Bin [2.8-3.1)
            pionLeftBefore = 0.986;
            pionRightBefore = 1.0128;
            kaonRightBefore = 1.003;
        } else if (i == 7) { // Bin [3.1-3.4)
            pionLeftBefore = 0.9868;
            pionRightBefore = 1.0121;
            kaonRightBefore = 1.0048;
        } else if (i == 8) { // Bin [3.4-3.7)
            pionLeftBefore = 0.9879;
            pionRightBefore = 1.0120;
            kaonRightBefore = 1.0049;
        } else if (i == 9) { // Bin [3.7-4.0)
            pionLeftBefore = 0.9880;
            pionRightBefore = 1.0120;
            kaonRightBefore = 1.012;
        } else {
            std::cerr << "Warning: Using default beta cuts for bin [" << pLow << "-" << pHigh << ") in before chi2pid section\n";
            pionLeftBefore = 0.9880;
            pionRightBefore = 1.0120;
            kaonRightBefore = 1.0;
        }
        if (pionLeftBefore >= pionRightBefore) {
            std::cerr << "Error: Invalid pion cuts for bin [" << pLow << "-" << pHigh << ") in before chi2pid section\n";
            csvFileBefore << pLow << "-" << pHigh << ",N/A,N/A\n";
            continue;
        }



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
            if (p >= pLow && p < pHigh ) {
                betaPionsBefore->Fill(beta);
            }
        }
        treeKaons->SetBranchAddress("p", &p);
        treeKaons->SetBranchAddress("beta", &beta);
        for (Long64_t j = 0; j < treeKaons->GetEntries(); j++) {
            treeKaons->GetEntry(j);
            if (p >= pLow && p < pHigh ) {
                betaKaonsBefore->Fill(beta);
            }
        }

        // Check for sufficient entries
        if (betaPionsBefore->GetEntries() < 10 || betaKaonsBefore->GetEntries() < 10) {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Too few beta entries for "
                      << (betaPionsBefore->GetEntries() < 10 ? "pions" : "kaons") << ": "
                      << (betaPionsBefore->GetEntries() < 10 ? betaPionsBefore->GetEntries() : betaKaonsBefore->GetEntries()) << "\n";
            csvFileBefore << pLow << "-" << pHigh << ",N/A,N/A\n";
            csvFilePionCut << pLow << "-" << pHigh << ",N/A,N/A\n";
            csvFileSigmaBin3 << pLow << "-" << pHigh << ",N/A,N/A\n";
            delete betaPionsBefore;
            delete betaKaonsBefore;
            continue;
        }

        // Fit beta histograms (before cut) for contamination
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


        bool logScale = false;
        double yMin = 0.0;
        double yMax = -1.0;
        // Plot beta histograms (before cut)
        canvasBeta->Clear();
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
            yMax = std::max(yMax, betaKaonsBefore->GetMaximum());
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
        canvasBeta->SetGrid();
        canvasBeta->Update();
        canvasBeta->Print("output/pdf/beta_before_cut_linear.pdf");
        delete leg;

        // Calculate contamination (before chi2pid cuts)
        double c1Before, c2Before;
        std::map<std::string, double> pFitBefore = {{"gaus_mean", pMeanBefore}, {"gaus_sigma", pSigmaBefore}, {"gaus_constant", pConstBefore}};
        std::map<std::string, double> kFitBefore = {{"gaus_mean", kMeanBefore}, {"gaus_sigma", kSigmaBefore}, {"gaus_constant", kConstBefore}};
        double contaminationBefore = calculateOverlapContamination(pFitBefore, kFitBefore, pionLeftBefore, pionRightBefore, kaonRightBefore, c1Before, c2Before);
        double kaonLeftBefore = -1.0;
        if (c1Before >= kMeanBefore && c1Before <= pMeanBefore) {
            kaonLeftBefore = c1Before;
        } else if (c2Before >= kMeanBefore && c2Before <= pMeanBefore) {
            kaonLeftBefore = c2Before;
        } else if (contaminationBefore >= 0) {
            // If no intersection lies between kMeanBefore and pMeanBefore, choose the point closest to kaon mean
            kaonLeftBefore = (std::abs(c1Before - kMeanBefore) < std::abs(c2Before - kMeanBefore)) ? c1Before : c2Before;
            std::cout << "No intersection between kMeanBefore = " << kMeanBefore << " and pMeanBefore = " << pMeanBefore 
                    << ", selecting point closest to kaon peak: " << kaonLeftBefore << std::endl;
        }

        // Save contamination (before chi2pid cuts)
        if (contaminationBefore < 0 || kaonLeftBefore == -1.0) {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Contamination (before chi2pid) calculation failed\n";
            csvFileBefore << pLow << "-" << pHigh << ",N/A,N/A\n";
        } else {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Before chi2pid - Pion range: [" << pionLeftBefore << ", " << pionRightBefore 
                      << "], Kaon range: [" << kaonLeftBefore << ", " << kaonRightBefore << "], Contamination: " 
                      << contaminationBefore << "%\n";
            csvFileBefore << pLow << "-" << pHigh << "," << kaonLeftBefore << "," << contaminationBefore << "\n";
        }
        delete betaPionsBefore;
        delete betaKaonsBefore;

        /* // Define beta cuts for "pion chi2pid cut only" section
        double pionLeftPionCut, pionRightPionCut, kaonRightPionCut;
        if (i == 0) { // Bin [1.0-1.3)
            pionLeftPionCut = 0.9880;
            pionRightPionCut = 1.0120;
            kaonRightPionCut = 0.0;
        } else if (i == 1) { // Bin [1.3-1.6)
            pionLeftPionCut = 0.9755;
            pionRightPionCut = 1.0135;
            kaonRightPionCut = 0.9721;
        } else if (i == 2) { // Bin [1.6-1.9)
            pionLeftPionCut = 0.98;
            pionRightPionCut = 1.0140;
            kaonRightPionCut = 0.986;
        } else if (i == 3) { // Bin [1.9-2.2)
            pionLeftPionCut = 0.9815;
            pionRightPionCut = 1.0145;
            kaonRightPionCut = 0.992;
        } else if (i == 4) { // Bin [2.2-2.5)
            pionLeftPionCut = 0.9880;
            pionRightPionCut = 1.0141;
            kaonRightPionCut = 0.9971;
        } else if (i == 5) { // Bin [2.5-2.8)
            pionLeftPionCut = 0.9848;
            pionRightPionCut = 1.013;
            kaonRightPionCut = 1.0;
        } else if (i == 6) { // Bin [2.8-3.1)
            pionLeftPionCut = 0.986;
            pionRightPionCut = 1.0128;
            kaonRightPionCut = 1.003;
        } else if (i == 7) { // Bin [3.1-3.4)
            pionLeftPionCut = 0.9868;
            pionRightPionCut = 1.0121;
            kaonRightPionCut = 1.0048;
        } else if (i == 8) { // Bin [3.4-3.7)
            pionLeftPionCut = 0.9879;
            pionRightPionCut = 1.0120;
            kaonRightPionCut = 1.0049;
        } else if (i == 9) { // Bin [3.7-4.0)
            pionLeftPionCut = 0.9880;
            pionRightPionCut = 1.0120;
            kaonRightPionCut = 1.012;
        } else {
            std::cerr << "Warning: Using default beta cuts for bin [" << pLow << "-" << pHigh << ") in pion chi2pid cut only section\n";
            pionLeftPionCut = 0.9880;
            pionRightPionCut = 1.0120;
        }
        if (pionLeftPionCut >= pionRightPionCut) {
            std::cerr << "Error: Invalid pion cuts for bin [" << pLow << "-" << pHigh << ") in pion chi2pid cut only section\n";
            csvFilePionCut << pLow << "-" << pHigh << ",N/A,N/A\n";
            continue;
        }

        // Create beta histograms (pion chi2pid cut only)
        TH1F* betaPionsAfterPionCut = new TH1F(TString::Format("beta_pions_after_pioncut_%d", i),
                                               TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                               70, betaHistRanges[i].first, betaHistRanges[i].second);
        TH1F* betaKaonsAfterPionCut = new TH1F(TString::Format("beta_kaons_after_pioncut_%d", i),
                                               TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                               70, betaHistRanges[i].first, betaHistRanges[i].second);

        betaPionsAfterPionCut->Sumw2();
        betaKaonsAfterPionCut->Sumw2();

        // Fill beta histograms (pion chi2pid cut only)
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("chi2pid", &chi2pid);
        treePions->SetBranchAddress("beta", &beta);
        for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh && chi2pid >= pionChi2Min && chi2pid <= pionChi2Max &&
                beta >= betaHistRanges[i].first && beta <= betaHistRanges[i].second) {
                betaPionsAfterPionCut->Fill(beta);
            }
        }
        treeKaons->SetBranchAddress("p", &p);
        treeKaons->SetBranchAddress("beta", &beta);
        for (Long64_t j = 0; j < treeKaons->GetEntries(); j++) {
            treeKaons->GetEntry(j);
            if (p >= pLow && p < pHigh && beta >= betaHistRanges[i].first && beta <= betaHistRanges[i].second) {
                betaKaonsAfterPionCut->Fill(beta);
            }
        }

        // Check for sufficient entries
        if (betaPionsAfterPionCut->GetEntries() < 10 && betaKaonsAfterPionCut->GetEntries() < 10) {
            csvFilePionCut << pLow << "-" << pHigh << ",N/A,N/A\n";
            delete betaPionsAfterPionCut;
            delete betaKaonsAfterPionCut;
            continue;
        }

        // Fit beta histograms (pion chi2pid cut only)
        double pMeanPionCut = 0, pSigmaPionCut = 0, pConstPionCut = 0, kMeanPionCut = 0, kSigmaPionCut = 0, kConstPionCut = 0;
        if (betaPionsAfterPionCut->GetEntries() >= 10) {
            TF1* gausPionPionCut = new TF1(TString::Format("gaus_%s", betaPionsAfterPionCut->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            gausPionPionCut->SetParameters(betaPionsAfterPionCut->GetMaximum(), 1.0, 0.005);
            betaPionsAfterPionCut->Fit(gausPionPionCut, "R", "", betaFitRanges[i].first, betaFitRanges[i].second);
            pMeanPionCut = gausPionPionCut->GetParameter(1);
            pSigmaPionCut = gausPionPionCut->GetParameter(2);
            pConstPionCut = gausPionPionCut->GetParameter(0);
            delete gausPionPionCut;
        }
        if (betaKaonsAfterPionCut->GetEntries() >= 10) {
            TF1* gausKaonPionCut = new TF1(TString::Format("gaus_%s", betaKaonsAfterPionCut->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            if (i == 9) { // Bin [3.7-4.0) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.1);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.986, 0.9945);
            } else if (i == 8) { // Bin [3.4-3.7) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.1);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.988, 0.9945);
            } else if (i == 10) { // Bin [4.0-4.3) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.1);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.987, 0.995);
            } else if (i == 11) { // Bin [4.3-4.6) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.996, 0.012);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.992, 0.996);
            } else if (i == 12) { // Bin [4.6-4.9) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.012);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.991, 0.996);
            } else if (i == 13) { // Bin [4.9-5.2) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.012);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.991, 0.995);
            } else if (i == 14) { // Bin [5.2-5.5) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.1);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.991, 0.996);
            } else if (i == 15) { // Bin [5.5-5.8) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.1);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.991, 0.996);
            } else if (i == 16) { // Bin [5.8-6.1) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.1);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.991, 0.997);
            } else if (i == 17) { // Bin [6.1-6.4) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.1);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.991, 0.997);
            } else if (i == 18) { // Bin [6.4-6.7) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.1);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.982, 0.997);
            } else if (i == 5) { // Bin [2.5-2.8) GeV/c
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.974, 0.992);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.974, 0.992);
            } else {
                gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.99, 0.005);
                betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R");
            }
            kMeanPionCut = gausKaonPionCut->GetParameter(1);
            kSigmaPionCut = gausKaonPionCut->GetParameter(2);
            kConstPionCut = gausKaonPionCut->GetParameter(0);
            delete gausKaonPionCut;
        }

        // Calculate contamination (pion chi2pid cut only)
        double c1PionCut, c2PionCut;
        std::map<std::string, double> pFitPionCut = {{"gaus_mean", pMeanPionCut}, {"gaus_sigma", pSigmaPionCut}, {"gaus_constant", pConstPionCut}};
        std::map<std::string, double> kFitPionCut = {{"gaus_mean", kMeanPionCut}, {"gaus_sigma", kSigmaPionCut}, {"gaus_constant", kConstPionCut}};
        double contaminationPionCut = calculateOverlapContamination(pFitPionCut, kFitPionCut, pionLeftPionCut, pionRightPionCut, kaonRightPionCut, c1PionCut, c2PionCut);
        double kaonLeftPionCut = -1.0;
        if (c1PionCut >= kMeanPionCut && c1PionCut <= pMeanPionCut) {
            kaonLeftPionCut = c1PionCut;
        } else if (c2PionCut >= kMeanPionCut && c2PionCut <= pMeanPionCut) {
            kaonLeftPionCut = c2PionCut;
        } else if (contaminationPionCut >= 0) {
            // If no intersection lies between kMeanPionCut and pMeanPionCut, choose the point closest to kaon mean
            kaonLeftPionCut = (std::abs(c1PionCut - kMeanPionCut) < std::abs(c2PionCut - kMeanPionCut)) ? c1PionCut : c2PionCut;
            std::cout << "No intersection between kMeanPionCut = " << kMeanPionCut << " and pMeanPionCut = " << pMeanPionCut 
                    << ", selecting point closest to kaon peak: " << kaonLeftPionCut << std::endl;
        }
        // Plot beta histograms (pion chi2pid cut only)
        canvasBetaPionCut->Clear();
        if (betaPionsAfterPionCut->GetEntries() >= 10) {
            betaPionsAfterPionCut->SetMarkerStyle(20);
            betaPionsAfterPionCut->SetMarkerSize(1);
            betaPionsAfterPionCut->SetMarkerColor(kBlue);
            betaPionsAfterPionCut->Draw("P");
            yMax = betaPionsAfterPionCut->GetMaximum();
            if (auto* fit = betaPionsAfterPionCut->GetFunction(TString::Format("gaus_%s", betaPionsAfterPionCut->GetName()))) {
                fit->SetLineColor(kBlue);
                fit->SetLineStyle(kDashed);
                fit->Draw("SAME");
            }
        }
        if (betaKaonsAfterPionCut->GetEntries() >= 10) {
            betaKaonsAfterPionCut->SetMarkerStyle(20);
            betaKaonsAfterPionCut->SetMarkerSize(1);
            betaKaonsAfterPionCut->SetMarkerColor(kGreen);
            betaKaonsAfterPionCut->Draw(betaPionsAfterPionCut->GetEntries() >= 10 ? "P SAME" : "P");
            yMax = std::max(yMax, betaKaonsAfterPionCut->GetMaximum());
            if (auto* fit = betaKaonsAfterPionCut->GetFunction(TString::Format("gaus_%s", betaKaonsAfterPionCut->GetName()))) {
                fit->SetLineColor(kGreen);
                fit->SetLineStyle(kDashed);
                fit->Draw("SAME");
                fit->SetRange(0.95, 1.03);
            }
        }
        yMax *= 1.2;
        if (betaPionsAfterPionCut->GetEntries() >= 10 || betaKaonsAfterPionCut->GetEntries() >= 10) {
            (betaPionsAfterPionCut->GetEntries() >= 10 ? betaPionsAfterPionCut : betaKaonsAfterPionCut)->GetYaxis()->SetRangeUser(yMin, yMax);
        }
        if (kaonLeftPionCut >= 0) {
            TLine* line = new TLine(kaonLeftPionCut, yMin, kaonLeftPionCut, yMax);
            line->SetLineColor(kRed);
            line->SetLineStyle(kSolid);
            line->SetLineWidth(2);
            line->Draw("SAME");
            delete line;
        }
        leg = new TLegend(0.15, 0.75, 0.35, 0.88);
        if (betaPionsAfterPionCut->GetEntries() >= 10) leg->AddEntry(betaPionsAfterPionCut, "Pions", "p");
        if (betaKaonsAfterPionCut->GetEntries() >= 10) leg->AddEntry(betaKaonsAfterPionCut, "Kaons", "p");
        leg->SetTextSize(0.03);
        leg->Draw();
        canvasBetaPionCut->SetGrid();
        canvasBetaPionCut->Update();
        canvasBetaPionCut->Print("output/pdf/beta_pion_cut_only_linear.pdf");
        delete leg;

        // Save contamination (pion chi2pid cut only)
        if (contaminationPionCut < 0 || kaonLeftPionCut == -1.0) {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Contamination (pion cut only) calculation failed\n";
            csvFilePionCut << pLow << "-" << pHigh << ",N/A,N/A\n";
        } else {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Pion cut only - Pion range: [" << pionLeftPionCut << ", " << pionRightPionCut 
                      << "], Kaon range: [" << kaonLeftPionCut << ", " << kaonRightPionCut << "], Contamination: " 
                      << contaminationPionCut << "%\n";
            csvFilePionCut << pLow << "-" << pHigh << "," << kaonLeftPionCut << "," << contaminationPionCut << "\n";
        }
        delete betaPionsAfterPionCut;
        delete betaKaonsAfterPionCut;

        // Skip bins above 4.9 GeV/c for the new case
        if (pLow >= 4.9) {
            csvFileSigmaBin3 << pLow << "-" << pHigh << ",N/A,N/A\n";
            continue;
        }

       // New case: Beta with chi2pid cut using sigma from bin 3 ([1.6, 1.9) GeV/c)

       double pionLeftBin3, pionRightBin3, kaonRightBin3;
        if (i == 0) { // Bin [1.0-1.3)
            pionLeftBin3 = 0.9880;
            pionRightBin3 = 1.0120;
            kaonRightBin3 = 0.0;
        } else if (i == 1) { // Bin [1.3-1.6)
            pionLeftBin3 = 0.9755;
            pionRightBin3 = 1.0135;
            kaonRightBin3 = 0.9721;
        } else if (i == 2) { // Bin [1.6-1.9)
            pionLeftBin3 = 0.98;
            pionRightBin3 = 1.0140;
            kaonRightBin3 = 0.986;
        } else if (i == 3) { // Bin [1.9-2.2)
            pionLeftBin3 = 0.9815;
            pionRightBin3 = 1.0145;
            kaonRightBin3 = 0.992;
        } else if (i == 4) { // Bin [2.2-2.5)
            pionLeftBin3 = 0.9880;
            pionRightBin3 = 1.0141;
            kaonRightBin3 = 0.9971;
        } else if (i == 5) { // Bin [2.5-2.8)
            pionLeftBin3 = 0.9848;
            pionRightBin3 = 1.013;
            kaonRightBin3 = 1.0;
        } else if (i == 6) { // Bin [2.8-3.1)
            pionLeftBin3 = 0.986;
            pionRightBin3 = 1.0128;
            kaonRightBin3 = 1.003;
        } else if (i == 7) { // Bin [3.1-3.4)
            pionLeftBin3 = 0.9868;
            pionRightBin3 = 1.0121;
            kaonRightBin3 = 1.0048;
        } else if (i == 8) { // Bin [3.4-3.7)
            pionLeftBin3 = 0.9879;
            pionRightBin3 = 1.0120;
            kaonRightBin3 = 1.0049;
        } else if (i == 9) { // Bin [3.7-4.0)
            pionLeftBin3 = 0.9880;
            pionRightBin3 = 1.0120;
            kaonRightBin3 = 1.012;
        } else {
            std::cerr << "Warning: Using default beta cuts for bin [" << pLow << "-" << pHigh << ") in pion chi2pid cut only section\n";
            pionLeftBin3 = 0.9880;
            pionRightBin3 = 1.0120;
        }
        if (pionLeftBin3 >= pionRightBin3) {
            std::cerr << "Error: Invalid pion cuts for bin [" << pLow << "-" << pHigh << ") in pion chi2pid cut only section\n";
            csvFilePionCut << pLow << "-" << pHigh << ",N/A,N/A\n";
            continue;
        }
        // Positive side: 3-sigma for p < 2.8 GeV/c, 2-sigma for 2.8 <= p < 4.9 GeV/c
        // Negative side: 3-sigma for all bins
        double nSigmaPositive = (pLow < 2.8) ? 3.0 : 2.0;
        double nSigmaNegative = 3.0;
        double pionChi2MinBin3 = fitPion->GetParameter(1) - nSigmaNegative * pionSigmaBin3;
        double pionChi2MaxBin3 = fitPion->GetParameter(1) + nSigmaPositive * pionSigmaBin3;

        // Store chi2pid cut values in a text file
        static std::ofstream chi2pidFile;
        if (!chi2pidFile.is_open()) {
            gSystem->MakeDirectory("output");
            gSystem->MakeDirectory("output/txt");
            chi2pidFile.open("output/txt/chi2pid_cuts_sigma_bin3.txt", std::ios::out);
            if (!chi2pidFile.is_open()) {
                std::cerr << "Error: Could not open output/txt/chi2pid_cuts_sigma_bin3.txt" << std::endl;
                continue;
            }
            chi2pidFile << "Momentum_Range,Chi2pid_Negative_Cut,Chi2pid_Positive_Cut\n";
        }
        chi2pidFile << pLow << "-" << pHigh << "," << pionChi2MinBin3 << "," << pionChi2MaxBin3 << "\n";

        // Create beta histograms (sigma from bin 3)
        TH1F* betaPionsSigmaBin3 = new TH1F(TString::Format("beta_pions_sigma_bin3_%d", i),
                                            TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                            70, betaHistRanges[i].first, betaHistRanges[i].second);
        TH1F* betaKaonsSigmaBin3 = new TH1F(TString::Format("beta_kaons_sigma_bin3_%d", i),
                                            TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                            70, betaHistRanges[i].first, betaHistRanges[i].second);

        betaPionsSigmaBin3->Sumw2();
        betaKaonsSigmaBin3->Sumw2();

        // Fill beta histograms (sigma from bin 3)
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("chi2pid", &chi2pid);
        treePions->SetBranchAddress("beta", &beta);
        for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh && chi2pid >= pionChi2MinBin3 && chi2pid <= pionChi2MaxBin3 ) {
                betaPionsSigmaBin3->Fill(beta);
            }
        }
        treeKaons->SetBranchAddress("p", &p);
        treeKaons->SetBranchAddress("beta", &beta);
        for (Long64_t j = 0; j < treeKaons->GetEntries(); j++) {
            treeKaons->GetEntry(j);
            if (p >= pLow && p < pHigh ) {
                betaKaonsSigmaBin3->Fill(beta);
            }
        }

        // Check for sufficient entries
        if (betaPionsSigmaBin3->GetEntries() < 10 && betaKaonsSigmaBin3->GetEntries() < 10) {
            csvFileSigmaBin3 << pLow << "-" << pHigh << ",N/A,N/A\n";
            delete betaPionsSigmaBin3;
            delete betaKaonsSigmaBin3;
            continue;
        }

        // Fit beta histograms (sigma from bin 3)
        double pMeanSigmaBin3 = 0, pSigmaSigmaBin3 = 0, pConstSigmaBin3 = 0, kMeanSigmaBin3 = 0, kSigmaSigmaBin3 = 0, kConstSigmaBin3 = 0;
        if (betaPionsSigmaBin3->GetEntries() >= 10) {
            TF1* gausPionSigmaBin3 = new TF1(TString::Format("gaus_%s", betaPionsSigmaBin3->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            gausPionSigmaBin3->SetParameters(betaPionsSigmaBin3->GetMaximum(), 1.0, 0.005);
            betaPionsSigmaBin3->Fit(gausPionSigmaBin3, "R", "", betaFitRanges[i].first, betaFitRanges[i].second);
            pMeanSigmaBin3 = gausPionSigmaBin3->GetParameter(1);
            pSigmaSigmaBin3 = gausPionSigmaBin3->GetParameter(2);
            pConstSigmaBin3 = gausPionSigmaBin3->GetParameter(0);
            delete gausPionSigmaBin3;
        }
        if (betaKaonsSigmaBin3->GetEntries() >= 10) {
            TF1* gausKaonSigmaBin3 = new TF1(TString::Format("gaus_%s", betaKaonsSigmaBin3->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            if (i == 9) { // Bin [3.7-4.0) GeV/c
                gausKaonSigmaBin3->SetParameters(betaKaonsSigmaBin3->GetMaximum(), 0.994, 0.1);
                betaKaonsSigmaBin3->Fit(gausKaonSigmaBin3, "R", "", 0.986, 0.9945);
            } else if (i == 8) { // Bin [3.4-3.7) GeV/c
                gausKaonSigmaBin3->SetParameters(betaKaonsSigmaBin3->GetMaximum(), 0.994, 0.1);
                betaKaonsSigmaBin3->Fit(gausKaonSigmaBin3, "R", "", 0.988, 0.9945);
            } else if (i == 10) { // Bin [4.0-4.3) GeV/c
                gausKaonSigmaBin3->SetParameters(betaKaonsSigmaBin3->GetMaximum(), 0.994, 0.1);
                betaKaonsSigmaBin3->Fit(gausKaonSigmaBin3, "R", "", 0.987, 0.995);
            } else if (i == 11) { // Bin [4.3-4.6) GeV/c
                gausKaonSigmaBin3->SetParameters(betaKaonsSigmaBin3->GetMaximum(), 0.996, 0.012);
                betaKaonsSigmaBin3->Fit(gausKaonSigmaBin3, "R", "", 0.992, 0.996);
            } else if (i == 12) { // Bin [4.6-4.9) GeV/c
                gausKaonSigmaBin3->SetParameters(betaKaonsSigmaBin3->GetMaximum(), 0.994, 0.012);
                betaKaonsSigmaBin3->Fit(gausKaonSigmaBin3, "R", "", 0.991, 0.996);
            } else if (i == 13) { // Bin [4.9-5.2) GeV/c
                gausKaonSigmaBin3->SetParameters(betaKaonsSigmaBin3->GetMaximum(), 0.994, 0.012);
                betaKaonsSigmaBin3->Fit(gausKaonSigmaBin3, "R", "", 0.991, 0.995);
            } else if (i == 14) { // Bin [5.2-5.5) GeV/c
                gausKaonSigmaBin3->SetParameters(betaKaonsSigmaBin3->GetMaximum(), 0.994, 0.1);
                betaKaonsSigmaBin3->Fit(gausKaonSigmaBin3, "R", "", 0.991, 0.996);
            } else if (i == 15) { // Bin [5.5-5.8) GeV/c
                gausKaonSigmaBin3->SetParameters(betaKaonsSigmaBin3->GetMaximum(), 0.994, 0.1);
                betaKaonsSigmaBin3->Fit(gausKaonSigmaBin3, "R", "", 0.991, 0.996);
            } else if (i == 16) { // Bin [5.8-6.1) GeV/c
                gausKaonSigmaBin3->SetParameters(betaKaonsSigmaBin3->GetMaximum(), 0.994, 0.1);
                betaKaonsSigmaBin3->Fit(gausKaonSigmaBin3, "R", "", 0.991, 0.997);
            } else if (i == 17) { // Bin [6.1-6.4) GeV/c
                gausKaonSigmaBin3->SetParameters(betaKaonsSigmaBin3->GetMaximum(), 0.994, 0.1);
                betaKaonsSigmaBin3->Fit(gausKaonSigmaBin3, "R", "", 0.991, 0.997);
            } else if (i == 18) { // Bin [6.4-6.7) GeV/c
                gausKaonSigmaBin3->SetParameters(betaKaonsSigmaBin3->GetMaximum(), 0.994, 0.1);
                betaKaonsSigmaBin3->Fit(gausKaonSigmaBin3, "R", "", 0.982, 0.997);
            } else if (i == 5) { // Bin [2.5-2.8) GeV/c
                gausKaonSigmaBin3->SetParameters(betaKaonsSigmaBin3->GetMaximum(), 0.974, 0.992);
                betaKaonsSigmaBin3->Fit(gausKaonSigmaBin3, "R", "", 0.974, 0.992);
            } else {
                gausKaonSigmaBin3->SetParameters(betaKaonsSigmaBin3->GetMaximum(), 0.99, 0.005);
                betaKaonsSigmaBin3->Fit(gausKaonSigmaBin3, "R");
            }
            kMeanSigmaBin3 = gausKaonSigmaBin3->GetParameter(1);
            kSigmaSigmaBin3 = gausKaonSigmaBin3->GetParameter(2);
            kConstSigmaBin3 = gausKaonSigmaBin3->GetParameter(0);
            delete gausKaonSigmaBin3;
        }

        // Calculate contamination (sigma from bin 3)
        double c1SigmaBin3, c2SigmaBin3;
        std::map<std::string, double> pFitSigmaBin3 = {{"gaus_mean", pMeanSigmaBin3}, {"gaus_sigma", pSigmaSigmaBin3}, {"gaus_constant", pConstSigmaBin3}};
        std::map<std::string, double> kFitSigmaBin3 = {{"gaus_mean", kMeanSigmaBin3}, {"gaus_sigma", kSigmaSigmaBin3}, {"gaus_constant", kConstSigmaBin3}};
        double contaminationSigmaBin3 = calculateOverlapContamination(pFitSigmaBin3, kFitSigmaBin3, pionLeftBin3, pionRightBin3, kaonRightBin3, c1SigmaBin3, c2SigmaBin3);
        double kaonLeftSigmaBin3 = -1.0;
        if (c1SigmaBin3 >= kMeanSigmaBin3 && c1SigmaBin3 <= pMeanSigmaBin3) {
            kaonLeftSigmaBin3 = c1SigmaBin3;
        } else if (c2SigmaBin3 >= kMeanSigmaBin3 && c2SigmaBin3 <= pMeanSigmaBin3) {
            kaonLeftSigmaBin3 = c2SigmaBin3;
        } else if (contaminationSigmaBin3 >= 0) {
            // If no intersection lies between kMeanSigmaBin3 and pMeanSigmaBin3, choose the point closest to kaon mean
            kaonLeftSigmaBin3 = (std::abs(c1SigmaBin3 - kMeanSigmaBin3) < std::abs(c2SigmaBin3 - kMeanSigmaBin3)) ? c1SigmaBin3 : c2SigmaBin3;
            std::cout << "No intersection between kMeanSigmaBin3 = " << kMeanSigmaBin3 << " and pMeanSigmaBin3 = " << pMeanSigmaBin3 
                    << ", selecting point closest to kaon peak: " << kaonLeftSigmaBin3 << std::endl;
        }
        // Plot beta histograms (sigma from bin 3)
        canvasBetaSigmaBin3->Clear();
        logScale = false;
        yMin = 0.0;
        yMax = -1.0;
        if (betaPionsSigmaBin3->GetEntries() >= 10) {
            betaPionsSigmaBin3->SetMarkerStyle(20);
            betaPionsSigmaBin3->SetMarkerSize(1);
            betaPionsSigmaBin3->SetMarkerColor(kBlue);
            betaPionsSigmaBin3->Draw("P");
            yMax = betaPionsSigmaBin3->GetMaximum();
            if (auto* fit = betaPionsSigmaBin3->GetFunction(TString::Format("gaus_%s", betaPionsSigmaBin3->GetName()))) {
                fit->SetLineColor(kBlue);
                fit->SetLineStyle(kDashed);
                fit->Draw("SAME");
            }
        }
        if (betaKaonsSigmaBin3->GetEntries() >= 10) {
            betaKaonsSigmaBin3->SetMarkerStyle(20);
            betaKaonsSigmaBin3->SetMarkerSize(1);
            betaKaonsSigmaBin3->SetMarkerColor(kGreen);
            betaKaonsSigmaBin3->Draw(betaPionsSigmaBin3->GetEntries() >= 10 ? "P SAME" : "P");
            yMax = std::max(yMax, betaKaonsSigmaBin3->GetMaximum());
            if (auto* fit = betaKaonsSigmaBin3->GetFunction(TString::Format("gaus_%s", betaKaonsSigmaBin3->GetName()))) {
                fit->SetLineColor(kGreen);
                fit->SetLineStyle(kDashed);
                fit->Draw("SAME");
                fit->SetRange(0.95, 1.03);
            }
        }
        yMax *= 1.2;
        if (betaPionsSigmaBin3->GetEntries() >= 10 || betaKaonsSigmaBin3->GetEntries() >= 10) {
            (betaPionsSigmaBin3->GetEntries() >= 10 ? betaPionsSigmaBin3 : betaKaonsSigmaBin3)->GetYaxis()->SetRangeUser(yMin, yMax);
        }
        if (kaonLeftSigmaBin3 >= 0) {
            TLine* line = new TLine(kaonLeftSigmaBin3, yMin, kaonLeftSigmaBin3, yMax);
            line->SetLineColor(kRed);
            line->SetLineStyle(kSolid);
            line->SetLineWidth(2);
            line->Draw("SAME");
            delete line;
        }
        leg = new TLegend(0.15, 0.75, 0.35, 0.88);
        if (betaPionsSigmaBin3->GetEntries() >= 10) leg->AddEntry(betaPionsSigmaBin3, "Pions", "p");
        if (betaKaonsSigmaBin3->GetEntries() >= 10) leg->AddEntry(betaKaonsSigmaBin3, "Kaons", "p");
        leg->SetTextSize(0.03);
        leg->Draw();
        canvasBetaSigmaBin3->SetGrid();
        canvasBetaSigmaBin3->Update();
        canvasBetaSigmaBin3->Print("output/pdf/beta_sigma_from_bin3_linear.pdf");
        delete leg;

        // Save contamination (sigma from bin 3)
        if (contaminationSigmaBin3 < 0 || kaonLeftSigmaBin3 == -1.0) {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Contamination (sigma from bin 3) calculation failed\n";
            csvFileSigmaBin3 << pLow << "-" << pHigh << ",N/A,N/A\n";
        } else {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Sigma from bin 3 - Pion range: [" << pionLeftPionCut << ", " << pionRightPionCut 
                      << "], Kaon range: [" << kaonLeftSigmaBin3 << ", " << kaonRightPionCut << "], Contamination: " 
                      << contaminationSigmaBin3 << "%\n";
            csvFileSigmaBin3 << pLow << "-" << pHigh << "," << kaonLeftSigmaBin3 << "," << contaminationSigmaBin3 << "\n";
        }
        delete betaPionsSigmaBin3;
        delete betaKaonsSigmaBin3; */
    }

    // Close PDFs and CSVs
    canvasChi2->Print("output/pdf/chi2pid_fits_linear.pdf]");
    canvasBeta->Print("output/pdf/beta_before_cut_linear.pdf]");
    /* canvasBetaPionCut->Print("output/pdf/beta_pion_cut_only_linear.pdf]");
    canvasBetaSigmaBin3->Print("output/pdf/beta_sigma_from_bin3_linear.pdf]"); */
    csvFileBefore.close();
    /* csvFilePionCut.close();
    csvFileSigmaBin3.close(); */

    // Clean up
    delete canvasChi2;
    delete canvasBeta;
    delete canvasBetaPionCut;
    delete canvasBetaSigmaBin3;
    file->Close();
    delete file;

    return 0;
}