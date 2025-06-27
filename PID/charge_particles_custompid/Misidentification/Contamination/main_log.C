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
    TFile* file = new TFile("../pkptreeCxC_4.root");
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
    double pMin = 1.0, pMax = 4.0;
    int nBins = 10;
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
            chi2FitRanges[i] = {-2.0, 2.0};
            chi2FitParams[i] = {1000.0, 0.0, 1.0};
        } else if (pLow < 3.0) {
            chi2HistRanges[i] = {-4.0, 4.0};
            chi2FitRanges[i] = {-1.5, 1.5};
            chi2FitParams[i] = {800.0, 0.1, 0.8};
        } else {
            chi2HistRanges[i] = {-3.0, 3.0};
            chi2FitRanges[i] = {-1.2, 1.2};
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
    canvasChi2->Print("output/pdf/chi2pid_fits.pdf[");
    canvasBeta->Print("output/pdf/beta_before_cut.pdf[");
    canvasBeta->Print("output/pdf/beta_after_cut.pdf[");
    canvasBetaPionCut->Print("output/pdf/beta_pion_cut_only.pdf[");
    std::ofstream csvFile("output/csv/contamination_1to4GeV.csv");
    std::ofstream csvFileBefore("output/csv/contamination_before_chi2pid_1to4GeV.csv");
    std::ofstream csvFilePionCut("output/csv/contamination_pion_cut_only_1to4GeV.csv");
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
        TLegend* legPions = new TLegend(0.65, 0.6, 0.9, 0.9);
        legPions->AddEntry(chi2Pions, "Pions", "l");
        legPions->Draw();
        canvasChi2->cd(2);
        chi2Kaons->Draw("HIST");
        if (chi2Kaons->GetFunction("gausFit")) {
            chi2Kaons->GetFunction("gausFit")->SetLineColor(kGreen);
            chi2Kaons->GetFunction("gausFit")->Draw("SAME");
        }
        TLegend* legKaons = new TLegend(0.65, 0.6, 0.9, 0.9);
        legKaons->AddEntry(chi2Kaons, "Kaons", "l");
        legKaons->Draw();
        canvasChi2->Print("output/pdf/chi2pid_fits.pdf");
        delete legPions;
        delete legKaons;
        delete chi2Pions;
        delete chi2Kaons;

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
        }
        if (betaKaonsBefore->GetEntries() >= 10) {
            //TF1* gausKaonBefore = new TF1(TString::Format("gaus_%s", betaKaonsBefore->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            TF1* gausKaonBefore = new TF1(TString::Format("gaus_%s", betaKaonsBefore->GetName()), "gaus", 0.95, 1.03);

            if (i == 9) { // Last bin [3.7-4.0) GeV/c
                //gausKaonBefore->SetRange(0.9845, 0.9945);
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R","",0.984,0.994);
                
            } 
             else if (i == 8) { // Last bin [3.7-4.0) GeV/c
                //gausKaonBefore->SetRange(0.9845, 0.9945);
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.994, 0.1);
                betaKaonsBefore->Fit(gausKaonBefore, "R","",0.984,0.9945);
                
            } 
            else if ( i==5 ){
               // gausKaonBefore->SetRange(0.974, 0.992);
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.974, 0.992);
                betaKaonsBefore->Fit(gausKaonBefore, "R", "", 0.974, 0.992);
            }
            
            else {
                gausKaonBefore->SetParameters(betaKaonsBefore->GetMaximum(), 0.99, 0.005);
                betaKaonsBefore->Fit(gausKaonBefore, "R","",0.95,1.0);
            }
            
            kMeanBefore = gausKaonBefore->GetParameter(1);
            kSigmaBefore = gausKaonBefore->GetParameter(2);
            kConstBefore = gausKaonBefore->GetParameter(0);
        }

        // Plot beta histograms (before cut)
        canvasBeta->Clear();
         bool logScale = false;
       // bool logScale = true;
       // canvasBeta->SetLogy();
        double yMin = logScale ? 0.1 : 0.0;

        double yMax = -1.0;
        if (betaPionsBefore->GetEntries() >= 10) {
            betaPionsBefore->SetMarkerStyle(20);
            betaPionsBefore->SetMarkerSize(1.5);
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
            betaKaonsBefore->SetMarkerSize(1.5);
            betaKaonsBefore->SetMarkerColor(kGreen);
            betaKaonsBefore->Draw(betaPionsBefore->GetEntries() >= 10 ? "P SAME" : "P");
            yMax = std::max(yMax, betaKaonsBefore->GetMaximum());
            if (auto* fit = betaKaonsBefore->GetFunction(TString::Format("gaus_%s", betaKaonsBefore->GetName()))) {
                fit->SetLineColor(kGreen);
                fit->SetLineStyle(kDashed);
                fit->Draw("SAME");
                fit->SetRange(0.95,1.03);
            }
        }
        yMax *= logScale ? 10 : 1.2;
        if (betaPionsBefore->GetEntries() >= 10 || betaKaonsBefore->GetEntries() >= 10) {
            (betaPionsBefore->GetEntries() >= 10 ? betaPionsBefore : betaKaonsBefore)->GetYaxis()->SetRangeUser(yMin, yMax);
        }
        TLegend* leg = new TLegend(0.65, 0.6, 0.9, 0.9);
        if (betaPionsBefore->GetEntries() >= 10) leg->AddEntry(betaPionsBefore, "Pions", "p");
        if (betaKaonsBefore->GetEntries() >= 10) leg->AddEntry(betaKaonsBefore, "Kaons", "p");
        leg->Draw();
        canvasBeta->Print("output/pdf/beta_before_cut.pdf");
        delete leg;

        // Calculate contamination (before chi2pid cuts)
        double c1Before, c2Before;
        std::map<std::string, double> pFitBefore = {{"gaus_mean", pMeanBefore}, {"gaus_sigma", pSigmaBefore}, {"gaus_constant", pConstBefore}};
        std::map<std::string, double> kFitBefore = {{"gaus_mean", kMeanBefore}, {"gaus_sigma", kSigmaBefore}, {"gaus_constant", kConstBefore}};
        double pionLeft, pionRight, kaonRight;
        if (i == 0) { // Bin [1.0-1.3)
            pionLeft = 0.9880;
            pionRight = 1.0120;
            kaonRight = 0.0;
        } else if (i == 1) { // Bin [1.3-1.6)
            pionLeft = 0.9755;
            pionRight = 1.0135;
            kaonRight = 0.9721;
        } else if (i == 2) { // Bin [1.6-1.9)
            pionLeft = 0.98;
            pionRight = 1.0140;
            kaonRight = 0.986;
        } else if (i == 3) { // Bin [1.9-2.2)
            pionLeft = 0.9815;
            pionRight = 1.0145;
            kaonRight = 0.992;
        } else if (i == 4) { // Bin [2.2-2.5)
            pionLeft = 0.9880;
            pionRight = 1.0141;
            kaonRight = 0.9971;
        } else if (i == 5) { // Bin [2.5-2.8)
            pionLeft = 0.9848;
            pionRight = 1.013;
            kaonRight = 1.0;
        } else if (i == 6) { // Bin [2.8-3.1)
            pionLeft = 0.986;
            pionRight = 1.0128;
            kaonRight = 1.003;
        } else if (i == 7) { // Bin [3.1-3.4)
            pionLeft = 0.9868;
            pionRight = 1.0121;
            kaonRight = 1.0048;
        } else if (i == 8) { // Bin [3.4-3.7)
            pionLeft = 0.9879;
            pionRight = 1.0120;
            kaonRight = 1.0049;
        } else if (i == 9) { // Bin [3.7-4.0)
            pionLeft = 0.9880;
            pionRight = 1.0120;
            kaonRight = 1.012;
        } else {
            pionLeft = 0.9880;
            pionRight = 1.0120;
            kaonRight = 1.0;
        }
        double contaminationBefore = calculateOverlapContamination(pFitBefore, kFitBefore, pionLeft, pionRight, kaonRight, c1Before, c2Before);
        double kaonLeftBefore = (c1Before >= kMeanBefore && c1Before <= pMeanBefore) ? c1Before : (c2Before >= kMeanBefore && c2Before <= pMeanBefore) ? c2Before : -1.0;

        // Save contamination (before chi2pid cuts)
        if (contaminationBefore < 0 || kaonLeftBefore == -1.0) {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Contamination (before chi2pid) calculation failed\n";
            csvFileBefore << pLow << "-" << pHigh << ",N/A,N/A\n";
        } else {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Before chi2pid - Pion range: [" << pionLeft << ", " << pionRight 
                      << "], Kaon range: [" << kaonLeftBefore << ", " << kaonRight << "], Contamination: " 
                      << contaminationBefore << "%\n";
            csvFileBefore << pLow << "-" << pHigh << "," << kaonLeftBefore << "," << contaminationBefore << "\n";
        }
        delete betaPionsBefore;
        delete betaKaonsBefore;

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
        } else {
            // Fit beta histograms (pion chi2pid cut only)
            double pMeanPionCut = 0, pSigmaPionCut = 0, pConstPionCut = 0, kMeanPionCut = 0, kSigmaPionCut = 0, kConstPionCut = 0;
            if (betaPionsAfterPionCut->GetEntries() >= 10) {
                TF1* gausPionPionCut = new TF1(TString::Format("gaus_%s", betaPionsAfterPionCut->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
                gausPionPionCut->SetParameters(betaPionsAfterPionCut->GetMaximum(), 1.0, 0.005);
                betaPionsAfterPionCut->Fit(gausPionPionCut, "R", "", betaFitRanges[i].first, betaFitRanges[i].second);
                pMeanPionCut = gausPionPionCut->GetParameter(1);
                pSigmaPionCut = gausPionPionCut->GetParameter(2);
                pConstPionCut = gausPionPionCut->GetParameter(0);
            }
            if (betaKaonsAfterPionCut->GetEntries() >= 10) {
                TF1* gausKaonPionCut = new TF1(TString::Format("gaus_%s", betaKaonsAfterPionCut->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
                if (i == 9) { // Last bin [3.7-4.0) GeV/c
                    gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.1);
                    betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.984, 0.9944);
                } else if ( i == 8 ){
                    gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.994, 0.1);
                    betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R", "", 0.9845, 0.994);
                }
                
                else {
                    gausKaonPionCut->SetParameters(betaKaonsAfterPionCut->GetMaximum(), 0.99, 0.005);
                    betaKaonsAfterPionCut->Fit(gausKaonPionCut, "R");
                }
                
                kMeanPionCut = gausKaonPionCut->GetParameter(1);
                kSigmaPionCut = gausKaonPionCut->GetParameter(2);
                kConstPionCut = gausKaonPionCut->GetParameter(0);
            }

            // Calculate contamination (pion chi2pid cut only)
            double c1PionCut, c2PionCut;
            std::map<std::string, double> pFitPionCut = {{"gaus_mean", pMeanPionCut}, {"gaus_sigma", pSigmaPionCut}, {"gaus_constant", pConstPionCut}};
            std::map<std::string, double> kFitPionCut = {{"gaus_mean", kMeanPionCut}, {"gaus_sigma", kSigmaPionCut}, {"gaus_constant", kConstPionCut}};
            double contaminationPionCut = calculateOverlapContamination(pFitPionCut, kFitPionCut, pionLeft, pionRight, kaonRight, c1PionCut, c2PionCut);
            double kaonLeftPionCut = (c1PionCut >= kMeanPionCut && c1PionCut <= pMeanPionCut) ? c1PionCut : (c2PionCut >= kMeanPionCut && c2PionCut <= pMeanPionCut) ? c2PionCut : -1.0;

            // Plot beta histograms (pion chi2pid cut only)
            canvasBetaPionCut->Clear();
            logScale = true; // Toggle log scale here
            canvasBetaPionCut->SetLogy(); // Uncomment to enable log scale
            yMin = logScale ? 0.1 : 0.0;
            yMax = -1.0;
            if (betaPionsAfterPionCut->GetEntries() >= 10) {
                betaPionsAfterPionCut->SetMarkerStyle(20);
                betaPionsAfterPionCut->SetMarkerSize(1.5);
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
                betaKaonsAfterPionCut->SetMarkerSize(1.5);
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
            yMax *= logScale ? 10 : 1.2;
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
            leg = new TLegend(0.65, 0.6, 0.9, 0.9);
            if (betaPionsAfterPionCut->GetEntries() >= 10) leg->AddEntry(betaPionsAfterPionCut, "Pions", "p");
            if (betaKaonsAfterPionCut->GetEntries() >= 10) leg->AddEntry(betaKaonsAfterPionCut, "Kaons", "p");
            leg->Draw();
            canvasBetaPionCut->Print("output/pdf/beta_pion_cut_only.pdf");
            delete leg;

            // Save contamination (pion chi2pid cut only)
            if (contaminationPionCut < 0 || kaonLeftPionCut == -1.0) {
                std::cout << "Bin [" << pLow << "-" << pHigh << "): Contamination (pion cut only) calculation failed\n";
                csvFilePionCut << pLow << "-" << pHigh << ",N/A,N/A\n";
            } else {
                std::cout << "Bin [" << pLow << "-" << pHigh << "): Pion cut only - Pion range: [" << pionLeft << ", " << pionRight 
                          << "], Kaon range: [" << kaonLeftPionCut << ", " << kaonRight << "], Contamination: " 
                          << contaminationPionCut << "%\n";
                csvFilePionCut << pLow << "-" << pHigh << "," << kaonLeftPionCut << "," << contaminationPionCut << "\n";
            }
        }
        delete betaPionsAfterPionCut;
        delete betaKaonsAfterPionCut;

        // Create beta histograms (after cut, both pions and kaons)
        TH1F* betaPionsAfter = new TH1F(TString::Format("beta_pions_after_%d", i),
                                        TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                        70, betaHistRanges[i].first, betaHistRanges[i].second);
        TH1F* betaKaonsAfter = new TH1F(TString::Format("beta_kaons_after_%d", i),
                                        TString::Format("p: [%.2f-%.2f) GeV/c;beta;Counts", pLow, pHigh),
                                        70, betaHistRanges[i].first, betaHistRanges[i].second);

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
        }
        if (betaKaonsAfter->GetEntries() >= 10) {
            TF1* gausKaon = new TF1(TString::Format("gaus_%s", betaKaonsAfter->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            if (i == 9) { // Last bin [3.7-4.0) GeV/c
                gausKaon->SetRange(0.97, 1.01);
                gausKaon->SetParameters(betaKaonsAfter->GetMaximum(), 0.994, 0.1);
            } else {
                gausKaon->SetParameters(betaKaonsAfter->GetMaximum(), 0.99, 0.005);
            }
            betaKaonsAfter->Fit(gausKaon, "R");
            kMean = gausKaon->GetParameter(1);
            kSigma = gausKaon->GetParameter(2);
            kConst = gausKaon->GetParameter(0);
        }

        // Calculate contamination (after chi2pid cuts, both pions and kaons)
        double c1, c2;
        std::map<std::string, double> pFit = {{"gaus_mean", pMean}, {"gaus_sigma", pSigma}, {"gaus_constant", pConst}};
        std::map<std::string, double> kFit = {{"gaus_mean", kMean}, {"gaus_sigma", kSigma}, {"gaus_constant", kConst}};
        double contamination = calculateOverlapContamination(pFit, kFit, pionLeft, pionRight, kaonRight, c1, c2);
        double kaonLeft = (c1 >= kMean && c1 <= pMean) ? c1 : (c2 >= kMean && c2 <= pMean) ? c2 : -1.0;

        // Plot beta histograms (after cut, both pions and kaons)
        canvasBeta->Clear();
        logScale = true; // Toggle log scale here
        canvasBeta->SetLogy(); // Uncomment to enable log scale
        yMin = logScale ? 0.1 : 0.0;
        yMax = -1.0;
        if (betaPionsAfter->GetEntries() >= 10) {
            betaPionsAfter->SetMarkerStyle(20);
            betaPionsAfter->SetMarkerSize(1.5);
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
            betaKaonsAfter->SetMarkerSize(1.5);
            betaKaonsAfter->SetMarkerColor(kGreen);
            betaKaonsAfter->Draw(betaPionsAfter->GetEntries() >= 10 ? "P SAME" : "P");
            yMax = std::max(yMax, betaKaonsAfter->GetMaximum());
            if (auto* fit = betaKaonsAfter->GetFunction(TString::Format("gaus_%s", betaKaonsAfter->GetName()))) {
                fit->SetLineColor(kGreen);
                fit->SetLineStyle(kDashed);
                fit->Draw("SAME");
            }
        }
        yMax *= logScale ? 10 : 1.2;
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
        leg = new TLegend(0.65, 0.6, 0.9, 0.9);
        if (betaPionsAfter->GetEntries() >= 10) leg->AddEntry(betaPionsAfter, "Pions", "p");
        if (betaKaonsAfter->GetEntries() >= 10) leg->AddEntry(betaKaonsAfter, "Kaons", "p");
        leg->Draw();
        canvasBeta->Print("output/pdf/beta_after_cut.pdf");
        delete leg;

        // Save contamination (after chi2pid cuts, both pions and kaons)
        if (contamination < 0 || kaonLeft == -1.0) {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Contamination (both cuts) calculation failed\n";
            csvFile << pLow << "-" << pHigh << ",N/A,N/A\n";
        } else {
            std::cout << "Bin [" << pLow << "-" << pHigh << "): Both cuts - Pion range: [" << pionLeft << ", " << pionRight 
                      << "], Kaon range: [" << kaonLeft << ", " << kaonRight << "], Contamination: " 
                      << contamination << "%\n";
            csvFile << pLow << "-" << pHigh << "," << kaonLeft << "," << contamination << "\n";
        }
        delete betaPionsAfter;
        delete betaKaonsAfter;
    }

    // Close PDFs and CSVs
    canvasChi2->Print("output/pdf/chi2pid_fits.pdf]");
    canvasBeta->Print("output/pdf/beta_before_cut.pdf]");
    canvasBeta->Print("output/pdf/beta_after_cut.pdf]");
    canvasBetaPionCut->Print("output/pdf/beta_pion_cut_only.pdf]");
    csvFile.close();
    csvFileBefore.close();
    csvFilePionCut.close();

    // Clean up
    delete canvasChi2;
    delete canvasBeta;
    delete canvasBetaPionCut;
    file->Close();
    delete file;

    return 0;
}