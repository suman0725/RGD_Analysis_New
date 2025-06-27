#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TF1.h>
#include <TH2F.h>
#include <TLine.h>
#include <TSystem.h>
#include <fstream>
#include <iostream>

void chi2pid_310_345() {
    // Enable batch mode
    gROOT->SetBatch(kTRUE);

    // Open the ROOT file
    TFile* file = new TFile("pkptreeCxC_4.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open pkptreeCxC_4.root" << std::endl;
        return;
    }

    // Read the trees
    TTree* treeEBPions = (TTree*)file->Get("EB_pid_pions");
    TTree* treeEBKaons = (TTree*)file->Get("EB_pid_kaons");
    TTree* treePionHypothesis = (TTree*)file->Get("pion_hypothesis");
    TTree* treeKaonHypothesis = (TTree*)file->Get("kaon_hypothesis");
    if (!treeEBPions || !treePionHypothesis || !treeEBKaons) {
        std::cerr << "Error: Cannot find one or more trees in pkptreeCxC_4.root" << std::endl;
        file->Close();
        return;
    }

    // Define the momentum bin
    const float pMin = 2.50;
    const float pMax = 2.80;

    // Define chi2pid histogram range
    const int nChi2Bins = 100;
    const float chi2Min = -15.0;
    const float chi2Max = 15.0;

    // Create histograms
    TH1F* histEBPions = new TH1F("h_chi2pid_ebpions", "Chi2pid for EB Pions (2.50 < p < 2.80 GeV)", nChi2Bins, chi2Min, chi2Max);
    TH1F* histPionHypothesis = new TH1F("h_chi2pid_pionhypothesis", "Chi2pid for Pion Hypothesis (2.50 < p < 2.80 GeV)", nChi2Bins, chi2Min, chi2Max);
    TH1F* histKaonHypothesis = new TH1F("h_chi2pid_kaonhypothesis", "Chi2pid for Kaon Hypothesis (2.50 < p < 2.80 GeV)", nChi2Bins, chi2Min, chi2Max);
    TH1F* histEBKaons = new TH1F("h_chi2pid_ebkaons", "Chi2pid for EB Kaons (2.50 < p < 2.80 GeV)", nChi2Bins, chi2Min, chi2Max);
    TH1F* histEBPions_beta_vs_p_bc = new TH1F("h_beta_vs_p_ebpions_bc", "beta for EB PID +ve Pions Vs momentum", nChi2Bins, 0.96, 1.03);
    TH1F* histEBKaons_beta_vs_p_bc = new TH1F("h_beta_vs_p_ebkaons_bc", "beta for EB PID +ve Kaon Vs momentum", nChi2Bins, 0.96, 1);
    TH1F* histEBPions_beta_vs_p_ac = new TH1F("h_beta_vs_p_ebpions_ac", "beta for EB PID +ve Pions Vs momentum", nChi2Bins, 0.96, 1.03);
    TH1F* histEBKaons_beta_vs_p_ac = new TH1F("h_beta_vs_p_ebkaons_ac", "beta for EB PID +ve Kaon Vs momentum", nChi2Bins, 0.96, 1);
    TH2F* histEBPions_mass_vs_p = new TH2F("h_calmass_vs_p_ebpions", "FTOF mass for EB PID +ve Pions Vs momentum;Momentum (GeV); Mass (GeV)", 200, 0, 10, 200, 0, 1);
    TH2F* histEBPions_mass_vs_p_zoom = new TH2F("h_calmass_vs_p_ebpions_zoom", "FTOF mass for EB PID +ve Pions Vs momentum;Momentum (GeV); Mass (GeV)", 200, 0, 8, 200, 0, 0.35);
    TH2F* histPHPions_mass_vs_p = new TH2F("h_calmass_vs_p_phpions", "FTOF mass for +ve charge particles (except positron) Vs momentum;Momentum (GeV); Mass (GeV)", 200, 0, 10, 200, 0, 1.5);
    TH2F* histPHPions_mass_vs_p_zoom = new TH2F("h_calmass_vs_p_phpions_zoom", "FTOF mass for +ve pions (hard mass cut at 0.35 GeV) (except positron) Vs momentum;Momentum (GeV); Mass (GeV)", 200, 0, 8, 200, 0, 0.35);

    // Prevent histograms from being owned by TFile
    histEBPions->SetDirectory(0);
    histPionHypothesis->SetDirectory(0);
    histEBKaons->SetDirectory(0);
    histEBPions_beta_vs_p_bc->SetDirectory(0);
    histEBKaons_beta_vs_p_bc->SetDirectory(0);
    histEBPions_beta_vs_p_ac->SetDirectory(0);
    histEBKaons_beta_vs_p_ac->SetDirectory(0);
    histEBPions_mass_vs_p->SetDirectory(0);
    histEBPions_mass_vs_p_zoom->SetDirectory(0);
    histPHPions_mass_vs_p->SetDirectory(0);
    histPHPions_mass_vs_p_zoom->SetDirectory(0);

    // Set up the trees to read chi2pid and p
    float chi2pidEBPions, pEBPions;
    float betaEBPions, betaEBKaons;
    float EBPionsmass, PHPionsmass, chi2pidEBKaons, pEBKaons;
    float chi2pidPionHypothesis, pPionHypothesis, chi2pidKaonHypothesis;
    int pidPionHypothesis;
    treeEBPions->SetBranchAddress("chi2pid", &chi2pidEBPions);
    treeEBPions->SetBranchAddress("p", &pEBPions);
    treeEBPions->SetBranchAddress("beta", &betaEBPions);
    treeEBPions->SetBranchAddress("mass", &EBPionsmass);
    treePionHypothesis->SetBranchAddress("chi2pid", &chi2pidPionHypothesis);
    treePionHypothesis->SetBranchAddress("p", &pPionHypothesis);
    treePionHypothesis->SetBranchAddress("pid", &pidPionHypothesis);
    treePionHypothesis->SetBranchAddress("mass", &PHPionsmass);
    treeEBKaons->SetBranchAddress("chi2pid", &chi2pidEBKaons);
    treeEBKaons->SetBranchAddress("p", &pEBKaons);
    treeEBKaons->SetBranchAddress("beta", &betaEBKaons);
    treeKaonHypothesis->SetBranchAddress("chi2pid_kaon", &chi2pidKaonHypothesis);

    // Fill histograms for EB pions
    Long64_t nEntriesEBPions = treeEBPions->GetEntries();
    for (Long64_t j = 0; j < nEntriesEBPions; j++) {
        treeEBPions->GetEntry(j);
        if (pEBPions >= pMin && pEBPions < pMax && chi2pidEBPions != 99999.0f) {
            histEBPions->Fill(chi2pidEBPions);
            histEBPions_beta_vs_p_bc->Fill(betaEBPions);
            if (chi2pidEBPions > -2.97 && chi2pidEBPions < 2.84) {
                histEBPions_beta_vs_p_ac->Fill(betaEBPions);
            }
        }
        histEBPions_mass_vs_p->Fill(pEBPions, EBPionsmass);
        if (EBPionsmass < 0.35) {
            histEBPions_mass_vs_p_zoom->Fill(pEBPions, EBPionsmass);
        }
    }

    // Fill histograms for EB kaons
    Long64_t nEntriesEBKaons = treeEBKaons->GetEntries();
    for (Long64_t j = 0; j < nEntriesEBKaons; j++) {
        treeEBKaons->GetEntry(j);
        if (pEBKaons >= pMin && pEBKaons < pMax && chi2pidEBKaons != 99999.0f) {
            histEBKaons->Fill(chi2pidEBKaons);
            histEBKaons_beta_vs_p_bc->Fill(betaEBKaons);
            if (chi2pidEBKaons > -3.49 && chi2pidEBKaons < 1.10) {
                histEBKaons_beta_vs_p_ac->Fill(betaEBKaons);
            }
        }
    }


    // Fill histograms for manual pion hypothesis
    Long64_t nEntriesPionHypothesis = treePionHypothesis->GetEntries();
    int nKaons = 0, nProtons = 0, nPions = 0;
    for (Long64_t j = 0; j < nEntriesPionHypothesis; j++) {
        treePionHypothesis->GetEntry(j);
        if (pPionHypothesis >= pMin && pPionHypothesis < pMax && chi2pidPionHypothesis != 99999.0f) {
            histPionHypothesis->Fill(chi2pidPionHypothesis);
            if (pidPionHypothesis == 321) nKaons++;
            if (pidPionHypothesis == 2212) nProtons++;
            if (pidPionHypothesis == 211) nPions++;
        }
        histPHPions_mass_vs_p->Fill(pPionHypothesis, PHPionsmass);
        if (PHPionsmass < 0.35) {
            histPHPions_mass_vs_p_zoom->Fill(pPionHypothesis, PHPionsmass);
        }
    }
    // Fill histograms for manual kaon hypothesis   
    Long64_t nEntriesKaonHypothesis = treeKaonHypothesis->GetEntries();
    for (Long64_t j = 0; j < nEntriesKaonHypothesis; j++) {
        treeKaonHypothesis->GetEntry(j);
        if (pPionHypothesis >= pMin && pPionHypothesis < pMax && chi2pidKaonHypothesis != 99999.0f) {
            histKaonHypothesis->Fill(chi2pidKaonHypothesis);
        }
    }
    // Fit the EB pions histogram with a Gaussian
    TF1* gausFit = new TF1("gausFit", "gaus", -5, 5);
    gausFit->SetParameters(6200, 0, 1);
    gausFit->SetLineStyle(2);
    gausFit->SetLineColor(kBlack);
    histEBPions->Fit(gausFit, "R");
    double mean = gausFit->GetParameter(1);
    double sigma = gausFit->GetParameter(2);
    double gausChi2 = gausFit->GetChisquare();
    double gausNDF = gausFit->GetNDF();

    // Fit with Crystal Ball
    TF1* cb = new TF1("cb", "crystalball", -5, 5);
    cb->SetParameters(gausFit->GetParameter(0), gausFit->GetParameter(1), gausFit->GetParameter(2), 1.5, 2.0);
    cb->SetLineColor(kMagenta + 2);
    cb->SetLineWidth(2);
    histEBPions->Fit(cb, "R+");
    double cbMean = cb->GetParameter(1);
    double cbSigma = cb->GetParameter(2);
    double cbChi2 = cb->GetChisquare();
    double cbNDF = cb->GetNDF();

    // Define a Double Gaussian
    TF1* doubleGausFit = new TF1("doubleGausFit", "gaus(0) + gaus(3)", -5, 9);
    doubleGausFit->SetParameters(6000, 0, 1, 900, 5, 1);
    doubleGausFit->SetLineColor(kGreen);
    doubleGausFit->SetLineWidth(1);

    // Define a Double Crystal Ball
    TF1* doubleCBFit = new TF1("doubleCBFit", "crystalball(0) + crystalball(5)", -5, 9);
    doubleCBFit->SetParameters(6000, 0, 1, 1.5, 2.0, 900, 5, 1, 1.5, 2.0);
    doubleCBFit->SetParNames("Const1", "Mean1", "Sigma1", "Alpha1", "n1", "Const2", "Mean2", "Sigma2", "Alpha2", "n2");
    doubleCBFit->SetLineColor(kBlack);
    doubleCBFit->SetLineWidth(1);

    // Calculate particles in 3σ range
    int nPionsInRange = 0, nKaonsInRange = 0, nProtonsInRange = 0;
    for (Long64_t j = 0; j < nEntriesPionHypothesis; j++) {
        treePionHypothesis->GetEntry(j);
        if (pPionHypothesis >= pMin && pPionHypothesis < pMax &&
            chi2pidPionHypothesis >= (mean - 3 * sigma) &&
            chi2pidPionHypothesis <= (mean + 3 * sigma)) {
            if (pidPionHypothesis == 211) nPionsInRange++;
            if (pidPionHypothesis == 321) nKaonsInRange++;
            if (pidPionHypothesis == 2212) nProtonsInRange++;
        }
    }

    // Calculate efficiency and contamination
    int nEBPionsInRange = histEBPions->Integral(histEBPions->FindBin(mean - 3 * sigma),
                                                histEBPions->FindBin(mean + 3 * sigma));
    double efficiency = nPions ? (double)nPionsInRange / nPions : 0;
    double kaon_contamination = nEBPionsInRange ? (double)nKaonsInRange / nEBPionsInRange : 0;
    double proton_contamination = nEBPionsInRange ? (double)nProtonsInRange / nEBPionsInRange : 0;

    // Create output directories
    TString parentDir = "output_1";
    gSystem->Exec(Form("rm -rf %s", parentDir.Data())); // Remove existing parent directory
    if (gSystem->mkdir(parentDir, kTRUE) != 0) {
        std::cerr << "Error: Failed to create parent directory " << parentDir << std::endl;
        return;
    }
    TString outputDir = Form("%s/%.2f_to_%.2f_GeV", parentDir.Data(), pMin, pMax);
    if (gSystem->mkdir(outputDir, kTRUE) != 0) {
        std::cerr << "Error: Failed to create directory " << outputDir << std::endl;
        return;
    }
    TString outputDir_1 = "output_1/mass_vs_p";
    if (gSystem->mkdir(outputDir_1, kTRUE) != 0) {
        std::cerr << "Error: Failed to create directory " << outputDir_1 << std::endl;
        return;
    }
    std::cout << "Created directories: " << parentDir << ", " << outputDir << ", " << outputDir_1 << std::endl;

    // Verify histogram content
    std::cout << "histEBPions entries: " << histEBPions->GetEntries() << std::endl;
    std::cout << "histPionHypothesis entries: " << histPionHypothesis->GetEntries() << std::endl;
    std::cout << "histEBKaons entries: " << histEBKaons->GetEntries() << std::endl;

    // Save results to CSV
    TString csvFileName = outputDir + "/efficiency_contamination.csv";
    gSystem->Unlink(csvFileName);
    std::ofstream csvFile(csvFileName.Data());
    csvFile << "Momentum Bin,Mean,Sigma,3Sigma Range,Total Pions,Total Pions in 3Sigma,Kaons in 3Sigma,Protons in 3Sigma,Efficiency,Kaon Contamination,Proton Contamination\n";
    csvFile << pMin << "-" << pMax << "," << mean << "," << sigma << ","
            << "[" << mean - 3 * sigma << "," << mean + 3 * sigma << "],"
            << nPions << "," << nPionsInRange << "," << nKaonsInRange << "," << nProtonsInRange << ","
            << efficiency << "," << kaon_contamination << "," << proton_contamination << "\n";
    csvFile.close();

    // First canvas for pions
    TCanvas* canvas = new TCanvas("canvas", "Chi2pid Overlay", 800, 600);
    histEBPions->SetLineColor(kBlue);
    histPionHypothesis->SetLineColor(kRed);
    histEBPions->Draw();
    histPionHypothesis->Draw("SAME");
    gausFit->Draw("SAME");
    cb->Draw("SAME");

    // Add vertical lines for mean ± 3σ
    TLine* lineMean = new TLine(mean, 0, mean, histEBPions->GetMaximum());
    TLine* lineMeanPlus3Sigma = new TLine(mean + 3 * sigma, 0, mean + 3 * sigma, histEBPions->GetMaximum());
    TLine* lineMeanMinus3Sigma = new TLine(mean - 3 * sigma, 0, mean - 3 * sigma, histEBPions->GetMaximum());
    lineMean->SetLineColor(kBlack);
    lineMeanPlus3Sigma->SetLineColor(kBlack);
    lineMeanMinus3Sigma->SetLineColor(kBlack);
    lineMean->Draw("SAME");
    lineMeanPlus3Sigma->Draw("SAME");
    lineMeanPlus3Sigma->SetLineStyle(2);
    lineMeanMinus3Sigma->Draw("SAME");
    lineMeanMinus3Sigma->SetLineStyle(2);

    TLine* lineCBMean = new TLine(cbMean, 0, cbMean, histEBPions->GetMaximum());
    TLine* lineCBMeanPlus3Sigma = new TLine(cbMean + 3 * cbSigma, 0, cbMean + 3 * cbSigma, histEBPions->GetMaximum());
    TLine* lineCBMeanMinus3Sigma = new TLine(cbMean - 3 * cbSigma, 0, cbMean - 3 * cbSigma, histEBPions->GetMaximum());
    lineCBMean->SetLineColor(kMagenta + 2);
    lineCBMeanPlus3Sigma->SetLineColor(kMagenta + 2);
    lineCBMeanMinus3Sigma->SetLineColor(kMagenta + 2);
    lineCBMean->Draw("SAME");
    lineCBMeanPlus3Sigma->Draw("SAME");
    lineCBMeanPlus3Sigma->SetLineStyle(2);
    lineCBMeanMinus3Sigma->Draw("SAME");
    lineCBMeanMinus3Sigma->SetLineStyle(2);

    // Add legend
    TLegend* legend = new TLegend(0.7, 0.5, 0.9, 0.9);
    legend->SetTextSize(0.02);
    legend->AddEntry(histEBPions, "EB Pions (pid == 211)", "l");
    legend->AddEntry(histPionHypothesis, "Manual Pion Hypothesis", "l");
    legend->AddEntry(gausFit, "Gaussian Fit (EB pions: Initial)", "l");
    legend->AddEntry((TObject*)0, Form("Constant: %.2f", gausFit->GetParameter(0)), "");
    legend->AddEntry((TObject*)0, Form("Mean: %.2f", gausFit->GetParameter(1)), "");
    legend->AddEntry((TObject*)0, Form("Sigma: %.2f", gausFit->GetParameter(2)), "");
    legend->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.1f / %d = %.2f", gausChi2, (int)gausNDF, gausChi2 / gausNDF), "");
    legend->AddEntry((TObject*)0, Form("Mean #pm 3#sigma: [%.2f, %.2f]", mean - 3 * sigma, mean + 3 * sigma), "");
    legend->AddEntry(lineMean, "Mean", "l");
    legend->AddEntry(cb, "Crystal Ball function (EB pions):Final", "l");
    legend->AddEntry((TObject*)0, Form("Constant: %.2f", cb->GetParameter(0)), "");
    legend->AddEntry((TObject*)0, Form("Mean: %.4f", cb->GetParameter(1)), "");
    legend->AddEntry((TObject*)0, Form("Sigma: %.4f", cb->GetParameter(2)), "");
    legend->AddEntry((TObject*)0, Form("Alpha: %.4f", cb->GetParameter(3)), "");
    legend->AddEntry((TObject*)0, Form("n: %.4f", cb->GetParameter(4)), "");
    legend->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.1f / %d = %.2f", cbChi2, (int)cbNDF, cbChi2 / cbNDF), "");
    legend->AddEntry((TObject*)0, Form("Mean #pm 3#sigma: [%.2f, %.2f]", cbMean - 3 * cbSigma, cbMean + 3 * cbSigma), "");
    legend->Draw();

    // Save chi2pid_overlay.pdf
    TString pdfFileName = outputDir + "/chi2pid_overlay.pdf";
    gSystem->Unlink(pdfFileName);
    canvas->Update();
    canvas->Print(pdfFileName + "[");
    canvas->Print(pdfFileName);

    // Second canvas for kaons
    // Second canvas for kaons
TCanvas* canvas2 = new TCanvas("canvas2", "Kaons", 800, 600);
if (/* histEBKaons && */ histEBKaons->GetEntries() > 0/*  && histKaonHypothesis *//*  && histKaonHypothesis->GetEntries() > 0 */) {
    histEBKaons->SetLineColor(kGreen);
    histEBKaons->Draw();
    histKaonHypothesis->SetLineColor(kRed);
    std::cout << "histKaonHypothesis: " << (histKaonHypothesis ? "Defined" : "Not Defined") << ", Entries: " << (histKaonHypothesis ? histKaonHypothesis->GetEntries() : 0) << std::endl;
    histKaonHypothesis->Draw("SAME");


    // Fit with Gaussian
    TF1* gausFitKaon = new TF1("gausFitKaon", "gaus", -5, 5);
    gausFitKaon->SetParameters(1000, 0, 1); // Initial guess: amplitude, mean, sigma
    gausFitKaon->SetLineStyle(2);
    gausFitKaon->SetLineColor(kBlack);
    histEBKaons->Fit(gausFitKaon, "R+");
    double kaonMean = gausFitKaon->GetParameter(1);
    double kaonSigma = gausFitKaon->GetParameter(2);
    double kaonGausChi2 = gausFitKaon->GetChisquare();
    double kaonGausNDF = gausFitKaon->GetNDF();

    // Fit with Crystal Ball, adjusted for right-side tail and extended range
    TF1* cbKaon = new TF1("cbKaon", "crystalball", -5, 15); // Extended range to capture the tail
    cbKaon->SetParameters(
        gausFitKaon->GetParameter(0), // Initial constant
        gausFitKaon->GetParameter(1), // Initial mean
        gausFitKaon->GetParameter(2), // Initial sigma
        -1.0,                         // Alpha < 0 for right-side tail
        3.0                           // Adjusted n for tail shape
    );
    cbKaon->SetParLimits(0, 0, 2000);        // Constrain constant to positive, reasonable range
    cbKaon->SetParLimits(1, -2, 2);          // Constrain mean around peak
    cbKaon->SetParLimits(2, 0.5, 2.0);       // Constrain sigma to reasonable width
    cbKaon->SetParLimits(3, -5.0, 0.0);      // Constrain alpha to negative (right tail)
    cbKaon->SetParLimits(4, 1.0, 10.0);      // Constrain n to reasonable range
    cbKaon->SetLineColor(kMagenta + 2);
    cbKaon->SetLineWidth(2);
    histEBKaons->Fit(cbKaon, "R+");
    double kaonCBMean = cbKaon->GetParameter(1);
    double kaonCBSigma = cbKaon->GetParameter(2);
    double kaonCBChi2 = cbKaon->GetChisquare();
    double kaonCBNDF = cbKaon->GetNDF();

    // Draw mean ± 3σ lines for Gaussian
    TLine* lineKaonMean = new TLine(kaonMean, 0, kaonMean, histEBKaons->GetMaximum());
    TLine* lineKaonMeanPlus3Sigma = new TLine(kaonMean + 3 * kaonSigma, 0, kaonMean + 3 * kaonSigma, histEBKaons->GetMaximum());
    TLine* lineKaonMeanMinus3Sigma = new TLine(kaonMean - 3 * kaonSigma, 0, kaonMean - 3 * kaonSigma, histEBKaons->GetMaximum());
    lineKaonMean->SetLineColor(kBlack);
    lineKaonMeanPlus3Sigma->SetLineColor(kBlack);
    lineKaonMeanMinus3Sigma->SetLineColor(kBlack);
    lineKaonMean->Draw("SAME");
    lineKaonMeanPlus3Sigma->Draw("SAME");
    lineKaonMeanPlus3Sigma->SetLineStyle(2);
    lineKaonMeanMinus3Sigma->Draw("SAME");
    lineKaonMeanMinus3Sigma->SetLineStyle(2);

    // Draw mean ± 3σ lines for Crystal Ball
    TLine* lineKaonCBMean = new TLine(kaonCBMean, 0, kaonCBMean, histEBKaons->GetMaximum());
    TLine* lineKaonCBMeanPlus3Sigma = new TLine(kaonCBMean + 3 * kaonCBSigma, 0, kaonCBMean + 3 * kaonCBSigma, histEBKaons->GetMaximum());
    TLine* lineKaonCBMeanMinus3Sigma = new TLine(kaonCBMean - 3 * kaonCBSigma, 0, kaonCBMean - 3 * kaonCBSigma, histEBKaons->GetMaximum());
    lineKaonCBMean->SetLineColor(kMagenta + 2);
    lineKaonCBMeanPlus3Sigma->SetLineColor(kMagenta + 2);
    lineKaonCBMeanMinus3Sigma->SetLineColor(kMagenta + 2);
    lineKaonCBMean->Draw("SAME");
    lineKaonCBMeanPlus3Sigma->Draw("SAME");
    lineKaonCBMeanPlus3Sigma->SetLineStyle(2);
    lineKaonCBMeanMinus3Sigma->Draw("SAME");
    lineKaonCBMeanMinus3Sigma->SetLineStyle(2);

    // Redraw Gaussian fit on top to ensure visibility
    gausFitKaon->Draw("SAME");

    // Add legend
    TLegend* legendKaon = new TLegend(0.7, 0.5, 0.9, 0.9);
    legendKaon->SetTextSize(0.02);
    legendKaon->AddEntry(histEBKaons, "EB Kaons (pid == 321)", "l");
    legendKaon->AddEntry(gausFitKaon, "Gaussian Fit (EB Kaons: Initial)", "l");
    legendKaon->AddEntry((TObject*)0, Form("Constant: %.2f", gausFitKaon->GetParameter(0)), "");
    legendKaon->AddEntry((TObject*)0, Form("Mean: %.2f", gausFitKaon->GetParameter(1)), "");
    legendKaon->AddEntry((TObject*)0, Form("Sigma: %.2f", gausFitKaon->GetParameter(2)), "");
    legendKaon->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.1f / %d = %.2f", kaonGausChi2, (int)kaonGausNDF, kaonGausChi2 / kaonGausNDF), "");
    legendKaon->AddEntry((TObject*)0, Form("Mean #pm 3#sigma: [%.2f, %.2f]", kaonMean - 3 * kaonSigma, kaonMean + 3 * kaonSigma), "");
    legendKaon->AddEntry(lineKaonMean, "Mean (Gaussian)", "l");
    legendKaon->AddEntry(cbKaon, "Crystal Ball function (EB Kaons): Final", "l");
    legendKaon->AddEntry((TObject*)0, Form("Constant: %.2f", cbKaon->GetParameter(0)), "");
    legendKaon->AddEntry((TObject*)0, Form("Mean: %.4f", cbKaon->GetParameter(1)), "");
    legendKaon->AddEntry((TObject*)0, Form("Sigma: %.4f", cbKaon->GetParameter(2)), "");
    legendKaon->AddEntry((TObject*)0, Form("Alpha: %.4f", cbKaon->GetParameter(3)), "");
    legendKaon->AddEntry((TObject*)0, Form("n: %.4f", cbKaon->GetParameter(4)), "");
    legendKaon->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.1f / %d = %.2f", kaonCBChi2, (int)kaonCBNDF, kaonCBChi2 / kaonCBNDF), "");
    legendKaon->AddEntry((TObject*)0, Form("Mean #pm 3#sigma: [%.2f, %.2f]", kaonCBMean - 3 * kaonCBSigma, kaonCBMean + 3 * kaonCBSigma), "");
    legendKaon->Draw();
} else {
    std::cerr << "Error: histEBKaons is empty or not defined!" << std::endl;
    TH1F* placeholder = new TH1F("placeholder", "No Kaon Data", 10, -10, 10);
    placeholder->Draw();
    delete placeholder;
}
canvas2->Update();
canvas2->Print(pdfFileName);    // Write second page
canvas2->Print(pdfFileName + "]"); // Close PDF

    // Clear fit functions
    histEBPions->GetListOfFunctions()->Clear();
    histPionHypothesis->GetListOfFunctions()->Clear();
    histEBKaons->GetListOfFunctions()->Clear();

    // Canvas for double Gaussian and CB fits
    TCanvas* canvas_1 = new TCanvas("canvas_1", "Chi2pid Overlay with Double Gaussian and CB Fits", 800, 600);
    histEBPions->SetLineColor(kBlue);
    histPionHypothesis->SetLineColor(kRed);
    histEBPions->Draw();
    histPionHypothesis->Draw("SAME");

    // Fit with double Gaussian
    histPionHypothesis->Fit(doubleGausFit, "R+");
    double const1 = doubleGausFit->GetParameter(0);
    double mean1 = doubleGausFit->GetParameter(1);
    double sigma1 = doubleGausFit->GetParameter(2);
    double const2 = doubleGausFit->GetParameter(3);
    double mean2 = doubleGausFit->GetParameter(4);
    double sigma2 = doubleGausFit->GetParameter(5);

    // Fit with double Crystal Ball
    doubleCBFit->SetParameters(const1, mean1, sigma1, 1.5, 2.0, const2, mean2, sigma2, 1.5, 2.0);
    histPionHypothesis->Fit(doubleCBFit, "R+");

    // Individual Gaussians
    TF1* gausPion = new TF1("gausPion", "gaus", -5, 9);
    gausPion->SetParameters(doubleGausFit->GetParameter(0), doubleGausFit->GetParameter(1), doubleGausFit->GetParameter(2));
    gausPion->SetLineColor(kGreen);
    gausPion->SetLineStyle(kDashed);

    TF1* gausKaon = new TF1("gausKaon", "gaus", -5, 9);
    gausKaon->SetParameters(doubleGausFit->GetParameter(3), doubleGausFit->GetParameter(4), doubleGausFit->GetParameter(5));
    gausKaon->SetLineColor(kGreen);
    gausKaon->SetLineStyle(kDotted);

    // Individual Crystal Ball functions
    TF1* cbPion = new TF1("cbPion", "crystalball", -5, 9);
    cbPion->SetParameters(const1, mean1, sigma1, doubleCBFit->GetParameter(3), doubleCBFit->GetParameter(4));
    cbPion->SetLineColor(kBlack);
    cbPion->SetLineStyle(kDashed);

    TF1* cbKaon = new TF1("cbKaon", "crystalball", -5, 9);
    cbKaon->SetParameters(const2, mean2, sigma2, doubleCBFit->GetParameter(8), doubleCBFit->GetParameter(9));
    cbKaon->SetLineColor(kBlack);
    cbKaon->SetLineStyle(kDotted);

    // Draw components
    gausPion->Draw("SAME");
    gausKaon->Draw("SAME");
    cbPion->Draw("SAME");
    cbKaon->Draw("SAME");

    // Calculate areas
    const double sqrt2pi = sqrt(2 * TMath::Pi());
    double pionAreaGaus = doubleGausFit->GetParameter(0) * doubleGausFit->GetParameter(2) * sqrt2pi;
    double kaonAreaGaus = doubleGausFit->GetParameter(3) * doubleGausFit->GetParameter(5) * sqrt2pi;
    double totalAreaGaus = pionAreaGaus + kaonAreaGaus;
    double pionFractionGaus = totalAreaGaus ? pionAreaGaus / totalAreaGaus : 0;
    double kaonFractionGaus = totalAreaGaus ? kaonAreaGaus / totalAreaGaus : 0;

    double pionAreaCB = cbPion->Integral(-5, 9);
    double kaonAreaCB = cbKaon->Integral(-5, 9);
    double totalAreaCB = pionAreaCB + kaonAreaCB;
    double pionFractionCB = totalAreaCB ? pionAreaCB / totalAreaCB : 0;
    double kaonFractionCB = totalAreaCB ? kaonAreaCB / totalAreaCB : 0;

    // Add legend for canvas_1
    TLegend* legend_1 = new TLegend(0.1, 0.3, 0.3, 0.9);
    legend_1->SetTextSize(0.02);
    legend_1->SetBorderSize(1);
    legend_1->AddEntry(histEBPions, "EB Pions (pid == 211)", "l");
    legend_1->AddEntry(histPionHypothesis, "Manual Pion Hypothesis", "l");
    legend_1->AddEntry(doubleGausFit, "Double Gaussian Fit (Pion Hypothesis)", "l");
    legend_1->AddEntry(gausPion, "Gaussian (Pion Peak)", "l");
    legend_1->AddEntry(gausKaon, "Gaussian (Kaon Peak)", "l");
    legend_1->AddEntry((TObject*)0, Form("Const1: %.2f", doubleGausFit->GetParameter(0)), "");
    legend_1->AddEntry((TObject*)0, Form("Mean1: %.2f", doubleGausFit->GetParameter(1)), "");
    legend_1->AddEntry((TObject*)0, Form("Sigma1: %.2f", doubleGausFit->GetParameter(2)), "");
    legend_1->AddEntry((TObject*)0, Form("Const2: %.2f", doubleGausFit->GetParameter(3)), "");
    legend_1->AddEntry((TObject*)0, Form("Mean2: %.2f", doubleGausFit->GetParameter(4)), "");
    legend_1->AddEntry((TObject*)0, Form("Sigma2: %.2f", doubleGausFit->GetParameter(5)), "");
    double dGausChi2 = doubleGausFit->GetChisquare();
    double dGausNDF = doubleGausFit->GetNDF();
    legend_1->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.1f / %d = %.2f", dGausChi2, (int)dGausNDF, dGausChi2 / dGausNDF), "");
    legend_1->AddEntry((TObject*)0, Form("Pion Fraction (Gaus): %.3f", pionFractionGaus), "");
    legend_1->AddEntry((TObject*)0, Form("Kaon Fraction (Gaus): %.3f", kaonFractionGaus), "");
    legend_1->AddEntry(doubleCBFit, "Double Crystal Ball Fit (Pion Hypothesis)", "l");
    legend_1->AddEntry(cbPion, "Crystal Ball (Pion Peak)", "l");
    legend_1->AddEntry(cbKaon, "Crystal Ball (Kaon Peak)", "l");
    legend_1->AddEntry((TObject*)0, Form("Const1: %.2f", cbPion->GetParameter(0)), "");
    legend_1->AddEntry((TObject*)0, Form("Mean1: %.2f", cbPion->GetParameter(1)), "");
    legend_1->AddEntry((TObject*)0, Form("Sigma1: %.2f", cbPion->GetParameter(2)), "");
    legend_1->AddEntry((TObject*)0, Form("Alpha1: %.2f", cbPion->GetParameter(3)), "");
    legend_1->AddEntry((TObject*)0, Form("n1: %.2f", cbPion->GetParameter(4)), "");
    legend_1->AddEntry((TObject*)0, Form("Const2: %.2f", cbKaon->GetParameter(0)), "");
    legend_1->AddEntry((TObject*)0, Form("Mean2: %.2f", cbKaon->GetParameter(1)), "");
    legend_1->AddEntry((TObject*)0, Form("Sigma2: %.2f", cbKaon->GetParameter(2)), "");
    legend_1->AddEntry((TObject*)0, Form("Alpha2: %.2f", cbKaon->GetParameter(3)), "");
    legend_1->AddEntry((TObject*)0, Form("n2: %.2f", cbKaon->GetParameter(4)), "");
    double dCBChi2 = doubleCBFit->GetChisquare();
    double dCBNDF = doubleCBFit->GetNDF();
    legend_1->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.1f / %d = %.2f", dCBChi2, (int)dCBNDF, dCBChi2 / dCBNDF), "");
    legend_1->AddEntry((TObject*)0, Form("Pion Fraction (CB): %.3f", pionFractionCB), "");
    legend_1->AddEntry((TObject*)0, Form("Kaon Fraction (CB): %.3f", kaonFractionCB), "");
    legend_1->Draw();

    // Save chi2pid_overlay_1.pdf
    TString pdfFileName_2 = outputDir + "/chi2pid_overlay_1.pdf";
    gSystem->Unlink(pdfFileName_2);
    canvas_1->Update();
    canvas_1->Print(pdfFileName_2);

    // Canvas for Gaussian and CB fits with kaons
    TCanvas* canvas_2 = new TCanvas("canvas_2", "Chi2pid Overlay with Gaussian and CB Fits", 800, 600);
    histEBPions->SetLineColor(kBlue);
    histEBPions->SetFillColorAlpha(kBlue, 0.3);
    histPionHypothesis->SetLineColor(kBlack);
    histEBPions->Draw();
    histPionHypothesis->Draw("SAME");
    histPionHypothesis->Fit(gausKaon, "RN", "", 3, 8);
    gausKaon->SetRange(3, 8);
    gausKaon->SetLineStyle(1);
    gausKaon->SetLineColor(kBlack);
    gausKaon->SetLineWidth(1);
    histEBKaons->SetLineColor(kGreen);
    histEBKaons->Draw("SAME");

    TH1F* shiftedKaons = (TH1F*)histEBKaons->Clone("shiftedKaons");
    shiftedKaons->Reset();
    int nbins = histEBKaons->GetNbinsX();
    double shift = 5.10;
    for (int i = 1; i <= nbins; ++i) {
        double x_old = histEBKaons->GetBinCenter(i);
        double y = histEBKaons->GetBinContent(i);
        int new_bin = shiftedKaons->FindBin(x_old + shift);
        if (new_bin >= 1 && new_bin <= nbins)
            shiftedKaons->SetBinContent(new_bin, shiftedKaons->GetBinContent(new_bin) + y);
    }
    shiftedKaons->SetLineColor(kGreen + 2);
    shiftedKaons->SetFillColorAlpha(kGreen, 0.3);
    shiftedKaons->SetLineColor(kGreen);
    shiftedKaons->Draw("HIST SAME");

    // Save chi2pid_overlay_2.pdf
    TString pdfFileName_3 = outputDir + "/chi2pid_overlay_2.pdf";
    gSystem->Unlink(pdfFileName_3);
    canvas_2->Update();
    canvas_2->Print(pdfFileName_3);

    
    // Canvas for beta vs p (Page 1: Scatter plots before cuts)
// Canvas for beta vs p (Page 1: Scatter plots before cuts)
TCanvas* betaCanvas = new TCanvas("betaCanvas", "Beta", 800, 600);
histEBPions_beta_vs_p_bc->SetMarkerStyle(20);
histEBPions_beta_vs_p_bc->SetMarkerColor(kBlue);
histEBPions_beta_vs_p_bc->SetMarkerSize(1);
histEBPions_beta_vs_p_bc->SetFillColorAlpha(kBlue, 0.3);
histEBKaons_beta_vs_p_bc->SetMarkerStyle(20);
histEBKaons_beta_vs_p_bc->SetMarkerColor(kGreen + 2);
histEBKaons_beta_vs_p_bc->SetMarkerSize(1);
histEBKaons_beta_vs_p_bc->SetFillColorAlpha(kGreen, 0.3);
gPad->SetLogy();
histEBPions_beta_vs_p_bc->Draw("P2");
histEBKaons_beta_vs_p_bc->Draw("P2 SAME");
histEBPions_beta_vs_p_bc->SetTitle("Beta EB Pions and Kaons");

// Save chi2pid_overlay_3.pdf (Page 1)
TString pdfFileName_4 = outputDir + "/chi2pid_overlay_3.pdf";
gSystem->Unlink(pdfFileName_4);
betaCanvas->Update();
betaCanvas->Print(pdfFileName_4 + "[");
betaCanvas->Print(pdfFileName_4); // Write first page

// Second page for beta after cut with Gaussian fits
betaCanvas->Clear();
gPad->SetLogy(); // Reapply log scale for the second page
histEBPions_beta_vs_p_ac->SetMarkerStyle(20);
histEBPions_beta_vs_p_ac->SetMarkerColor(kBlue);
histEBPions_beta_vs_p_ac->SetMarkerSize(1);
histEBPions_beta_vs_p_ac->SetFillColorAlpha(kBlue, 0.3);
histEBKaons_beta_vs_p_ac->SetMarkerStyle(20);
histEBKaons_beta_vs_p_ac->SetMarkerColor(kGreen + 2);
histEBKaons_beta_vs_p_ac->SetMarkerSize(1);
histEBKaons_beta_vs_p_ac->SetFillColorAlpha(kGreen, 0.3);
histEBPions_beta_vs_p_ac->Draw("P2");
histEBKaons_beta_vs_p_ac->Draw("P2 SAME");
histEBPions_beta_vs_p_ac->SetTitle("Beta AC Pions and Kaons After Cut");

// Fit pions histogram with Gaussian
TF1* gausFitPions = new TF1("gausFitPions", "gaus", 0.984, 1.026); // Beta range around 1
gausFitPions->SetParameters(histEBPions_beta_vs_p_ac->GetMaximum(), 0, 1); // Initial guess: amplitude, mean, sigma
gausFitPions->SetLineColor(kBlue);
gausFitPions->SetLineStyle(2);
histEBPions_beta_vs_p_ac->Fit(gausFitPions, "R+");
double pionMean = gausFitPions->GetParameter(1);
double pionSigma = gausFitPions->GetParameter(2);
double pionChi2 = gausFitPions->GetChisquare();
double pionNDF = gausFitPions->GetNDF();

// Fit kaons histogram with Gaussian
TF1* gausFitKaons = new TF1("gausFitKaons", "gaus", 0.965, 0.995); // Beta range around 1
gausFitKaons->SetParameters(histEBKaons_beta_vs_p_ac->GetMaximum(), 0, 1); // Initial guess: amplitude, mean, sigma
gausFitKaons->SetLineColor(kGreen + 2);
gausFitKaons->SetLineStyle(2);
histEBKaons_beta_vs_p_ac->Fit(gausFitKaons, "R+");
double kaonMean = gausFitKaons->GetParameter(1);
double kaonSigma = gausFitKaons->GetParameter(2);
double kaonChi2 = gausFitKaons->GetChisquare();
double kaonNDF = gausFitKaons->GetNDF();

double overlapLow = TMath::Min(kaonMean, pionMean);
double overlapHigh = TMath::Max(kaonMean, pionMean);
// Integrate kaon Gaussian over the overlap region
double kaonIntegralOverlap = gausFitKaons->Integral(overlapLow, overlapHigh);

// Total integral of kaon Gaussian (over the full fit range)
double totalKaonIntegral = gausFitKaons->Integral(0.965, 0.995);

// Fraction of kaons in the overlap region
double kaonFraction = kaonIntegralOverlap / totalKaonIntegral;

// Total number of kaons from the histogram
double totalKaons = histEBKaons_beta_vs_p_ac->Integral();

// Number of kaons in the overlap region
double kaonsInOverlap = totalKaons * kaonFraction;

// Compute total pions from gausFitPions and scale to match histogram
double pionIntegral = gausFitPions->Integral(0.984, 1.026); // Corrected range
double totalPionsHist = histEBPions_beta_vs_p_ac->Integral();
double scaleFactorPions = totalPionsHist / pionIntegral;
double totalPions = pionIntegral * scaleFactorPions;

// Compute the ratio of kaons in overlap to total pions
double ratio = kaonsInOverlap / totalPions;

// Print the results
std::cout << "Kaons in overlap region: " << kaonsInOverlap << std::endl;
std::cout << "Total pions (from gausFitPions, scaled): " << totalPions << std::endl;
std::cout << "Ratio (kaons in overlap / total pions): " << ratio << std::endl;

// Add legend with fit parameters
TLegend* legendBeta = new TLegend(0.7, 0.5, 0.9, 0.9);
legendBeta->SetTextSize(0.02);
legendBeta->AddEntry(histEBPions_beta_vs_p_ac, "EB Pions Beta After Cut", "p");
legendBeta->AddEntry(gausFitPions, "Gaussian Fit (Pions)", "l");
legendBeta->AddEntry((TObject*)0, Form("Constant: %.2f", gausFitPions->GetParameter(0)), "");
legendBeta->AddEntry((TObject*)0, Form("Mean: %.4f", gausFitPions->GetParameter(1)), "");
legendBeta->AddEntry((TObject*)0, Form("Sigma: %.4f", gausFitPions->GetParameter(2)), "");
legendBeta->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.1f / %d = %.2f", pionChi2, (int)pionNDF, pionChi2 / pionNDF), "");
legendBeta->AddEntry(histEBKaons_beta_vs_p_ac, "EB Kaons Beta After Cut", "p");
legendBeta->AddEntry(gausFitKaons, "Gaussian Fit (Kaons)", "l");
legendBeta->AddEntry((TObject*)0, Form("Constant: %.2f", gausFitKaons->GetParameter(0)), "");
legendBeta->AddEntry((TObject*)0, Form("Mean: %.4f", gausFitKaons->GetParameter(1)), "");
legendBeta->AddEntry((TObject*)0, Form("Sigma: %.4f", gausFitKaons->GetParameter(2)), "");
legendBeta->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.1f / %d = %.2f", kaonChi2, (int)kaonNDF, kaonChi2 / kaonNDF), "");
legendBeta->Draw();

betaCanvas->Update();
betaCanvas->Print(pdfFileName_4); // Write second page
betaCanvas->Print(pdfFileName_4 + "]"); // Close PDF

// Clean up temporary objects
delete gausFitPions;
delete gausFitKaons;
    // Canvas for mass vs p (EB pions)
    TCanvas* mass = new TCanvas("mass", "massVsp for EB PID +ve pions", 1200, 600);
    mass->Divide(2, 1);
    mass->cd(1);
    gStyle->SetOptStat(1);
    histEBPions_mass_vs_p->Draw("COLZ");
    gPad->SetLogz();
    mass->cd(2);
    histEBPions_mass_vs_p_zoom->Draw("COLZ");
    gPad->SetLogz();

    // Save massVsP.pdf
    TString pdfFileName_1 = outputDir_1 + "/massVsP.pdf";
    gSystem->Unlink(pdfFileName_1);
    mass->Update();
    mass->Print(pdfFileName_1);

    // Canvas for mass vs p (pion hypothesis)
    TCanvas* massPH = new TCanvas("massPH", "massVsp for +ve pions hypothesis", 1200, 600);
    massPH->Divide(2, 1);
    massPH->cd(1);
    histPHPions_mass_vs_p->SetStats(kTRUE);
    histPHPions_mass_vs_p->Draw("COLZ");
    gPad->SetLogz();
    massPH->cd(2);
    histPHPions_mass_vs_p_zoom->Draw("COLZ");
    gPad->SetLogz();

    // Save massVsP_PH.pdf
    TString pdfFileName_PH = outputDir_1 + "/massVsP_PH.pdf";
    gSystem->Unlink(pdfFileName_PH);
    massPH->Update();
    massPH->Print(pdfFileName_PH);

    // Clean up
    delete histEBPions;
    delete histPionHypothesis;
    delete histEBKaons;
    delete histEBPions_beta_vs_p_bc;
    delete histEBKaons_beta_vs_p_bc;
    delete histEBPions_beta_vs_p_ac;
    delete histEBKaons_beta_vs_p_ac;
    delete histEBPions_mass_vs_p;
    delete histEBPions_mass_vs_p_zoom;
    delete histPHPions_mass_vs_p;
    delete histPHPions_mass_vs_p_zoom;
    delete gausFit;
    delete cb;
    delete doubleGausFit;
    delete doubleCBFit;
    delete gausPion;
    delete gausKaon;
    delete cbPion;
    delete cbKaon;
    delete shiftedKaons;
    file->Close();
}

int main() {
    chi2pid_310_345();
    return 0;
}