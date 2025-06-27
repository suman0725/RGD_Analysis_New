// using double crystal ball function for chipid and beta
// change output/PID_Cuts_P_P/pdf to output/PID_Cuts_P_P/pdf_1.1

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

int main()
{
    // Step 1: Open the ROOT file and load the trees for pions and kaons
    TFile *file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification//Skim/pkptreeCxC_9_test_modified.root");
    // TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification//Skim/pkptreeCxC_9_test.root");
    if (!file || file->IsZombie())
    {
        cerr << "Error: Cannot open file pkptreeCxC_9_test_modified.root\n";
        return 1;
    }
    TTree *treePions = (TTree *)file->Get("EB_pid_pions");

    if (!treePions)
    {
        cerr << "Error: Cannot load trees\n";
        file->Close();
        delete file;
        return 1;
    }

    // Step 2: Create output directories for saving PDFs and CSVs
    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/PID_Cuts_P_P/pdf", kTRUE);
    gSystem->mkdir("output/PID_Cuts_P_P/csv", kTRUE);

    // Step 3: Set up momentum bins (1 to 7 GeV, 20 bins)
    vector<double> pBins;
    double pMin = 1.0, pMax = 7.0;
    int nBins = 20;
    double binWidth = (pMax - pMin) / nBins;
    for (int i = 0; i <= nBins; i++)
        pBins.push_back(pMin + i * binWidth);

    // Step 4: Set chi2pid histogram ranges, fit ranges, and fit parameters for each bin
    vector<pair<double, double>> chi2HistRanges(nBins);
    vector<pair<double, double>> chi2FitRanges(nBins);
    vector<tuple<double, double, double>> chi2FitParams(nBins);
    for (int i = 0; i < nBins; i++)
    {
        double pLow = pBins[i];
        if (pLow < 2.0)
        {
            chi2HistRanges[i] = {-5.0, 5.0};
            chi2FitRanges[i] = {-3, 3};
            chi2FitParams[i] = {1000.0, 0.0, 1.0};
        }
        else if (pLow < 3.0)
        {
            chi2HistRanges[i] = {-4.0, 4.0};
            chi2FitRanges[i] = {-1.5, 1.5};
            chi2FitParams[i] = {800.0, 0.1, 0.8};
        }
        else
        {
            chi2HistRanges[i] = {-3.0, 3.0};
            chi2FitRanges[i] = {-1, 1};
            chi2FitParams[i] = {500.0, -0.1, 0.5};
        }
    }

    TCanvas *canvasChi2Before = new TCanvas("canvasChi2Before", "Chi2pid Fits (Before)", 1200, 800);
    canvasChi2Before->SetMargin(0.1, 0.05, 0.15, 0.1);
    TCanvas *canvasChi2After_bin_1 = new TCanvas("canvasChi2After_bin_1", "Chi2pid Fits (After)", 1200, 800);
    canvasChi2After_bin_1->SetMargin(0.1, 0.05, 0.15, 0.1);
    TCanvas *canvasChi2After_bin_all = new TCanvas("canvasChi2After_bin_all", "Chi2pid (After)", 1200, 800);
    canvasChi2After_bin_all->SetMargin(0.1, 0.05, 0.15, 0.1);

    canvasChi2Before->Print("output/PID_Cuts_P_P/pdf/chi2pid_fits_linear_before.pdf[");
    canvasChi2After_bin_1->Print("output/PID_Cuts_P_P/pdf/chi2pid_after_bin_1.pdf[");
    canvasChi2After_bin_all->Print("output/PID_Cuts_P_P/pdf/chi2pid_after_bin_all.pdf[");
    // Step 8: Open CSV files to save results
    ofstream csvFileCuts("output/PID_Cuts_P_P/csv/chi2pid_cuts.csv");
    csvFileCuts << "Momentum Bin (GeV/c),Pion Mean,Pion Sigma,Chi2 Min (Mean-5σ),Chi2 Max (Mean+5σ)\n";

    // Step 9: Loop over each momentum bin to process data
    for (int i = 0; i < nBins; i++)
    {
        double pLow = pBins[i], pHigh = pBins[i + 1];

        // --- Before Chi2pid Cut ---
        // Create chi2pid histograms for pions and kaons (before cut)
        TH1F *chi2PionsBefore = new TH1F(TString::Format("chi2_EBpions_before_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                         100, -5, 5);

        // Fill chi2pid histograms
        float p, chi2pid, beta, orig_chi2pid, recomputed_chi2pid;
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("orig_chi2pid", &orig_chi2pid);
        for (Long64_t j = 0; j < treePions->GetEntries(); j++)
        {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh && orig_chi2pid != 9999.0)
                chi2PionsBefore->Fill(orig_chi2pid);
        }

        // Fit chi2pid histograms with a Gaussian to determine the ±5σ cut
        double pionMeanBefore, pionSigmaBefore, pionChi2Min, pionChi2Max, pionAmpBefore, pionChi2Min_bin_1, pionChi2Max_bin_1;
        auto [pionAmp, pionMeanInit, pionSigmaInit] = chi2FitParams[i];

        // commented out this gaussian fit to chi2pid pions before cut
        TF1 *fitPion = new TF1("gausFit", "gaus", chi2FitRanges[i].first, chi2FitRanges[i].second);
        fitPion->SetParameters(pionAmp, pionMeanInit, pionSigmaInit);
        chi2PionsBefore->Fit(fitPion, "R", "", chi2FitRanges[i].first, chi2FitRanges[i].second);
        pionMeanBefore = fitPion->GetParameter(1);
        pionSigmaBefore = fitPion->GetParameter(2);
        if (i == 0)
        {
            pionChi2Min_bin_1 = pionMeanBefore - 5 * pionSigmaBefore;
            pionChi2Max_bin_1 = pionMeanBefore + 5 * pionSigmaBefore;
        }
        pionChi2Min = pionMeanBefore - 5 * pionSigmaBefore;
        pionChi2Max = pionMeanBefore + 5 * pionSigmaBefore;

        // Save chi2pid fit parameters to CSV
        csvFileCuts << pLow << "-" << pHigh << "," << pionMeanBefore << "," << pionSigmaBefore << ","
                    << pionChi2Min << "," << pionChi2Max << "\n";

        canvasChi2Before->Clear();
        chi2PionsBefore->Draw("HIST");
        fitPion->SetLineColor(kBlue);
        fitPion->SetLineWidth(1);
        fitPion->Draw("SAME");
        double chi2 = fitPion->GetChisquare();
        double ndf = fitPion->GetNDF();
        double chi2_ndf = (ndf != 0) ? chi2 / ndf : 0;
        TLegend *legPions = new TLegend(0.6, 0.6, 0.9, 0.9);
        legPions->AddEntry(chi2PionsBefore, "Pions", "l");
        legPions->AddEntry(fitPion, "Gaus Fit", "l");
        legPions->AddEntry((TObject *)0, TString::Format("A_{left}: %.2f", fitPion->GetParameter(0)), "");
        legPions->AddEntry((TObject *)0, TString::Format("#mu: %.2f", fitPion->GetParameter(1)), "");
        legPions->AddEntry((TObject *)0, TString::Format("#sigma: %.2f", fitPion->GetParameter(2)), "");
        legPions->AddEntry((TObject *)0, TString::Format("#chi^{2}/ndf = %.1f / %.0f = %.2f ", chi2_ndf, chi2, ndf), "");
        legPions->SetBorderSize(0);
        legPions->SetFillStyle(0);
        legPions->SetTextSize(0.02);
        legPions->Draw();

        canvasChi2Before->Update();
        canvasChi2Before->Print("output/PID_Cuts_P_P/pdf/chi2pid_fits_linear_before.pdf");
        delete legPions;

        // --- After Chi2pid Cut (±5σ) ---

        TH1F *chi2PionsAfter_bin_1 = new TH1F(TString::Format("chi2_pions_after%d", i),
                                              TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                              100, -5, 5);
        TH1F *chi2PionsAfter_bin_all = new TH1F(TString::Format("chi2_pions_after_%d", i),
                                                TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                                100, -5, 5);
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("orig_chi2pid", &orig_chi2pid);
        for (Long64_t j = 0; j < treePions->GetEntries(); j++)
        {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh && orig_chi2pid >= pionChi2Min_bin_1 && orig_chi2pid <= pionChi2Max_bin_1)
            {
                chi2PionsAfter_bin_1->Fill(orig_chi2pid);
            }
            if (p >= pLow && p < pHigh && orig_chi2pid >= pionChi2Min && orig_chi2pid <= pionChi2Max)
            {
                chi2PionsAfter_bin_all->Fill(orig_chi2pid);
            }
        }

        canvasChi2After_bin_1->Clear();
        chi2PionsAfter_bin_1->Draw("HIST");
        TLegend *legPionsAfter_bin_1 = new TLegend(0.1, 0.75, 0.3, 0.9);
        legPionsAfter_bin_1->AddEntry(chi2PionsAfter_bin_1, "Pions", "l");
        legPionsAfter_bin_1->SetBorderSize(0);
        legPionsAfter_bin_1->SetFillStyle(0);
        legPionsAfter_bin_1->SetTextSize(0.025);
        legPionsAfter_bin_1->Draw();
        canvasChi2After_bin_1->Update();
        canvasChi2After_bin_1->Print("output/PID_Cuts_P_P/pdf/chi2pid_after_bin_1.pdf");
        delete legPionsAfter_bin_1;

        canvasChi2After_bin_all->Clear();
        chi2PionsAfter_bin_all->Draw("HIST");
        TLegend *legPionsAfter_bin_all = new TLegend(0.1, 0.75, 0.3, 0.9);
        legPionsAfter_bin_all->AddEntry(chi2PionsAfter_bin_1, "Pions", "l");
        legPionsAfter_bin_all->SetBorderSize(0);
        legPionsAfter_bin_all->SetFillStyle(0);
        legPionsAfter_bin_all->SetTextSize(0.025);
        legPionsAfter_bin_all->Draw();
        canvasChi2After_bin_all->Update();
        canvasChi2After_bin_all->Print("output/PID_Cuts_P_P/pdf/chi2pid_after_bin_all.pdf");
        delete legPionsAfter_bin_all;
     

        // Clean up histograms and fits
        delete chi2PionsBefore;
        delete chi2PionsAfter_bin_1;
        delete chi2PionsAfter_bin_all;
    }

    // Step 10: Close PDF and CSV files

    canvasChi2Before->Print("output/PID_Cuts_P_P/pdf/chi2pid_fits_linear_before.pdf]");
    canvasChi2After_bin_1->Print("output/PID_Cuts_P_P/pdf/chi2pid_after_bin_1.pdf]");
    canvasChi2After_bin_all->Print("output/PID_Cuts_P_P/pdf/chi2pid_after_bin_all.pdf]");
    csvFileCuts.close();

    // Step 11: Clean up

    delete canvasChi2Before;
    delete canvasChi2After_bin_1;
    file->Close();
    delete file;

    return 0;
}