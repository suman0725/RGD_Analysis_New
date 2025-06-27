#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>  // Added for 2D histograms
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <vector>
#include <iostream>

void plot_chi2pid() {
    // Enable batch mode to avoid graphical display issues
    gROOT->SetBatch(kTRUE);

    // Open the ROOT file
    TFile* file = new TFile("piontreeCxC.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open piontree.root" << std::endl;
        return;
    }

    // Read the trees
    TTree* treeEBPions = (TTree*)file->Get("EB_pid_pions");
    TTree* treePionHypothesis = (TTree*)file->Get("pion_hypothesis");
    if (!treeEBPions || !treePionHypothesis) {
        std::cerr << "Error: Cannot find one or both trees in piontree.root" << std::endl;
        file->Close();
        return;
    }

    // Define momentum bins (1 to 8 GeV, 0.35 GeV bin size)
    const int nBins = 20; // (8 - 1) / 0.35 = 20 bins
    float pBins[nBins + 1];
    for (int i = 0; i <= nBins; i++) {
        pBins[i] = 1.0 + i * 0.35; // 1.00, 1.35, 1.70, ..., 8.00
    }

    // Define chi2pid histogram range
    const int nChi2Bins = 100;
    const float chi2Min = -15.0;
    const float chi2Max = 15.0;

    // Vectors to store 1D histograms for each tree
    std::vector<TH1F*> histEBPions(nBins); 
    std::vector<TH1F*> histPionHypothesis(nBins);

    // Create canvases for each momentum bin (for 1D histograms)
    std::vector<TCanvas*> canvases(nBins);

    // Set up the trees to read chi2pid and p
    float chi2pidEBPions, pEBPions;
    float chi2pidPionHypothesis, pPionHypothesis;
    treeEBPions->SetBranchAddress("chi2pid", &chi2pidEBPions);
    treeEBPions->SetBranchAddress("p", &pEBPions);
    treePionHypothesis->SetBranchAddress("chi2pid", &chi2pidPionHypothesis);
    treePionHypothesis->SetBranchAddress("p", &pPionHypothesis);

    // Create 1D histograms for each momentum bin
    for (int i = 0; i < nBins; i++) {
        histEBPions[i] = new TH1F(Form("h_chi2pid_ebpions_bin%d", i),
                                  Form("Chi2pid (%.2f < p < %.2f GeV)", pBins[i], pBins[i+1]),
                                  nChi2Bins, chi2Min, chi2Max);
        histPionHypothesis[i] = new TH1F(Form("h_chi2pid_pionhypothesis_bin%d", i),
                                         Form("Chi2pid (%.2f < p < %.2f GeV)", pBins[i], pBins[i+1]),
                                         nChi2Bins, chi2Min, chi2Max);

        // Style the 1D histograms
        histEBPions[i]->SetLineColor(kBlue);
        histEBPions[i]->SetFillColorAlpha(kBlue, 0.3);
        histEBPions[i]->SetXTitle("chi2pid (EB pid 211)");
        histEBPions[i]->SetYTitle("Entries");

        histPionHypothesis[i]->SetLineColor(kRed);
        histPionHypothesis[i]->SetFillColorAlpha(kRed, 0.3);
        histPionHypothesis[i]->SetXTitle("chi2pid (Pion Hypothesis)");
        histPionHypothesis[i]->SetYTitle("Entries");
    }

    // Fill 1D histograms for EB pions (pid == 211)
    Long64_t nEntriesEBPions = treeEBPions->GetEntries();
    for (Long64_t j = 0; j < nEntriesEBPions; j++) {
        treeEBPions->GetEntry(j);
        for (int i = 0; i < nBins; i++) {
            if (pEBPions >= pBins[i] && pEBPions < pBins[i+1]) {
                 if (chi2pidEBPions != 99999.0f)histEBPions[i]->Fill(chi2pidEBPions);
                break;
            }
        }
    }

    // Fill 1D histograms for particles under manual pion hypothesis
    Long64_t nEntriesPionHypothesis = treePionHypothesis->GetEntries();
    for (Long64_t j = 0; j < nEntriesPionHypothesis; j++) {
        treePionHypothesis->GetEntry(j);
        for (int i = 0; i < nBins; i++) {
            if (pPionHypothesis >= pBins[i] && pPionHypothesis < pBins[i+1]) {
                if(chi2pidPionHypothesis != 99999.0f) histPionHypothesis[i]->Fill(chi2pidPionHypothesis);
                break;
            }
        }
    }

    // Create canvases and overlay 1D histograms
    for (int i = 0; i < nBins; i++) {
        canvases[i] = new TCanvas(Form("canvas_bin%d", i),
                                  Form("Chi2pid for %.2f < p < %.2f GeV", pBins[i], pBins[i+1]),
                                  800, 600);

        // Normalize 1D histograms to unit area for better comparison (optional)
       // if (histEBPions[i]->Integral() > 0) histEBPions[i]->Scale(1.0 / histEBPions[i]->Integral());
        //if (histPionHypothesis[i]->Integral() > 0) histPionHypothesis[i]->Scale(1.0 / histPionHypothesis[i]->Integral());

        // Draw the 1D histograms overlaid
        histEBPions[i]->Draw("");
        histPionHypothesis[i]->Draw("SAME");

        // Set the y-axis range to fit both 1D histograms
        float maxY = std::max(histEBPions[i]->GetMaximum(), histPionHypothesis[i]->GetMaximum());
        histEBPions[i]->SetMaximum(maxY * 1.2); // Add 20% headroom

        // Add a legend with updated labels
        TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg->AddEntry(histEBPions[i], "EB Pions (pid == 211)", "f");
        leg->AddEntry(histPionHypothesis[i], "Manual Pion Hypothesis", "f");
        leg->Draw();
    }

    // Save all 1D histogram canvases to a single PDF file
    canvases[0]->Print("chi2pid_overlaid.pdf["); // Open PDF
    for (int i = 0; i < nBins; i++) {
        canvases[i]->Print("chi2pid_overlaid.pdf"); // Add each canvas to PDF
    }
    canvases[nBins-1]->Print("chi2pid_overlaid.pdf]"); // Close PDF

    // --- Add 2D Histograms for chi2pid vs p ---

    // Define 2D histogram ranges
    const int nPBins = 20; // Momentum bins from 1 to 8 GeV
    const float pMin = 1.0;
    const float pMax = 8.0;

    // Create 2D histograms
    TH2F* hist2DEBPions = new TH2F("h_chi2pid_vs_p_ebpions", "EB Pions (pid == 211): chi2pid vs p",
                                   nPBins, pMin, pMax, nChi2Bins, chi2Min, chi2Max);
    TH2F* hist2DPionHypothesis = new TH2F("h_chi2pid_vs_p_pionhypothesis", "Manual Pion Hypothesis: chi2pid vs p",
                                          nPBins, pMin, pMax, nChi2Bins, chi2Min, chi2Max);

    // Fill 2D histogram for EB pions (pid == 211)
    for (Long64_t j = 0; j < nEntriesEBPions; j++) {
        treeEBPions->GetEntry(j);
        hist2DEBPions->Fill(pEBPions, chi2pidEBPions);
    }

    // Fill 2D histogram for particles under manual pion hypothesis
    for (Long64_t j = 0; j < nEntriesPionHypothesis; j++) {
        treePionHypothesis->GetEntry(j);
        hist2DPionHypothesis->Fill(pPionHypothesis, chi2pidPionHypothesis);
    }

    // Create canvases and draw 2D histograms
    TCanvas* canvas2DEBPions = new TCanvas("canvas_2d_ebpions", "EB Pions (pid == 211): chi2pid vs p", 800, 600);
    hist2DEBPions->SetXTitle("Momentum p (GeV)");
    hist2DEBPions->SetYTitle("chi2pid (EB pid 211)");
    hist2DEBPions->Draw("COLZ");

    TCanvas* canvas2DPionHypothesis = new TCanvas("canvas_2d_pionhypothesis", "Manual Pion Hypothesis: chi2pid vs p", 800, 600);
    hist2DPionHypothesis->SetXTitle("Momentum p (GeV)");
    hist2DPionHypothesis->SetYTitle("chi2pid (Pion Hypothesis)");
    hist2DPionHypothesis->Draw("COLZ");

    // Save 2D histograms to a separate PDF file
    
    canvas2DEBPions->Print("chi2pid_vs_p_2D_EBpion.pdf");  // Add EB Pions 2D plot
    canvas2DPionHypothesis->Print("chi2pid_vs_p_2D_pionhypo.pdf"); // Add Manual Pion Hypothesis 2D plot
   

    // Clean up
    file->Close();
    delete file;
}

int main() {
    plot_chi2pid();
    return 0;
}