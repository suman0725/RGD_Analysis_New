#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

using namespace std;
using namespace ROOT;

int main() {

     // Start total program timer
     auto start = chrono::high_resolution_clock::now();
    // Step 1: Enable multi-threading
    ROOT::EnableImplicitMT(4); // Adjust based on your CPU (e.g., 4, 8)
    cout << "Multi-threading enabled for RDataFrame processing." << endl;

    // Step 2: Open the ROOT file and create RDataFrame
    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/Skim/pkptreeCxC_9_test_modified.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file pkptreeCxC_9_test_modified.root\n";
        ROOT::DisableImplicitMT();
        return 1;
    }
    ROOT::RDataFrame df("EB_pid_pions", file);

    // Step 3: Create output directories
    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/PID_Cuts_P_P/pdf_1.1", kTRUE);
    gSystem->mkdir("output/PID_Cuts_P_P/csv", kTRUE);

    // Step 4: Set up momentum bins (1 to 7 GeV, 20 bins)
    vector<double> pBins;
    double pMin = 0.1, pMax = 7.0;
    int nBins = 23;
    double binWidth = (pMax - pMin) / nBins;
    for (int i = 0; i <= nBins; i++) {
        pBins.push_back(pMin + i * binWidth);
    }

    // Step 5: Set chi2pid histogram ranges, fit ranges, and initial fit parameters
    vector<pair<double, double>> chi2HistRanges(nBins);
    vector<pair<double, double>> chi2FitRanges(nBins);
    vector<tuple<double, double, double>> chi2FitParams(nBins); // For Gaussian: amplitude, mean, sigma
    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i];
        if (pLow <= 3.10) { // For bins [1.0-2.0), [2.0-3.0), [3.0-4.0), [4.0-5.0), [5.0-6.0){
            chi2HistRanges[i] = {-5.0, 5.0};
            chi2FitRanges[i] = {-1.6, 1.6};
            chi2FitParams[i] = {167283.78, 0.0, 1.0};
        } else if(pLow <=3.40) {
            chi2HistRanges[i] = {-5.0, 5.0};
            chi2FitRanges[i] = {-1.6, 1.0};
            chi2FitParams[i] = {500.0, -0.5, 0.5};
        }else if(pLow <=3.40) {
            chi2HistRanges[i] = {-5.0, 5.0};
            chi2FitRanges[i] = {-1.6, 1.0};
            chi2FitParams[i] = {500.0, -0.5, 0.5};
        }else if(pLow <=3.70) {
            chi2HistRanges[i] = {-5.0, 5.0};
            chi2FitRanges[i] = {-1.6, 0.8};
            chi2FitParams[i] = {500.0, -0.5, 0.5};
        }else if(pLow <=4) {
            chi2HistRanges[i] = {-5.0, 5.0};
            chi2FitRanges[i] = {-1.6, 0.4};
            chi2FitParams[i] = {500.0, -0.5, 0.5};
        }else if(pLow <=4.9) {
            chi2HistRanges[i] = {-5.0, 5.0};
            chi2FitRanges[i] = {-1.6, 0.5};
            chi2FitParams[i] = {500.0, -0.5, 0.5};
        }
        else {
            chi2HistRanges[i] = {-5.0, 5.0};
            chi2FitRanges[i] = {-0.5, 0.5};
            chi2FitParams[i] = {500.0, -0.1, 0.5};
        }
    }

    // Step 6: Initialize canvases
    TCanvas* canvasChi2Before = new TCanvas("canvasChi2Before", "Chi2pid Fits (Before)", 1200, 800);
    canvasChi2Before->SetMargin(0.1, 0.05, 0.15, 0.1);
    TCanvas* canvasChi2After_bin_1 = new TCanvas("canvasChi2After_bin_1", "Chi2pid Fits (After, Bin 1)", 1200, 800);
    canvasChi2After_bin_1->SetMargin(0.1, 0.05, 0.15, 0.1);
    TCanvas* canvasChi2After_bin_all = new TCanvas("canvasChi2After_bin_all", "Chi2pid (After, All Bins)", 1200, 800);
    canvasChi2After_bin_all->SetMargin(0.1, 0.05, 0.15, 0.1);

    canvasChi2Before->Print("output/PID_Cuts_P_P/pdf_1.1/chi2pid_fits_linear_before.pdf[");
    canvasChi2After_bin_1->Print("output/PID_Cuts_P_P/pdf_1.1/chi2pid_after_bin_1.pdf[");
    canvasChi2After_bin_all->Print("output/PID_Cuts_P_P/pdf_1.1/chi2pid_after_bin_all.pdf[");

    // Step 7: Open CSV file
    ofstream csvFileCuts("output/PID_Cuts_P_P/csv/chi2pid_cuts.csv");
    csvFileCuts << "Momentum Bin (GeV/c),Pion Mean,Pion Sigma,Chi2 Min (Mean-5σ),Chi2 Max (Mean+5σ)\n";

    // Step 8: Process each momentum bin
    double pionChi2Min_bin_1 = 0, pionChi2Max_bin_1 = 0;
    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];

        // Filter for momentum range and valid chi2pid
        auto filtered_df = df.Filter([pLow, pHigh](float p, float orig_chi2pid) {
            return p >= pLow && p < pHigh ;
        }, {"p", "orig_chi2pid"});

        // Create histogram model for before cut
        ROOT::RDF::TH1DModel modelChi2PionsBefore(
            TString::Format("chi2_EBpions_before_%d", i),
            TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
            100, chi2HistRanges[i].first, chi2HistRanges[i].second
        );

        // Fill histogram
        auto histo = filtered_df.Histo1D(modelChi2PionsBefore, "orig_chi2pid");

        // Fit with Gaussian
        auto [pionAmp, pionMeanInit, pionSigmaInit] = chi2FitParams[i];
        TF1* fitPion = new TF1("gausFit", "gaus", chi2FitRanges[i].first, chi2FitRanges[i].second);
        fitPion->SetParameters(pionAmp, pionMeanInit, pionSigmaInit);
        histo->Fit(fitPion, "R", "", chi2FitRanges[i].first, chi2FitRanges[i].second);

        // Extract fit parameters
        double pionMeanBefore = fitPion->GetParameter(1);
        double pionSigmaBefore = fitPion->GetParameter(2);
        double pionChi2Min = pionMeanBefore - 5 * pionSigmaBefore;
        double pionChi2Max = pionMeanBefore + 5 * pionSigmaBefore;
        if (i == 0) {
            pionChi2Min_bin_1 = pionChi2Min;
            pionChi2Max_bin_1 = pionChi2Max;
        }

        // Save to CSV
        csvFileCuts << pLow << "-" << pHigh << "," << pionMeanBefore << "," << pionSigmaBefore << ","
                    << pionChi2Min << "," << pionChi2Max << "\n";

        // Plot before cut
        canvasChi2Before->Clear();
        histo->Draw("HIST");
        fitPion->SetLineColor(kBlue);
        fitPion->SetLineWidth(1);
        fitPion->SetLineColor(kRed);
        fitPion->Draw("SAME");
        double chi2 = fitPion->GetChisquare();
        double ndf = fitPion->GetNDF();
        double chi2_ndf = (ndf != 0) ? chi2 / ndf : 0;
        TLegend* legPions = new TLegend(0.6, 0.6, 0.9, 0.9);
        legPions->AddEntry(histo.GetPtr(), "Pions", "l");
        legPions->AddEntry(fitPion, "Gaus Fit", "l");
        legPions->AddEntry((TObject*)0, TString::Format("A: %.4f #pm %.4f", fitPion->GetParameter(0), fitPion->GetParError(0)), "");
        legPions->AddEntry((TObject*)0, TString::Format("#mu: %.4f #pm %.4f", fitPion->GetParameter(1), fitPion->GetParError(1)), "");
        legPions->AddEntry((TObject*)0, TString::Format("#sigma: %.4f #pm %.4f", fitPion->GetParameter(2), fitPion->GetParError(2)), "");
        
        legPions->AddEntry((TObject*)0, TString::Format("#chi^{2}/ndf = %.4f / %.0f = %.4f", chi2, ndf, chi2_ndf), "");
        legPions->SetBorderSize(0);
        legPions->SetFillStyle(0);
        legPions->SetTextSize(0.02);
        legPions->Draw();
        canvasChi2Before->Update();
        canvasChi2Before->Print("output/PID_Cuts_P_P/pdf_1.1/chi2pid_fits_linear_before.pdf");
        delete legPions;

        // After chi2pid cut (±5σ)
        // Filter for bin 1 cut
        auto filtered_after_bin_1 = df.Filter([pLow, pHigh, pionChi2Min_bin_1, pionChi2Max_bin_1](float p, float orig_chi2pid) {
            return p >= pLow && p < pHigh && orig_chi2pid >= pionChi2Min_bin_1 && orig_chi2pid <= pionChi2Max_bin_1 ;
        }, {"p", "orig_chi2pid"});

        // Filter for individual bin cut
        auto filtered_after_bin_all = df.Filter([pLow, pHigh, pionChi2Min, pionChi2Max](float p, float orig_chi2pid) {
            return p >= pLow && p < pHigh && orig_chi2pid >= pionChi2Min && orig_chi2pid <= pionChi2Max ;
        }, {"p", "orig_chi2pid"});

        // Histogram models after cuts
        ROOT::RDF::TH1DModel modelChi2PionsAfter_bin_1(
            TString::Format("chi2_pions_after_%d", i),
            TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
            100, chi2HistRanges[i].first, chi2HistRanges[i].second
        );
        ROOT::RDF::TH1DModel modelChi2PionsAfter_bin_all(
            TString::Format("chi2_pions_after_all_%d", i),
            TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
            100, chi2HistRanges[i].first, chi2HistRanges[i].second
        );

        // Fill histograms
        auto histo_after_bin_1 = filtered_after_bin_1.Histo1D(modelChi2PionsAfter_bin_1, "orig_chi2pid");
        auto histo_after_bin_all = filtered_after_bin_all.Histo1D(modelChi2PionsAfter_bin_all, "orig_chi2pid");

        // Plot after cut (bin 1)
        canvasChi2After_bin_1->Clear();
        histo_after_bin_1->Draw("HIST");
        TLegend* legPionsAfter_bin_1 = new TLegend(0.1, 0.75, 0.3, 0.9);
        legPionsAfter_bin_1->AddEntry(histo_after_bin_1.GetPtr(), "Pions", "l");
        legPionsAfter_bin_1->SetBorderSize(0);
        legPionsAfter_bin_1->SetFillStyle(0);
        legPionsAfter_bin_1->SetTextSize(0.025);
        legPionsAfter_bin_1->Draw();
        canvasChi2After_bin_1->Update();
        canvasChi2After_bin_1->Print("output/PID_Cuts_P_P/pdf_1.1/chi2pid_after_bin_1.pdf");
        delete legPionsAfter_bin_1;

        // Plot after cut (all bins)
        canvasChi2After_bin_all->Clear();
        histo_after_bin_all->Draw("HIST");
        TLegend* legPionsAfter_bin_all = new TLegend(0.1, 0.75, 0.3, 0.9);
        legPionsAfter_bin_all->AddEntry(histo_after_bin_all.GetPtr(), "Pions", "l");
        legPionsAfter_bin_all->SetBorderSize(0);
        legPionsAfter_bin_all->SetFillStyle(0);
        legPionsAfter_bin_all->SetTextSize(0.025);
        legPionsAfter_bin_all->Draw();
        canvasChi2After_bin_all->Update();
        canvasChi2After_bin_all->Print("output/PID_Cuts_P_P/pdf_1.1/chi2pid_after_bin_all.pdf");
        delete legPionsAfter_bin_all;

        // Clean up fit
        delete fitPion;
    }

    // Step 9: Close PDF and CSV files
    canvasChi2Before->Print("output/PID_Cuts_P_P/pdf_1.1/chi2pid_fits_linear_before.pdf]");
    canvasChi2After_bin_1->Print("output/PID_Cuts_P_P/pdf_1.1/chi2pid_after_bin_1.pdf]");
    canvasChi2After_bin_all->Print("output/PID_Cuts_P_P/pdf_1.1/chi2pid_after_bin_all.pdf]");
    csvFileCuts.close();

    // Step 10: Clean up
    delete canvasChi2Before;
    delete canvasChi2After_bin_1;
    delete canvasChi2After_bin_all;
    file->Close();
    delete file;

    // Disable multi-threading
    ROOT::DisableImplicitMT();
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Total program time: " << duration.count() / 1000.0 << " seconds" << endl;


    cout << "Program completed successfully." << endl;
    return 0;
}