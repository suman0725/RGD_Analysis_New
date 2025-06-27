#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
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
    ROOT::RDataFrame df("EB_all_pion_assumed", file);

    // Step 3: Create output directories
    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/pdf_1.1", kTRUE);
    gSystem->mkdir("output/csv", kTRUE);

    // Step 4: Set up momentum bins (1 to 7 GeV, 20 bins)
    vector<double> pBins;
    double pMin = 0.4, pMax = 3.0;
    int nBins = 26;
    double binWidth = (pMax - pMin) / nBins;
    for (int i = 0; i <= nBins; i++) {
        pBins.push_back(pMin + i * binWidth);
    }

    // Step 6: Initialize canvases
    TCanvas* canvasChi2pid = new TCanvas("canvasChi2pid", "Chi2pid Fits (Before)", 1200, 800);
    canvasChi2pid->SetMargin(0.1, 0.05, 0.15, 0.1);

    canvasChi2pid->Print("output/pdf_1.1/chi2pid_gaus_321_2212_under_pionHypothesis.pdf[");

    // Step 7: Open CSV file
    ofstream csvFileCuts("output/csv/chi2pid_cuts.csv");
    csvFileCuts << "Momentum Bin (GeV/c),Pion Mean,Pion Sigma,Chi2 Min (Mean-5σ),Chi2 Max (Mean+5σ)\n";

    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];
    
        // Filter for momentum range and valid chi2pid
        auto filtered_df = df.Filter([pLow, pHigh](float p, float recomputed_chi2pid) {
            return p >= pLow && p < pHigh;
        }, {"p", "recomputed_chi2pid"});
    
        // Create histogram model for before cut
        ROOT::RDF::TH1DModel modelChi2pid(
            TString::Format("chi2_EBpions_before_%d", i),
            TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
            100, -15, 15
        );
    
        // Fill histogram
        auto histo = filtered_df.Histo1D(modelChi2pid, "recomputed_chi2pid");
        histo->SetStats(0); // Disable the stat box
    
        // Fit with Gaussian
        TF1* fitPion = new TF1("gausFit", "gaus", -1.5, 1.5);
        fitPion->SetParameters(590000, 0, 1);
        histo->Fit(fitPion, "R", "", -1.5, 1.5);
    
        // Extract fit parameters with validation
        double mean = fitPion->GetParameter(1);
        double sigma = fitPion->GetParameter(2);
        if (sigma <= 0 || isnan(mean) || isnan(sigma)) {
            cout << "Warning: Invalid fit for bin " << pLow << "-" << pHigh << ": mean = " << mean << ", sigma = " << sigma << ". Skipping 5σ lines." << endl;
            sigma = 1.0; // Default to 1 if invalid
        }
        double pionChi2Min = mean - 5 * sigma;
        double pionChi2Max = mean + 5 * sigma;
    
        // Save to CSV
        csvFileCuts << pLow << "-" << pHigh << "," << mean << "," << sigma << ","
                    << pionChi2Min << "," << pionChi2Max << "\n";
    
        // Plot before cut
        canvasChi2pid->Clear();
        histo->Draw("HIST");
        fitPion->SetLineColor(kRed);
        fitPion->SetLineWidth(1);
        fitPion->Draw("SAME");
    
        // Add vertical 5σ lines in blue (from y=0 to max count)
        double maxCount = histo->GetMaximum();
        TLine* lineMinus5Sigma = new TLine(pionChi2Min, 0, pionChi2Min, maxCount);
        TLine* linePlus5Sigma = new TLine(pionChi2Max, 0, pionChi2Max, maxCount);
        lineMinus5Sigma->SetLineColor(kBlack);
        linePlus5Sigma->SetLineColor(kBlack);
        lineMinus5Sigma->SetLineWidth(1); // Increased width for visibility
        linePlus5Sigma->SetLineWidth(1);
        lineMinus5Sigma->Draw("SAME");
        linePlus5Sigma->Draw("SAME");
    
        // Debug output
        cout << "Bin " << pLow << "-" << pHigh << ": Mean = " << mean << ", Sigma = " << sigma
             << ", 5σ Min = " << pionChi2Min << ", 5σ Max = " << pionChi2Max << endl;
    
        double chi2 = fitPion->GetChisquare();
        double ndf = fitPion->GetNDF();
        double chi2_ndf = (ndf != 0) ? chi2 / ndf : 0;
        TLegend* legPions = new TLegend(0.7, 0.6, 1.0, 0.9);
        legPions->AddEntry(histo.GetPtr(), "Pions", "l");
        legPions->AddEntry(fitPion, "Gaus Fit", "l");
        legPions->AddEntry(lineMinus5Sigma, "#mu #pm 5#sigma", "l"); // Added legend entry for 5σ lines

        legPions->AddEntry(histo.GetPtr(), TString::Format("Entries: %.0f", histo->GetEntries()), "");
        legPions->AddEntry((TObject*)0, TString::Format("A: %.4f #pm %.4f", fitPion->GetParameter(0), fitPion->GetParError(0)), "");
        legPions->AddEntry((TObject*)0, TString::Format("#mu: %.4f #pm %.4f", fitPion->GetParameter(1), fitPion->GetParError(1)), "");
        legPions->AddEntry((TObject*)0, TString::Format("#sigma: %.4f #pm %.4f", fitPion->GetParameter(2), fitPion->GetParError(2)), "");
        legPions->AddEntry((TObject*)0, TString::Format("#chi^{2}/ndf = %.4f / %.0f = %.4f", chi2, ndf, chi2_ndf), "");
        legPions->SetBorderSize(0);
        legPions->SetFillStyle(0);
        legPions->SetTextSize(0.02);
        legPions->Draw();
        canvasChi2pid->Update();
        canvasChi2pid->Print("output/pdf_1.1/chi2pid_gaus_321_2212_under_pionHypothesis.pdf");
        delete legPions;
        delete lineMinus5Sigma;
        delete linePlus5Sigma;
    
        // Clean up fit
        delete fitPion;
    }

    // Step 9: Close PDF and CSV files
    canvasChi2pid->Print("output/pdf_1.1/chi2pid_gaus_321_2212_under_pionHypothesis.pdf]");
    csvFileCuts.close();

    // Step 10: Clean up
    delete canvasChi2pid;
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