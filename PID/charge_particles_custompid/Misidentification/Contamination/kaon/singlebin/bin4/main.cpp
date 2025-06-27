#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLegend.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <TRandom.h>

using namespace std;

// Double Crystal Ball function
Double_t doubleCrystalBall(Double_t* x, Double_t* par) {
    Double_t t = (x[0] - par[1]) / par[2]; // t = (x - mean) / sigma
    Double_t absAlphaLeft = TMath::Abs(par[3]);
    Double_t absAlphaRight = TMath::Abs(par[6]);
/* par[0]: Amplitude of the left Gaussian/tail ($ A_{\text{left}} $)
par[1]: Mean ($ \mu $)
par[2]: Width ($ \sigma $)
par[3]: Left tail parameter ($ \alpha_{\text{left}} $)
par[4]: Left tail exponent ($ n_{\text{left}} $)
par[5]: Amplitude of the right Gaussian/tail ($ A_{\text{right}} $)
par[6]: Right tail parameter ($ \alpha_{\text{right}} $)
par[7]: Right tail exponent ($ n_{\text{right}} $)*/
    Double_t left = (t < -absAlphaLeft) ? 
        par[0] * TMath::Power(par[4] / absAlphaLeft, par[4]) * TMath::Exp(-0.5 * absAlphaLeft * absAlphaLeft) /
        TMath::Power(par[4] / absAlphaLeft - absAlphaLeft - t, par[4]) : 
        par[0] * TMath::Exp(-0.5 * t * t);

    Double_t right = (t > absAlphaRight) ? 
        par[5] * TMath::Power(par[7] / absAlphaRight, par[7]) * TMath::Exp(-0.5 * absAlphaRight * absAlphaRight) /
        TMath::Power(par[7] / absAlphaRight - absAlphaRight + t, par[7]) : 
        par[5] * TMath::Exp(-0.5 * t * t);

    return left + right;
}

// Triple Gaussian function
Double_t tripleGaussian(Double_t* x, Double_t* par) {
    /* par[0]: Amplitude of first Gaussian ($ A_1 $)
       par[1]: Mean of first Gaussian ($ \mu_1 $)
       par[2]: Width of first Gaussian ($ \sigma_1 $)
       par[3]: Amplitude of second Gaussian ($ A_2 $)
       par[4]: Mean of second Gaussian ($ \mu_2 $)
       par[5]: Width of second Gaussian ($ \sigma_2 $)
       par[6]: Amplitude of third Gaussian ($ A_3 $)
       par[7]: Mean of third Gaussian ($ \mu_3 $)
       par[8]: Width of third Gaussian ($ \sigma_3 $) */

    Double_t g1 = par[0] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[1]) / par[2], 2));
    Double_t g2 = par[3] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[4]) / par[5], 2));
    Double_t g3 = par[6] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[7]) / par[8], 2));

    return g1 + g2 + g3;
}

Double_t gaussianPlusPoly(Double_t *x, Double_t *par) {
    // par[0] = Amplitude (A)
    // par[1] = Mean (mu)
    // par[2] = Sigma (sigma)
    // par[3] = Quadratic coefficient (a)
    // par[4] = Linear coefficient (b)
    Double_t gaussian = par[0] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[1]) / par[2], 2));
    Double_t poly = par[3] * TMath::Power(x[0], 2) + par[4] * x[0]; // Quadratic + linear
    return gaussian + poly;
}


// Function to find the intersection point between two TF1 functions using binary search
double findIntersection(TF1* f1, TF1* f2, double xMin, double xMax, double tolerance = 1e-6) {
    double xLeft = xMin, xRight = xMax;
    double f1Val, f2Val, xMid;

    // Ensure the functions cross within the range
    f1Val = f1->Eval(xLeft) - f2->Eval(xLeft);
    if (f1Val * (f1->Eval(xRight) - f2->Eval(xRight)) >= 0) {
        cout << "Warning: No intersection found in range [" << xMin << ", " << xMax << "]. Using midpoint.\n";
        return -1.0;
    }

    // Binary search for the intersection
    while (xRight - xLeft > tolerance) {
        xMid = (xLeft + xRight) / 2.0;
        f1Val = f1->Eval(xMid) - f2->Eval(xMid);
        if (fabs(f1Val) < tolerance) {
            return xMid;
        }
        if (f1Val * (f1->Eval(xLeft) - f2->Eval(xLeft)) < 0) {
            xRight = xMid;
        } else {
            xLeft = xMid;
        }
    }
    return (xLeft + xRight) / 2.0;
}

int main() {
    // Step 1: Open the ROOT file and load the trees for pions and kaons
    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification//Skim/pkptreeCxC_9_test_modified.root");
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
    gSystem->mkdir("output/pdf/kaon_1.1", kTRUE);
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
        if (pLow < .0) {
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

    // Step 6: Set beta cuts for pions and kaons (needed for integration limits)
    vector<double> pionLeftAfter(nBins), pionRightAfter(nBins), kaonRightAfter(nBins);
    for (int i = 0; i < nBins; i++) {
        if (i == 0) {
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  kaonRightAfter[i] = 0.0;
        } else if (i == 1) {
            pionLeftAfter[i] = 0.9755;  pionRightAfter[i] = 1.0135;  kaonRightAfter[i] = 0.9721;
        } else if (i == 2) {
            pionLeftAfter[i] = 0.98;    pionRightAfter[i] = 1.0140;  kaonRightAfter[i] = 0.986;
        } else if (i == 3) {
            pionLeftAfter[i] = 0.9815;  pionRightAfter[i] = 1.0145;  kaonRightAfter[i] = 0.992;
        } else if (i == 4) {
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0141;  kaonRightAfter[i] = 0.9971;
        } else if (i == 5) { // Bin [2.5-2.8)
            pionLeftAfter[i] = 0.9848;  pionRightAfter[i] = 1.013;   kaonRightAfter[i] = 1.0;
        } else if (i == 6) {
            pionLeftAfter[i] = 0.986;   pionRightAfter[i] = 1.0128;  kaonRightAfter[i] = 1.003;
        } else if (i == 7) {
            pionLeftAfter[i] = 0.9868;  pionRightAfter[i] = 1.0121;  kaonRightAfter[i] = 1.0048;
        } else if (i == 8) {
            pionLeftAfter[i] = 0.9879;  pionRightAfter[i] = 1.0120;  kaonRightAfter[i] = 1.0049;
        } else if (i == 9) {
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  kaonRightAfter[i] = 1.012;
        } else {
            cerr << "Warning: Using default beta cuts for bin [" << pBins[i] << "-" << pBins[i+1] << ")\n";
            pionLeftAfter[i] = 0.9880;  pionRightAfter[i] = 1.0120;  kaonRightAfter[i] = 1.0;
        }
    }

    // Step 7: Create canvases for chi2pid (before) and beta (after)
    TCanvas* canvasChi2Before = new TCanvas("canvasChi2Before", "Chi2pid Fits (Before)", 1900, 900);
    canvasChi2Before->SetMargin(0.1, 0.05, 0.15, 0.1);
    TCanvas* canvasBetaAfter = new TCanvas("canvasBetaAfter", "Beta Plots (After)", 1200, 800);
    canvasBetaAfter->SetMargin(0.1, 0.05, 0.15, 0.1);

    // Open PDF files
    canvasChi2Before->Print("output/pdf/kaon_1.1/chi2pid_fits_before.pdf[");
    canvasBetaAfter->Print("output/pdf/kaon_1.1/beta_after.pdf[");

    // Step 8: Open CSV file to save chi2pid cuts
    ofstream csvFileCuts("output/csv/kaon/chi2pid_cuts.csv");
    csvFileCuts << "Momentum Bin (GeV/c),Pion Mean,Pion Sigma,Chi2 Min (Mean-5σ),Chi2 Max (Mean+5σ)\n";

    // Step 9: Open CSV file to save contamination results
    ofstream csvFileContamination("output/csv/kaon/contamination_after_chi2pid.csv");
    csvFileContamination << "Momentum Bin (GeV/c),Intersection Point,Kaon-to-Pion Ratio (%)\n";

    // Step 10: Process only bin 5
    for (int i = 4; i <= 4; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];

        // --- Before Chi2pid Cut (Chi2pid Histograms) ---
        TH1F* chi2PionsBefore = new TH1F(TString::Format("chi2_EBpions_before_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                         100, -10, 10);
        TH1F* chi2KaonsBefore = new TH1F(TString::Format("chi2_EBkaons_before_%d", i),
                                         TString::Format("p: [%.2f-%.2f) GeV/c;chi2pid;Counts", pLow, pHigh),
                                         100, -10, 10);

        // Fill chi2pid histograms
        float p, chi2pid, beta, orig_chi2pid, recomputed_chi2pid;
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

        // Fit chi2pid histograms with a Double Crystal Ball
        double pMeanBefore, pSigmaBefore, pionChi2Min, pionChi2Max, pAmpBefore;
        auto [pionAmp, pionMeanInit, pionSigmaInit] = chi2FitParams[i];

        TF1* fitPion = new TF1("gausFit", "gaus", -10, 10);
        fitPion->SetParameters(pionAmp, pionMeanInit, pionSigmaInit);
        chi2PionsBefore->Fit(fitPion, "R", "", -10, 10);
        pMeanBefore = fitPion->GetParameter(1);
        pSigmaBefore = fitPion->GetParameter(2);
        pAmpBefore = fitPion->GetParameter(0);
        delete fitPion;

        TF1* doubleCbPionBefore = new TF1("cbFit", doubleCrystalBall, -10, 10, 8);
        doubleCbPionBefore->SetParameters(
            pAmpBefore * 0.7, // A_left
            pMeanBefore,      // mu
            pSigmaBefore,     // sigma
            1.5,              // alpha_left
            2.0,              // n_left
            pAmpBefore * 0.3, // A_right
            1.5,              // alpha_right
            2.0               // n_right
        );
        doubleCbPionBefore->SetParLimits(0, 0, pAmpBefore * 1.2);
        doubleCbPionBefore->SetParLimits(1, pMeanBefore - 2 * pSigmaBefore, pMeanBefore + 2 * pSigmaBefore);
        doubleCbPionBefore->SetParLimits(2, pSigmaBefore * 0.5, pSigmaBefore * 2.0);
        doubleCbPionBefore->SetParLimits(3, 0.5, 5.0);
        doubleCbPionBefore->SetParLimits(4, 1.0, 30.0);
        doubleCbPionBefore->SetParLimits(5, 0, pAmpBefore * 0.6);
        doubleCbPionBefore->SetParLimits(6, 0.5, 5.0);
        doubleCbPionBefore->SetParLimits(7, 1.0, 30.0);
        chi2PionsBefore->Fit(doubleCbPionBefore, "R", "", -10, 10);

        pMeanBefore = doubleCbPionBefore->GetParameter(1);
        pSigmaBefore = doubleCbPionBefore->GetParameter(2);
        pAmpBefore = doubleCbPionBefore->GetParameter(0) + doubleCbPionBefore->GetParameter(5);

        pionChi2Min = pMeanBefore - 5 * pSigmaBefore;
        pionChi2Max = pMeanBefore + 5 * pSigmaBefore;

        cout << "Bin [" << pLow << "-" << pHigh << "): Pions chi2 range (before): [" 
             << pionChi2Min << ", " << pionChi2Max << "]\n";

        // Save chi2pid fit parameters to CSV
        csvFileCuts << pLow << "-" << pHigh << "," << pMeanBefore << "," << pSigmaBefore << "," 
                    << pionChi2Min << "," << pionChi2Max << "\n";

        // Plot chi2pid histograms (before cut)
        canvasChi2Before->Clear();
        canvasChi2Before->Divide(2, 1);
        canvasChi2Before->cd(1);

        // Plot chi2pid histograms (before cut)
        canvasChi2Before->Clear();
        canvasChi2Before->Divide(2, 1);
        canvasChi2Before->cd(1);     // Disable fill color to avoid bars
        chi2PionsBefore->Draw("HIST");
        doubleCbPionBefore->SetLineColor(kBlue);
        doubleCbPionBefore->SetLineWidth(1);
        doubleCbPionBefore->Draw("SAME");
        TLegend* legPions = new TLegend(0.6, 0.6, 0.9, 0.9);
        legPions->AddEntry(chi2PionsBefore, "Pions", "l");
        legPions->AddEntry(doubleCbPionBefore, "Double CB Fit", "l");
        legPions->AddEntry((TObject*)0, TString::Format("A_{left}: %.2f", doubleCbPionBefore->GetParameter(0)), "");
        legPions->AddEntry((TObject*)0, TString::Format("#mu: %.2f", doubleCbPionBefore->GetParameter(1)), "");
        legPions->AddEntry((TObject*)0, TString::Format("#sigma: %.2f", doubleCbPionBefore->GetParameter(2)), "");
        legPions->AddEntry((TObject*)0, TString::Format("#alpha_{left}: %.2f", doubleCbPionBefore->GetParameter(3)), "");
        legPions->AddEntry((TObject*)0, TString::Format("n_{left}: %.2f", doubleCbPionBefore->GetParameter(4)), "");
        legPions->AddEntry((TObject*)0, TString::Format("A_{right}: %.2f", doubleCbPionBefore->GetParameter(5)), "");
        legPions->AddEntry((TObject*)0, TString::Format("#alpha_{right}: %.2f", doubleCbPionBefore->GetParameter(6)), "");
        legPions->AddEntry((TObject*)0, TString::Format("n_{right}: %.2f", doubleCbPionBefore->GetParameter(7)), "");
        // Add chi2/NDF to the legend
    double chi2 = doubleCbPionBefore->GetChisquare();
    int ndf = doubleCbPionBefore->GetNDF();
    double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0; // Avoid division by zero
    legPions->AddEntry((TObject*)0, TString::Format("#chi^{2}/NDF: %.2f", chi2ndf), "");
        legPions->SetBorderSize(0);
        legPions->SetFillStyle(0);
        legPions->SetTextSize(0.02);
        legPions->Draw();
        
        canvasChi2Before->cd(2);
        chi2KaonsBefore->Draw();
        TLegend* legKaons = new TLegend(0.1, 0.75, 0.3, 0.9);
        legKaons->AddEntry(chi2KaonsBefore, "EB Kaons (under pion hypothesis)", "l");
        legKaons->SetBorderSize(0);
        legKaons->SetFillStyle(0);
        legKaons->SetTextSize(0.025);
        legKaons->Draw();
        canvasChi2Before->Update();
        canvasChi2Before->Print("output/pdf/kaon_1.1/chi2pid_fits_before.pdf");
        delete legPions;
        delete legKaons;

        // --- After Chi2pid Cut (±5σ) for Beta Histograms ---
        TH1F* betaPionsAfter = new TH1F(TString::Format("beta_pions_after_%d", i),
                                        TString::Format("p: [%.2f-%.2f) GeV/c; beta;Counts", pLow, pHigh),
                                        100, betaHistRanges[i].first, betaHistRanges[i].second);
        TH1F* betaKaonsAfter = new TH1F(TString::Format("beta_kaons_after_%d", i),
                                        TString::Format("p: [%.2f-%.2f) GeV/c; beta;Counts", pLow, pHigh),
                                        100, betaHistRanges[i].first, betaHistRanges[i].second);
        betaPionsAfter->Sumw2();
        betaKaonsAfter->Sumw2();

        // Fill beta histograms (after cut) with ±5σ cut
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("beta", &beta);
        treePions->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);
        for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= -3 && recomputed_chi2pid <= 3) {
                betaPionsAfter->Fill(beta);
            }
        }
        treeKaons->SetBranchAddress("p", &p);
        treeKaons->SetBranchAddress("beta", &beta);
        treeKaons->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);
        for (Long64_t j = 0; j < treeKaons->GetEntries(); j++) {
            treeKaons->GetEntry(j);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= -3 && recomputed_chi2pid <= 3) {
                betaKaonsAfter->Fill(beta);
            }
        }

        // Fit beta histograms (after cut)
        double pMeanAfter = 0, pSigmaAfter = 0, pConstAfter = 0, kMeanAfter = 0, kSigmaAfter = 0, kConstAfter = 0;
        TF1* cbPionBetaAfter = nullptr;
        TF1* gausPion = nullptr;
        TF1* gausKaonAfter = nullptr;

       /*  if (betaPionsAfter->GetEntries() >= 10) {
            // Estimate peak height, mean, and width directly from histogram
            double peak = betaPionsAfter->GetMaximum();
            int binMax = betaPionsAfter->GetMaximumBin();
            double meanGuess = betaPionsAfter->GetBinCenter(binMax);
            double sigmaGuess = betaPionsAfter->GetRMS() * 0.5;  // Can adjust multiplier
        
            // Create Double Crystal Ball function
            cbPionBetaAfter = new TF1("cbPionBetaAfter", doubleCrystalBall, 0.93, 1.03, 8); // Extended range to 0.93
        
            // Initial parameters: A_left, mu, sigma, alpha_left, n_left, A_right, alpha_right, n_right
            cbPionBetaAfter->SetParameters(
                peak * 0.9,     // A_left
                meanGuess,      // mu
                sigmaGuess,     // sigma
                1,            // alpha_left
                2.0,            // n_left
                peak * 0.1,     // A_right
                0.9,            // alpha_right
                10.0            // n_right
            );
        
            // Parameter limits to constrain the fit
            cbPionBetaAfter->SetParLimits(0, 0.0, peak * 1.2); // A_left
            cbPionBetaAfter->SetParLimits(1, meanGuess - 0.01, meanGuess + 0.01); // mu
            cbPionBetaAfter->SetParLimits(2, sigmaGuess * 0.5, sigmaGuess * 2.0); // sigma
            cbPionBetaAfter->SetParLimits(3, 0.8, 1.2);   // alpha_left
            cbPionBetaAfter->SetParLimits(4, 1.5, 3.0);   // n_left
            cbPionBetaAfter->SetParLimits(5, 0.0, peak);  // A_right
            cbPionBetaAfter->SetParLimits(6, 0.5, 1.5);   // alpha_right
            cbPionBetaAfter->SetParLimits(7, 5.0, 15.0);  // n_right
        
            // Fit the histogram with CB function
            betaPionsAfter->Fit(cbPionBetaAfter, "R", "", 0.93, 1.03);
        
            // Extract final fit parameters
            pMeanAfter = cbPionBetaAfter->GetParameter(1);
            pSigmaAfter = cbPionBetaAfter->GetParameter(2);
            pConstAfter = cbPionBetaAfter->GetParameter(0) + cbPionBetaAfter->GetParameter(5); // A_left + A_right
        }  */

        /* if (betaPionsAfter->GetEntries() >= 10) {
            // Estimate peak height, mean, and width directly from histogram
            double peak = betaPionsAfter->GetMaximum();
            int binMax = betaPionsAfter->GetMaximumBin();
            double meanGuess = betaPionsAfter->GetBinCenter(binMax);
            double sigmaGuess = betaPionsAfter->GetRMS() * 0.5;
        
            // Create Triple Gaussian function
            cbPionBetaAfter = new TF1("cbPionBetaAfter", tripleGaussian, 0.93, 1.03, 9);
        
            // Initial parameters: A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3
            cbPionBetaAfter->SetParameters(
                peak * 0.8,     // A1 (main peak amplitude)
                meanGuess,      // mu1 (mean of main peak)
                sigmaGuess,     // sigma1 (width of main peak)
                peak * 0.1,     // A2 (left tail amplitude)
                meanGuess - 0.01, // mu2 (left tail mean, slightly left of peak)
                sigmaGuess * 2.0, // sigma2 (wider for left tail)
                peak * 0.1,     // A3 (right tail amplitude)
                meanGuess + 0.01, // mu3 (right tail mean, slightly right of peak)
                sigmaGuess * 1.5  // sigma3 (wider for right tail)
            );
        
            // Parameter limits to constrain the fit
            cbPionBetaAfter->SetParLimits(0, 0.0, peak * 1.2); // A1
            cbPionBetaAfter->SetParLimits(1, meanGuess - 0.02, meanGuess + 0.02); // mu1
            cbPionBetaAfter->SetParLimits(2, sigmaGuess * 0.5, sigmaGuess * 2.0); // sigma1
            cbPionBetaAfter->SetParLimits(3, 0.0, peak * 0.5); // A2
            cbPionBetaAfter->SetParLimits(4, 0.93, meanGuess - 0.005); // mu2 (left of peak)
            cbPionBetaAfter->SetParLimits(5, sigmaGuess, sigmaGuess * 3.0); // sigma2
            cbPionBetaAfter->SetParLimits(6, 0.0, peak * 0.5); // A3
            cbPionBetaAfter->SetParLimits(7, meanGuess + 0.005, 1.03); // mu3 (right of peak)
            cbPionBetaAfter->SetParLimits(8, sigmaGuess, sigmaGuess * 3.0); // sigma3
        
            // Fit the histogram with triple Gaussian
            betaPionsAfter->Fit(cbPionBetaAfter, "R", "", 0.986, 1.03);
        
            // Extract final fit parameters
            pMeanAfter = (cbPionBetaAfter->GetParameter(1) * cbPionBetaAfter->GetParameter(0) +
                          cbPionBetaAfter->GetParameter(4) * cbPionBetaAfter->GetParameter(3) +
                          cbPionBetaAfter->GetParameter(7) * cbPionBetaAfter->GetParameter(6)) /
                         (cbPionBetaAfter->GetParameter(0) + cbPionBetaAfter->GetParameter(3) + cbPionBetaAfter->GetParameter(6));
            pSigmaAfter = TMath::Sqrt(
                (cbPionBetaAfter->GetParameter(2) * cbPionBetaAfter->GetParameter(2) * cbPionBetaAfter->GetParameter(0) +
                 cbPionBetaAfter->GetParameter(5) * cbPionBetaAfter->GetParameter(5) * cbPionBetaAfter->GetParameter(3) +
                 cbPionBetaAfter->GetParameter(8) * cbPionBetaAfter->GetParameter(8) * cbPionBetaAfter->GetParameter(6)) /
                (cbPionBetaAfter->GetParameter(0) + cbPionBetaAfter->GetParameter(3) + cbPionBetaAfter->GetParameter(6)));
            pConstAfter = cbPionBetaAfter->GetParameter(0) + cbPionBetaAfter->GetParameter(3) + cbPionBetaAfter->GetParameter(6);
        
            // Debug: Check fit status
            int fitStatus = betaPionsAfter->GetFunction("cbPionBetaAfter")->GetNumberFitPoints();
            if (fitStatus <= 0) {
                cout << "Fit failed or no points used. Check parameters and range." << endl;
            }
        }   */ 

       
        if (betaPionsAfter->GetEntries() >= 10) {
            // Estimate peak height, mean, and width directly from histogram
            double peak = betaPionsAfter->GetMaximum();
            int binMax = betaPionsAfter->GetMaximumBin();
            double meanGuess = betaPionsAfter->GetBinCenter(binMax);
            double sigmaGuess = betaPionsAfter->GetRMS() * 0.5;
            
            // Create Triple Gaussian function
            cbPionBetaAfter = new TF1("cbPionBetaAfter", tripleGaussian, 0.986, 1.03, 9);
            
            // Initial parameters: A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3
            cbPionBetaAfter->SetParameters(
            peak * 0.8,     // A1 (main peak amplitude)
            meanGuess,      // mu1 (mean of main peak)
            sigmaGuess,     // sigma1 (width of main peak)
            peak * 0.05,    // A2 (left tail amplitude, reduced to focus on pion tail)
            0.99,           // mu2 (left tail mean, set to ~0.99 near the peak)
            0.002,          // sigma2 (left tail width, narrowed to 0.002)
            peak * 0.1,     // A3 (right tail amplitude, keep as is since right tail is good)
            meanGuess + 0.01, // mu3 (right tail mean, slightly right of peak)
            sigmaGuess * 1.5  // sigma3 (wider for right tail)
            );
            
            // Parameter limits to constrain the fit
            cbPionBetaAfter->SetParLimits(0, 0.0, peak * 1.2); // A1
            cbPionBetaAfter->SetParLimits(1, meanGuess - 0.01, meanGuess + 0.01); // mu1
            cbPionBetaAfter->SetParLimits(2, sigmaGuess * 0.5, sigmaGuess * 2.0); // sigma1
            cbPionBetaAfter->SetParLimits(3, 0.0, peak * 0.1); // A2 (reduced range for pion tail)
            cbPionBetaAfter->SetParLimits(4, 0.986, 0.995); // mu2 (tighter range near 0.99)
            cbPionBetaAfter->SetParLimits(5, 0.001, 0.003); // sigma2 (narrower range)
            cbPionBetaAfter->SetParLimits(6, 0.0, peak * 0.5); // A3
            cbPionBetaAfter->SetParLimits(7, meanGuess + 0.005, 1.03); // mu3
            cbPionBetaAfter->SetParLimits(8, sigmaGuess, sigmaGuess * 3.0); // sigma3
            
            // Fit the histogram with triple Gaussian
            betaPionsAfter->Fit(cbPionBetaAfter, "R", "", 0.986, 1.03);
            
            // Extract final fit parameters
            pMeanAfter = (cbPionBetaAfter->GetParameter(1) * cbPionBetaAfter->GetParameter(0) +
            cbPionBetaAfter->GetParameter(4) * cbPionBetaAfter->GetParameter(3) +
            cbPionBetaAfter->GetParameter(7) * cbPionBetaAfter->GetParameter(6)) /
            (cbPionBetaAfter->GetParameter(0) + cbPionBetaAfter->GetParameter(3) + cbPionBetaAfter->GetParameter(6));
            pSigmaAfter = TMath::Sqrt(
            (cbPionBetaAfter->GetParameter(2) * cbPionBetaAfter->GetParameter(2) * cbPionBetaAfter->GetParameter(0) +
            cbPionBetaAfter->GetParameter(5) * cbPionBetaAfter->GetParameter(5) * cbPionBetaAfter->GetParameter(3) +
            cbPionBetaAfter->GetParameter(8) * cbPionBetaAfter->GetParameter(8) * cbPionBetaAfter->GetParameter(6)) /
            (cbPionBetaAfter->GetParameter(0) + cbPionBetaAfter->GetParameter(3) + cbPionBetaAfter->GetParameter(6)));
            pConstAfter = cbPionBetaAfter->GetParameter(0) + cbPionBetaAfter->GetParameter(3) + cbPionBetaAfter->GetParameter(6);
            
            // Debug: Check fit status
            int fitStatus = betaPionsAfter->GetFunction("cbPionBetaAfter")->GetNumberFitPoints();
            if (fitStatus <= 0) {
            cout << "Fit failed or no points used. Check parameters and range." << endl;
            }
            } 

            
        if (betaKaonsAfter->GetEntries() >= 10) {
            gausKaonAfter = new TF1(TString::Format("gaus_%s", betaKaonsAfter->GetName()), "gaus", betaFitRanges[i].first, betaFitRanges[i].second);
            gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.987, 0.01);
            betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.982, 0.993); // Adjusted range for bin 5
            kMeanAfter = gausKaonAfter->GetParameter(1);
            kSigmaAfter = gausKaonAfter->GetParameter(2);
            kConstAfter = gausKaonAfter->GetParameter(0);
        }

        // Calculate contamination (after finding intersection)
        double contaminationAfter = -1.0;
        double intersectionPoint = -1.0;
        if (betaPionsAfter->GetEntries() >= 10 && betaKaonsAfter->GetEntries() >= 10) {
            // Find the intersection point between kaon Gaussian and pion double CB
            intersectionPoint = findIntersection(gausKaonAfter, cbPionBetaAfter, kMeanAfter, pMeanAfter);

            // Integrate kaon Gaussian from intersection to the right
            double kaonIntegral = gausKaonAfter->Integral(intersectionPoint, betaHistRanges[i].second);

            // Integrate pion double CB from left to right
            double pionIntegral = cbPionBetaAfter->Integral(pionLeftAfter[i], pionRightAfter[i]);

            // Compute contamination ratio
            if (pionIntegral > 0) {
                contaminationAfter = (kaonIntegral / pionIntegral) * 100.0;
            } else {
                cout << "Error: Pion integral is zero or negative for bin " << i << "\n";
                contaminationAfter = -1.0;
            }

            cout << "Bin [" << pLow << "-" << pHigh << "): Intersection at beta = " << intersectionPoint
                 << ", Kaon Integral = " << kaonIntegral << ", Pion Integral = " << pionIntegral
                 << ", Contamination = " << contaminationAfter << "%\n";

            // Save to CSV
            csvFileContamination << pLow << "-" << pHigh << "," << intersectionPoint << "," << contaminationAfter << "\n";
        } else {
            cout << "Bin [" << pLow << "-" << pHigh << "): Not enough entries to calculate contamination\n";
            csvFileContamination << pLow << "-" << pHigh << ",N/A,N/A\n";
        }

        // Plot beta histograms (after cut)
        double yMin = 0.0;
        double yMax = -1.0;
        canvasBetaAfter->Clear();
        if (betaPionsAfter->GetEntries() >= 10) {
            betaPionsAfter->SetMarkerStyle(20);
            betaPionsAfter->SetMarkerSize(1);
            betaPionsAfter->SetMarkerColor(kBlue);
            betaPionsAfter->Draw("P");
            yMax = betaPionsAfter->GetMaximum();
            if (auto* fit = betaPionsAfter->GetFunction("cbPionBetaAfter")) {
                fit->SetLineColor(kBlue);
                fit->SetLineStyle(kDashed);
                fit->Draw("SAME");
                TLegend* leg = new TLegend(0.6, 0.6, 0.9, 0.9);
                leg->AddEntry(betaPionsAfter, "Pions", "p");
                leg->AddEntry(fit, "Triple Gaussian Fit", "l");
                leg->AddEntry((TObject*)0, TString::Format("A1: %.5f", fit->GetParameter(0)), "");
                leg->AddEntry((TObject*)0, TString::Format("#mu1: %.5f", fit->GetParameter(1)), "");
                leg->AddEntry((TObject*)0, TString::Format("#sigma1: %.5f", fit->GetParameter(2)), "");
                leg->AddEntry((TObject*)0, TString::Format("A2: %.5f", fit->GetParameter(3)), "");
                leg->AddEntry((TObject*)0, TString::Format("#mu2: %.5f", fit->GetParameter(4)), "");
                leg->AddEntry((TObject*)0, TString::Format("#sigma2: %.5f", fit->GetParameter(5)), "");
                leg->AddEntry((TObject*)0, TString::Format("A3: %.5f", fit->GetParameter(6)), "");
                leg->AddEntry((TObject*)0, TString::Format("#mu3: %.5f", fit->GetParameter(7)), "");
                leg->AddEntry((TObject*)0, TString::Format("#sigma3: %.5f", fit->GetParameter(8)), "");
                leg->AddEntry((TObject*)0, TString::Format("chi2/NDF: %.2f", fit->GetChisquare() / fit->GetNDF()), "");
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextSize(0.02);
                leg->Draw();
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
                fit->SetRange(0.95, 1.03);
                fit->Draw("SAME");
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
        canvasBetaAfter->Print("output/pdf/kaon_1.1/beta_after.pdf");
        delete legAfter;

        // Clean up histograms and functions
        delete chi2PionsBefore;
        delete chi2KaonsBefore;
        delete betaPionsAfter;
        delete betaKaonsAfter;
        if (cbPionBetaAfter) delete cbPionBetaAfter;
        if (gausPion) delete gausPion;
        if (gausKaonAfter) delete gausKaonAfter;
    }

    // Close PDF and CSV files
    canvasChi2Before->Print("output/pdf/kaon_1.1/chi2pid_fits_before.pdf]");
    canvasBetaAfter->Print("output/pdf/kaon_1.1/beta_after.pdf]");
    csvFileCuts.close();
    csvFileContamination.close();

    // Final cleanup
    delete canvasChi2Before;
    delete canvasBetaAfter;
    file->Close();
    delete file;

    return 0;
}