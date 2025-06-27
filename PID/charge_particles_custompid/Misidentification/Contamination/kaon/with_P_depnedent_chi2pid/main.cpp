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
    Double_t g1 = par[0] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[1]) / par[2], 2));
    Double_t g2 = par[3] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[4]) / par[5], 2));
    Double_t g3 = par[6] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[7]) / par[8], 2));
    return g1 + g2 + g3;
}

// Function to compute chi2pid cut
double getChi2Cut(double p) {
    const double C = 0.9795; // Scaling factor
    if (p <= 2.75) {
        return 3.0 * C; // 3.0 * 0.9795 ≈ 2.9385
    } else {
        return C * (-0.20 + 3.00 * exp(-p / 4.51) + 37.26 * exp(-p / 0.87));
    }
}

// Function to find the intersection point between two TF1 functions using binary search
double findIntersection(TF1* f1, TF1* f2, double xMin, double xMax, double tolerance = 1e-6) {
    double xLeft = xMin, xRight = xMax;
    double f1Val, f2Val, xMid;
    f1Val = f1->Eval(xLeft) - f2->Eval(xLeft);
    if (f1Val * (f1->Eval(xRight) - f2->Eval(xRight)) >= 0) {
        cout << "Warning: No intersection found in range [" << xMin << ", " << xMax << "]. Using midpoint.\n";
        return (xLeft + xRight) / 2.0;
    }
    while (xRight - xLeft > tolerance) {
        xMid = (xLeft + xRight) / 2.0;
        f1Val = f1->Eval(xMid) - f2->Eval(xMid);
        if (fabs(f1Val) < tolerance) return xMid;
        if (f1Val * (f1->Eval(xLeft) - f2->Eval(xLeft)) < 0) xRight = xMid;
        else xLeft = xMid;
    }
    return (xLeft + xRight) / 2.0;
}

int main() {
    // Step 1: Open the ROOT file and load the trees
    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification//Skim/pkptreeCxC_9_test_modified.root");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file pkptreeCxC_9_test_modified.root\n";
        return 1;
    }
    TTree* treePions = (TTree*)file->Get("EB_pid_pions");
    TTree* treeKaons = (TTree*)file->Get("EB_pid_kaons");
    if (!treePions || !treeKaons) {
        cerr << "Error: Cannot load trees\n";
        file->Close();
        delete file;
        return 1;
    }

    // Step 2: Create output directories
    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/pdf/kaon_1.1", kTRUE);
    gSystem->mkdir("output/csv/kaon", kTRUE);

    // Step 3: Set up momentum bins (1 to 7 GeV, 20 bins)
    vector<double> pBins;
    double pMin = 1.0, pMax = 7.0;
    int nBins = 20;
    double binWidth = (pMax - pMin) / nBins;
    for (int i = 0; i <= nBins; i++) pBins.push_back(pMin + i * binWidth);

    // Step 4: Set beta histogram ranges
    vector<pair<double, double>> betaHistRanges(nBins);
    for (int i = 0; i < nBins; i++) betaHistRanges[i] = {0.95, 1.03};

    // Step 5: Set beta cuts for pions and kaons
    vector<double> pionLeftAfter(nBins), pionRightAfter(nBins), kaonRightAfter(nBins);
    for (int i = 0; i < nBins; i++) {
        if (i == 5) {
            pionLeftAfter[i] = 0.9848; pionRightAfter[i] = 1.013; kaonRightAfter[i] = 1.0;
        } else {
            pionLeftAfter[i] = 0.9880; pionRightAfter[i] = 1.0120; kaonRightAfter[i] = 1.0;
        }
    }

    // Step 6: Create canvas and open PDF file
    TCanvas* canvasBetaAfter = new TCanvas("canvasBetaAfter", "Beta Plots (After)", 1200, 800);
    canvasBetaAfter->SetMargin(0.1, 0.05, 0.15, 0.1);
    canvasBetaAfter->Print("output/pdf/kaon_1.1/beta_after.pdf[");

    // Step 7: Open CSV file for contamination results
    ofstream csvFileContamination("output/csv/kaon/contamination_after_chi2pid.csv");
    csvFileContamination << "Momentum Bin (GeV/c),Intersection Point,Kaon-to-Pion Ratio (%)\n";

    // Step 8: Process all bins
    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];

        // Create beta histograms after chi2pid cut
        TH1F* betaPionsAfter = new TH1F(TString::Format("beta_pions_after_%d", i),
                                        TString::Format("p: [%.2f-%.2f) GeV/c; beta;Counts", pLow, pHigh),
                                        100, betaHistRanges[i].first, betaHistRanges[i].second);
        TH1F* betaKaonsAfter = new TH1F(TString::Format("beta_kaons_after_%d", i),
                                        TString::Format("p: [%.2f-%.2f) GeV/c; beta;Counts", pLow, pHigh),
                                        100, betaHistRanges[i].first, betaHistRanges[i].second);
        betaPionsAfter->Sumw2();
        betaKaonsAfter->Sumw2();

        // Fill beta histograms with momentum-dependent chi2pid cut
        float p, beta, recomputed_chi2pid;
        const double C = 0.9795;
        const double chi2NegCut = -3.0 * C; // Negative cut for pion selection
        treePions->SetBranchAddress("p", &p);
        treePions->SetBranchAddress("beta", &beta);
        treePions->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);
        for (Long64_t j = 0; j < treePions->GetEntries(); j++) {
            treePions->GetEntry(j);
            double chi2PosCut = getChi2Cut(p);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= chi2NegCut && recomputed_chi2pid <= chi2PosCut) {
                betaPionsAfter->Fill(beta);
            }
        }
        treeKaons->SetBranchAddress("p", &p);
        treeKaons->SetBranchAddress("beta", &beta);
        treeKaons->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);
        for (Long64_t j = 0; j < treeKaons->GetEntries(); j++) {
            treeKaons->GetEntry(j);
            double chi2PosCut = getChi2Cut(p);
            if (p >= pLow && p < pHigh && recomputed_chi2pid >= chi2NegCut && recomputed_chi2pid <= chi2PosCut) {
                betaKaonsAfter->Fill(beta);
            }
        }

        // Fit beta histograms
        TF1* cbPionBetaAfter = nullptr;
        TF1* gausKaonAfter = nullptr;
        double pMeanAfter = 0, pConstAfter = 0, kMeanAfter = 0, kConstAfter = 0;
        if (betaPionsAfter->GetEntries() >= 10) {
            double peak = betaPionsAfter->GetMaximum();
            int binMax = betaPionsAfter->GetMaximumBin();
            double meanGuess = betaPionsAfter->GetBinCenter(binMax);
            double sigmaGuess = betaPionsAfter->GetRMS() * 0.5;
            cbPionBetaAfter = new TF1("cbPionBetaAfter", tripleGaussian, 0.986, 1.03, 9);
            cbPionBetaAfter->SetParameters(peak * 0.8, meanGuess, sigmaGuess, peak * 0.05, 0.99, 0.002, peak * 0.1, meanGuess + 0.01, sigmaGuess * 1.5);
            cbPionBetaAfter->SetParLimits(0, 0.0, peak * 1.2);
            cbPionBetaAfter->SetParLimits(1, meanGuess - 0.01, meanGuess + 0.01);
            cbPionBetaAfter->SetParLimits(2, sigmaGuess * 0.5, sigmaGuess * 2.0);
            cbPionBetaAfter->SetParLimits(3, 0.0, peak * 0.1);
            cbPionBetaAfter->SetParLimits(4, 0.986, 0.995);
            cbPionBetaAfter->SetParLimits(5, 0.001, 0.003);
            cbPionBetaAfter->SetParLimits(6, 0.0, peak * 0.5);
            cbPionBetaAfter->SetParLimits(7, meanGuess + 0.005, 1.03);
            cbPionBetaAfter->SetParLimits(8, sigmaGuess, sigmaGuess * 3.0);
            betaPionsAfter->Fit(cbPionBetaAfter, "R", "", 0.986, 1.03);
            pMeanAfter = (cbPionBetaAfter->GetParameter(1) * cbPionBetaAfter->GetParameter(0) +
                          cbPionBetaAfter->GetParameter(4) * cbPionBetaAfter->GetParameter(3) +
                          cbPionBetaAfter->GetParameter(7) * cbPionBetaAfter->GetParameter(6)) /
                         (cbPionBetaAfter->GetParameter(0) + cbPionBetaAfter->GetParameter(3) + cbPionBetaAfter->GetParameter(6));
            pConstAfter = cbPionBetaAfter->GetParameter(0) + cbPionBetaAfter->GetParameter(3) + cbPionBetaAfter->GetParameter(6);
        }
        if (betaKaonsAfter->GetEntries() >= 10) {
            gausKaonAfter = new TF1(TString::Format("gaus_%s", betaKaonsAfter->GetName()), "gaus", 0.95, 1.03);
            gausKaonAfter->SetParameters(betaKaonsAfter->GetMaximum(), 0.987, 0.01);
            betaKaonsAfter->Fit(gausKaonAfter, "R", "", 0.982, 0.993);
            kMeanAfter = gausKaonAfter->GetParameter(1);
        }

        // Calculate contamination using intersection
        double contaminationAfter = -1.0;
        double intersectionPoint = -1.0;
        if (betaPionsAfter->GetEntries() >= 10 && betaKaonsAfter->GetEntries() >= 10) {
            // Find the intersection point between kaon Gaussian and pion triple Gaussian
            intersectionPoint = findIntersection(gausKaonAfter, cbPionBetaAfter, kMeanAfter, pMeanAfter);

            // Integrate kaon Gaussian from intersection to the right
            double kaonIntegral = gausKaonAfter->Integral(intersectionPoint, betaHistRanges[i].second);

            // Integrate pion triple Gaussian from left to right
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

        // Plot beta histograms
        double yMax = -1.0;
        canvasBetaAfter->Clear();
        if (betaPionsAfter->GetEntries() >= 10) {
            betaPionsAfter->SetMarkerStyle(20);
            betaPionsAfter->SetMarkerSize(1);
            betaPionsAfter->SetMarkerColor(kBlue);
            betaPionsAfter->Draw("P");
            yMax = betaPionsAfter->GetMaximum();
            if (cbPionBetaAfter) {
                cbPionBetaAfter->SetLineColor(kBlue);
                cbPionBetaAfter->SetLineStyle(kDashed);
                cbPionBetaAfter->Draw("SAME");
            }
        }
        if (betaKaonsAfter->GetEntries() >= 10) {
            betaKaonsAfter->SetMarkerStyle(20);
            betaKaonsAfter->SetMarkerSize(1);
            betaKaonsAfter->SetMarkerColor(kGreen);
            betaKaonsAfter->Draw(betaPionsAfter->GetEntries() >= 10 ? "P SAME" : "P");
            yMax = max(yMax, betaKaonsAfter->GetMaximum());
            if (gausKaonAfter) {
                gausKaonAfter->SetLineColor(kGreen);
                gausKaonAfter->SetLineStyle(kDashed);
                gausKaonAfter->Draw("SAME");
            }
        }
        yMax *= 1.2;
        if (betaPionsAfter->GetEntries() >= 10 || betaKaonsAfter->GetEntries() >= 10) {
            (betaPionsAfter->GetEntries() >= 10 ? betaPionsAfter : betaKaonsAfter)->GetYaxis()->SetRangeUser(0, yMax);
        }
        TLegend* legAfter = new TLegend(0.15, 0.75, 0.35, 0.88);
        if (betaPionsAfter->GetEntries() >= 10) legAfter->AddEntry(betaPionsAfter, "Pions", "p");
        if (betaKaonsAfter->GetEntries() >= 10) legAfter->AddEntry(betaKaonsAfter, "Kaons (Misidentified)", "p");
        legAfter->SetTextSize(0.03);
        legAfter->Draw();
        canvasBetaAfter->SetGrid();
        canvasBetaAfter->Update();
        canvasBetaAfter->Print("output/pdf/kaon_1.1/beta_after.pdf");

        // Clean up histograms and functions for this bin
        delete betaPionsAfter;
        delete betaKaonsAfter;
        if (cbPionBetaAfter) delete cbPionBetaAfter;
        if (gausKaonAfter) delete gausKaonAfter;
        delete legAfter;
    }

    // Finalize PDF and CSV
    canvasBetaAfter->Print("output/pdf/kaon_1.1/beta_after.pdf]");
    csvFileContamination.close();

    // Clean up
    delete canvasBetaAfter;
    file->Close();
    delete file;

    return 0;
}