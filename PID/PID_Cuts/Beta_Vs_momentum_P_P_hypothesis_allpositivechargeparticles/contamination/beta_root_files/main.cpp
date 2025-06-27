#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TLine.h>
#include <TF1.h>
#include <iostream>
#include <vector>

using namespace std;
using namespace ROOT;

// Calculate theoretical beta
double getTheoreticalBeta(double p, double mass) {
    return p / sqrt(p * p + mass * mass);
}

// Single Gaussian
Double_t singleGaussian(Double_t* x, Double_t* par) {
    return par[0] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[1]) / par[2], 2)) / (par[2] * TMath::Sqrt(2 * TMath::Pi()));
}

// Two Gaussians (pion + kaon)
Double_t twoGaussians(Double_t* x, Double_t* par) {
    Double_t pi = par[0] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[1]) / par[2], 2)) / (par[2] * TMath::Sqrt(2 * TMath::Pi()));
    Double_t ka = par[3] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[4]) / par[5], 2)) / (par[5] * TMath::Sqrt(2 * TMath::Pi()));
    return pi + ka;
}

// Three Gaussians (pion + kaon + proton)
Double_t threeGaussians(Double_t* x, Double_t* par) {
    Double_t pi = par[0] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[1]) / par[2], 2)) / (par[2] * TMath::Sqrt(2 * TMath::Pi()));
    Double_t ka = par[3] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[4]) / par[5], 2)) / (par[5] * TMath::Sqrt(2 * TMath::Pi()));
    Double_t pr = par[6] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[7]) / par[8], 2)) / (par[8] * TMath::Sqrt(2 * TMath::Pi()));
    return pi + ka + pr;
}

int main() {
    TFile* infile = new TFile("beta_histograms.root", "READ");
    if (!infile || infile->IsZombie()) {
        cerr << "Error: Cannot open file beta_histograms.root\n";
        return 1;
    }
    gSystem->mkdir("output", kTRUE);

    vector<double> pBins;
    double pMin = 1.0, pMax = 7.0;
    int nBins = static_cast<int>((pMax - pMin) / 0.3);
    for (int i = 0; i <= nBins; i++) {
        pBins.push_back(pMin + i * 0.3);
    }

    TCanvas* canvas = new TCanvas("canvas", "Beta Distribution Fits", 1200, 800);
    canvas->Print("output/beta_fits.pdf[");

    const double pionMass = 0.1396;
    const double kaonMass = 0.4937;
    const double protonMass = 0.9383;
    const double sigmaBetaInitial = 0.0114;

    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];
        double pMid = (pLow + pHigh) / 2.0;

        TH1D* histo_before = (TH1D*)infile->Get(TString::Format("beta_before_%d", i));
        TH1D* histo_after = (TH1D*)infile->Get(TString::Format("beta_after_%d", i));
        if (!histo_before || !histo_after) {
            cerr << "Error: Histogram not found for bin " << i << endl;
            continue;
        }

        double betaPion = getTheoreticalBeta(pMid, pionMass);
        double betaKaon = getTheoreticalBeta(pMid, kaonMass);
        double betaProton = getTheoreticalBeta(pMid, protonMass);

        canvas->Clear();
        histo_before->Draw("P");
        histo_after->Draw("P SAME");

        TF1* fit;
        double initialAmp = histo_before->GetMaximum() * 0.8; // Base amplitude

        if (pMid <= 1.6) {
            // Region 1: Pion only
            fit = new TF1("fit", singleGaussian, 0.94, 1.03, 3);
            fit->SetParameters(initialAmp, betaPion, sigmaBetaInitial);
            fit->SetParLimits(1, betaPion - 0.05, betaPion + 0.05); // Tune mu range if needed
            fit->SetParLimits(2, 0.005, 0.03); // Tune sigma range if needed
            fit->SetParLimits(0, 0, 2 * initialAmp);
        } else if (pMid <= 2.8) {
            // Region 2: Pion + Kaon
            fit = new TF1("fit", twoGaussians, 0.94, 1.03, 6);
            fit->SetParameters(initialAmp, betaPion, sigmaBetaInitial, initialAmp * 0.2, betaKaon, sigmaBetaInitial);
            fit->SetParLimits(1, betaPion - 0.05, betaPion + 0.05); // Adjust if pion peak shifts
            fit->SetParLimits(4, betaKaon - 0.05, betaKaon + 0.05); // Adjust for kaon position
            fit->SetParLimits(2, 0.005, 0.03); // Shared sigma, tune if tails are wide
            fit->SetParLimits(5, 0.005, 0.03);
            fit->SetParLimits(0, 0, 2 * initialAmp);
            fit->SetParLimits(3, 0, 2 * initialAmp * 0.2); // Tune kaon amplitude ratio
        } else {
            // Region 3: Pion + Kaon + Proton
            fit = new TF1("fit", threeGaussians, 0.94, 1.03, 9);
            fit->SetParameters(initialAmp, betaPion, sigmaBetaInitial, initialAmp * 0.2, betaKaon, sigmaBetaInitial, initialAmp * 0.1, betaProton, sigmaBetaInitial);
            fit->SetParLimits(1, betaPion - 0.05, betaPion + 0.05);
            fit->SetParLimits(4, betaKaon - 0.05, betaKaon + 0.05);
            fit->SetParLimits(7, betaProton - 0.05, betaProton + 0.05); // Tune proton range
            fit->SetParLimits(2, 0.005, 0.03);
            fit->SetParLimits(5, 0.005, 0.03);
            fit->SetParLimits(8, 0.005, 0.03);
            fit->SetParLimits(0, 0, 2 * initialAmp);
            fit->SetParLimits(3, 0, 2 * initialAmp * 0.2);
            fit->SetParLimits(6, 0, 2 * initialAmp * 0.1); // Tune proton amplitude
        }

        histo_before->Fit(fit, "R0");

        // Extract fit parameters
        double N_pi = fit->GetParameter(0);
        double N_K = (pMid <= 1.6) ? 0 : fit->GetParameter(3);
        double N_p = (pMid <= 2.8) ? 0 : fit->GetParameter(6);
        double chi2 = fit->GetChisquare();
        int ndf = fit->GetNDF();

        // Draw fit and theoretical lines
        fit->SetLineColor(kMagenta);
        fit->Draw("SAME");
        double maxHeight = max(histo_before->GetMaximum(), histo_after->GetMaximum());
        TLine* linePion = new TLine(betaPion, 0, betaPion, maxHeight);
        linePion->SetLineColor(kBlue);
        linePion->SetLineStyle(2);
        linePion->Draw();
        if (pMid > 1.6) {
            TLine* lineKaon = new TLine(betaKaon, 0, betaKaon, maxHeight);
            lineKaon->SetLineColor(kGreen);
            lineKaon->SetLineStyle(2);
            lineKaon->Draw();
        }
        if (pMid > 2.8) {
            TLine* lineProton = new TLine(betaProton, 0, betaProton, maxHeight);
            lineProton->SetLineColor(kRed);
            lineProton->SetLineStyle(2);
            lineProton->Draw();
        }

        // Add legend
        TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg->AddEntry(histo_before, "Before Cut", "p");
        leg->AddEntry(histo_after, "After Cut", "p");
        leg->AddEntry(fit, "Fit", "l");
        leg->AddEntry(linePion, "Pion", "l");
        if (pMid > 1.6) leg->AddEntry((TObject*)0x1, "Kaon", "l");
        if (pMid > 2.8) leg->AddEntry((TObject*)0x1, "Proton", "l");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw();

        cout << "Momentum range [" << pLow << ", " << pHigh << ") GeV/c:" << endl;
        cout << "Fitted Yields - N_pi: " << N_pi << ", N_K: " << N_K << ", N_p: " << N_p << endl;
        cout << "Chi2/NDf: " << chi2 << "/" << ndf << " = " << (ndf > 0 ? chi2 / ndf : 0) << endl;

        canvas->Print(TString::Format("output/beta_fits.pdf"));
        delete leg;
        delete linePion;
       
        delete fit;
    }

    canvas->Print("output/beta_fits.pdf]");
    delete canvas;
    infile->Close();
    delete infile;

    return 0;
}