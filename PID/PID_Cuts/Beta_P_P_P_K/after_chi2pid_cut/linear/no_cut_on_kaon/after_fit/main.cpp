#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TF1.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <TMath.h>

using namespace std;
using namespace ROOT;

// Triple Gaussian function
Double_t tripleGaussian(Double_t* x, Double_t* par) {
    Double_t g1 = par[0] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[1]) / par[2], 2));
    Double_t g2 = par[3] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[4]) / par[5], 2));
    Double_t g3 = par[6] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[7]) / par[8], 2));
    return g1 + g2 + g3;
}

// Gaussian function
Double_t gaussian(Double_t* x, Double_t* par) {
    return par[0] * TMath::Exp(-0.5 * TMath::Power((x[0] - par[1]) / par[2], 2));
}

// Function to find the intersection point between two TF1 functions
double findIntersection(TF1* f1, TF1* f2, double xMin, double xMax, double tolerance = 1e-6) {
    double xLeft = xMin, xRight = xMax;
    double f1Val, f2Val, xMid;
    f1Val = f1->Eval(xLeft) - f2->Eval(xLeft);
    if (f1Val * (f1->Eval(xRight) - f2->Eval(xRight)) >= 0) {
        cout << "Warning: No intersection found in range [" << xMin << ", " << xMax << "]. Using midpoint.\n";
        return -1.0;
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

// Function to compute chi2pid cut
double getChi2Cut(double p) {
    const double C = 0.9795;
    if (p <= 2.75) return 3.0 * C;
    else return C * (-0.20 + 3.00 * exp(-p / 4.51) + 37.26 * exp(-p / 0.87));
}

int main() {
    auto start = chrono::high_resolution_clock::now();
    ROOT::EnableImplicitMT(4);
    cout << "Multi-threading enabled for RDataFrame processing." << endl;

    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/Skim/pkptreeCxC_9_test_modified.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file pkptreeCxC_9_test_modified.root\n";
        ROOT::DisableImplicitMT();
        return 1;
    }
    ROOT::RDataFrame df_pions("EB_pid_pions", file);
    ROOT::RDataFrame df_kaons("EB_pid_kaons", file);

    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/beta_plots", kTRUE);

    vector<double> pBins;
    double pMin = 1.0, pMax = 7.0;
    int nBins = static_cast<int>((pMax - pMin) / 0.3); // 6 bins (1.0-1.3, 1.3-1.6, ..., 6.7-7.0)
    for (int i = 0; i <= nBins; i++) pBins.push_back(pMin + i * 0.3);

    TCanvas* canvas = new TCanvas("canvas", "Beta Distribution per Momentum Bin", 1200, 800);
    canvas->Print("output/beta_plots/beta_distributions.pdf["); // Updated output file name to avoid confusion
    const double chi2NegCut = -3.0 * 0.9795;

    ofstream txtFile("output/beta_plots/contamination_table.txt");
    txtFile << "\nTABLE X: Data-driven estimation for kaon contamination in pion sample\n";
    txtFile << "Bin\tp (GeV)\tN_pi/N_pi_tot (%)\tN_Ki/N_pi (%)\tN_Ki/N_pi_tot (%)\n";

    double total_pions = *df_pions.Count(); // Total uncut pions as N_pi_tot
    vector<double> n_pi_tot(nBins, 0.0), n_k_in_pi(nBins, 0.0), n_k_in_pi_tot(nBins, 0.0);

    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];

        auto filtered_df_pions_cut = df_pions.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"})
                                            .Filter([chi2NegCut](float p, float chi2) {
                                                double chi2Cut = getChi2Cut(p);
                                                return chi2 > chi2NegCut && chi2 < chi2Cut;
                                            }, {"p", "orig_chi2pid"});

        auto filtered_df_kaons = df_kaons.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});

        // Define histograms
        std::string pionHistName = "beta_pions_cut_" + to_string(i);
        std::string pionHistTitle = "p: [" + to_string(pLow) + "-" + to_string(pHigh) + ") GeV/c; #beta; Counts";
        ROOT::RDF::TH1DModel modelBetaPionsCut(pionHistName.c_str(), pionHistTitle.c_str(), 100, 0.96, 1.02);

        std::string kaonHistName = "beta_kaons_" + to_string(i);
        std::string kaonHistTitle = "p: [" + to_string(pLow) + "-" + to_string(pHigh) + ") GeV/c; #beta; Counts";
        ROOT::RDF::TH1DModel modelBetaKaons(kaonHistName.c_str(), kaonHistTitle.c_str(), 100, 0.96, 1.02);

        auto histo_pions_cut = filtered_df_pions_cut.Histo1D(modelBetaPionsCut, "beta");
        auto histo_kaons = filtered_df_kaons.Histo1D(modelBetaKaons, "beta");

        histo_pions_cut->SetStats(0);
        histo_kaons->SetStats(0);

        TF1* fit_pion = nullptr;
        TF1* fit_kaon = nullptr;
        double pMean = 0, pSigma = 0, pConst = 0, kMean = 0, kSigma = 0, kConst = 0;
        double intersectionPoint = -1.0, contamination = -1.0;

        if (histo_pions_cut->GetEntries() >= 10) {
            double peak = histo_pions_cut->GetMaximum();
            int binMax = histo_pions_cut->GetMaximumBin();
            double meanGuess = histo_pions_cut->GetBinCenter(binMax);
            double sigmaGuess = histo_pions_cut->GetRMS() * 0.5;

            fit_pion = new TF1("fit_pion", tripleGaussian, 0.986, 1.03, 9);
            fit_pion->SetParameters(peak * 0.8, meanGuess, sigmaGuess, peak * 0.05, 0.99, 0.002,
                                   peak * 0.1, meanGuess + 0.01, sigmaGuess * 1.5);
            fit_pion->SetParLimits(0, 0.0, peak * 1.2);
            fit_pion->SetParLimits(1, meanGuess - 0.01, meanGuess + 0.01);
            fit_pion->SetParLimits(2, sigmaGuess * 0.5, sigmaGuess * 2.0);
            fit_pion->SetParLimits(3, 0.0, peak * 0.1);
            fit_pion->SetParLimits(4, 0.986, 0.995);
            fit_pion->SetParLimits(5, 0.001, 0.003);
            fit_pion->SetParLimits(6, 0.0, peak * 0.5);
            fit_pion->SetParLimits(7, meanGuess + 0.005, 1.03);
            fit_pion->SetParLimits(8, sigmaGuess, sigmaGuess * 3.0);
            histo_pions_cut->Fit(fit_pion, "R", "", 0.986, 1.03);

            pMean = (fit_pion->GetParameter(1) * fit_pion->GetParameter(0) +
                     fit_pion->GetParameter(4) * fit_pion->GetParameter(3) +
                     fit_pion->GetParameter(7) * fit_pion->GetParameter(6)) /
                    (fit_pion->GetParameter(0) + fit_pion->GetParameter(3) + fit_pion->GetParameter(6));
            pSigma = TMath::Sqrt((fit_pion->GetParameter(2) * fit_pion->GetParameter(2) * fit_pion->GetParameter(0) +
                                  fit_pion->GetParameter(5) * fit_pion->GetParameter(5) * fit_pion->GetParameter(3) +
                                  fit_pion->GetParameter(8) * fit_pion->GetParameter(8) * fit_pion->GetParameter(6)) /
                                 (fit_pion->GetParameter(0) + fit_pion->GetParameter(3) + fit_pion->GetParameter(6)));
            pConst = fit_pion->GetParameter(0) + fit_pion->GetParameter(3) + fit_pion->GetParameter(6);
        }

        if (histo_kaons->GetEntries() >= 10) {
            double peak = histo_kaons->GetMaximum();
            int binMax = histo_kaons->GetMaximumBin();
            double meanGuess = histo_kaons->GetBinCenter(binMax);
            double sigmaGuess = histo_kaons->GetRMS() * 0.5;

            fit_kaon = new TF1("fit_kaon", gaussian, 0.96, 1.02, 3);
            fit_kaon->SetParameters(peak, meanGuess, sigmaGuess);
            fit_kaon->SetParLimits(0, 0.0, peak * 1.2);
            fit_kaon->SetParLimits(1, 0.96, 1.0);
            fit_kaon->SetParLimits(2, sigmaGuess * 0.5, sigmaGuess * 2.0);
            histo_kaons->Fit(fit_kaon, "R", "", 0.96, 1.02);

            kMean = fit_kaon->GetParameter(1);
            kSigma = fit_kaon->GetParameter(2);
            kConst = fit_kaon->GetParameter(0);
        }

        if (fit_pion && fit_kaon) {
            intersectionPoint = findIntersection(fit_kaon, fit_pion, kMean, pMean);
            double pionIntegral = fit_pion->Integral(0.96, 1.02);
            double kaonIntegral = fit_kaon->Integral(intersectionPoint, 1.02);
            if (pionIntegral > 0) contamination = (kaonIntegral / pionIntegral) * 100.0;
            else cout << "Error: Pion integral is zero or negative for bin [" << pLow << "-" << pHigh << "]\n";

            double n_pi = histo_pions_cut->GetEntries(); // N_pi for this bin
            double n_pi_fraction = n_pi / total_pions;
            double n_k_in_pi_fraction = contamination / 100.0; // N_Ki/N_pi as a fraction
            double n_k_in_pi_tot_fraction = n_k_in_pi_fraction * n_pi_fraction; // N_Ki/N_pi_tot

            n_pi_tot[i] = n_pi_fraction * 100.0; // Convert to percentage
            n_k_in_pi[i] = n_k_in_pi_fraction * 100.0; // Convert to percentage
            n_k_in_pi_tot[i] = n_k_in_pi_tot_fraction * 100.0; // Convert to percentage

            txtFile << (i + 1) << "\t[" << fixed << setprecision(2) << pLow << "," << pHigh << ")\t"
                    << fixed << setprecision(1) << n_pi_tot[i] << "\t"
                    << fixed << setprecision(2) << n_k_in_pi[i] << "\t"
                    << fixed << setprecision(2) << n_k_in_pi_tot[i] << endl;
        } else {
            txtFile << (i + 1) << "\t[" << fixed << setprecision(2) << pLow << "," << pHigh << ")\tN/A\tN/A\tN/A\n";
        }

        canvas->Clear();
        histo_pions_cut->SetMarkerStyle(20);
        histo_pions_cut->SetMarkerSize(1.2);
        histo_pions_cut->SetMarkerColor(kRed);
        histo_pions_cut->Draw("P");
        if (fit_pion) {
            fit_pion->SetLineColor(kRed);
            fit_pion->SetLineStyle(kDashed);
            fit_pion->Draw("SAME");
        }
        histo_kaons->SetMarkerStyle(20);
        histo_kaons->SetMarkerSize(1.2);
        histo_kaons->SetMarkerColor(kGreen);
        histo_kaons->Draw("P SAME");
        if (fit_kaon) {
            fit_kaon->SetLineColor(kGreen);
            fit_kaon->SetLineStyle(kDashed);
            fit_kaon->Draw("SAME");
        }

        // Removed canvas->SetLogy(); to use linear scale

        TLegend* leg = new TLegend(0.6, 0.6, 0.9, 0.9);
        if (fit_pion) {
            leg->AddEntry(&(*histo_pions_cut), "Pions (Cut)", "p");
            leg->AddEntry(fit_pion, "Triple Gaussian Fit", "l");
            leg->AddEntry((TObject*)0, TString::Format("A1: %.5f", fit_pion->GetParameter(0)), "");
            leg->AddEntry((TObject*)0, TString::Format("#mu1: %.5f", fit_pion->GetParameter(1)), "");
            leg->AddEntry((TObject*)0, TString::Format("#sigma1: %.5f", fit_pion->GetParameter(2)), "");
            leg->AddEntry((TObject*)0, TString::Format("A2: %.5f", fit_pion->GetParameter(3)), "");
            leg->AddEntry((TObject*)0, TString::Format("#mu2: %.5f", fit_pion->GetParameter(4)), "");
            leg->AddEntry((TObject*)0, TString::Format("#sigma2: %.5f", fit_pion->GetParameter(5)), "");
            leg->AddEntry((TObject*)0, TString::Format("A3: %.5f", fit_pion->GetParameter(6)), "");
            leg->AddEntry((TObject*)0, TString::Format("#mu3: %.5f", fit_pion->GetParameter(7)), "");
            leg->AddEntry((TObject*)0, TString::Format("#sigma3: %.5f", fit_pion->GetParameter(8)), "");
            leg->AddEntry((TObject*)0, TString::Format("#chi^{2}/NDF: %.2f", fit_pion->GetChisquare() / fit_pion->GetNDF()), "");
        }
        if (fit_kaon) {
            leg->AddEntry(&(*histo_kaons), "Kaons", "p");
            leg->AddEntry(fit_kaon, "Gaussian Fit", "l");
            leg->AddEntry((TObject*)0, TString::Format("A: %.5f", fit_kaon->GetParameter(0)), "");
            leg->AddEntry((TObject*)0, TString::Format("#mu: %.5f", fit_kaon->GetParameter(1)), "");
            leg->AddEntry((TObject*)0, TString::Format("#sigma: %.5f", fit_kaon->GetParameter(2)), "");
            leg->AddEntry((TObject*)0, TString::Format("#chi^{2}/NDF: %.2f", fit_kaon->GetChisquare() / fit_kaon->GetNDF()), "");
        }
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.02);
        leg->Draw();

        canvas->Update();
        canvas->Print("output/beta_plots/beta_distributions.pdf");
        delete leg;
        if (fit_pion) delete fit_pion;
        if (fit_kaon) delete fit_kaon;
    }

    canvas->Print("output/beta_plots/beta_distributions.pdf]");
    txtFile.close();
    delete canvas;
    file->Close();
    delete file;
    ROOT::DisableImplicitMT();
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Total program time: " << duration.count() / 1000.0 << " seconds" << endl;
    cout << "Program completed successfully." << endl;
    return 0;
}