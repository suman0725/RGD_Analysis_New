#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TLine.h>
#include <TF1.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cmath>

using namespace std;
using namespace ROOT;

// chi2 PIC cut functions (momentum-dependent)
double getChi2CutNeg(double p) {
    return -4.606 + (-4.796) * exp(-2.766 * p); 
}

double getChi2CutPos(double p) {
    return 4.449 + 4.899 * exp(-2.366 * p); 
}

// Calculate theoretical beta
double getTheoreticalBeta(double p, double mass) {
    return p / sqrt(p * p + mass * mass);
}

// Function to get manual kaon fit range based on momentum bin
pair<double, double> getKaonFitRange(double pMid) {
    if (pMid >= 4.0 && pMid < 4.3) return {0.988, 0.994};   // 4.0-4.3 GeV/c
    if (pMid >= 4.3 && pMid < 4.6) return {0.989, 0.995};   // 4.3-4.6 GeV/c
    if (pMid >= 4.6 && pMid < 4.9) return {0.990, 0.996};   // 4.6-4.9 GeV/c
    if (pMid >= 4.9 && pMid < 5.2) return {0.991, 0.997};   // 4.9-5.2 GeV/c
    if (pMid >= 5.2 && pMid < 5.5) return {0.992, 0.997};   // 5.2-5.5 GeV/c
    if (pMid >= 5.5 && pMid < 5.8) return {0.993, 0.998};   // 5.5-5.8 GeV/c
    if (pMid >= 5.8 && pMid < 6.1) return {0.994, 0.998};   // 5.8-6.1 GeV/c
    if (pMid >= 6.1 && pMid < 6.4) return {0.994, 0.999};   // 6.1-6.4 GeV/c
    if (pMid >= 6.4 && pMid < 6.7) return {0.995, 0.999};   // 6.4-6.7 GeV/c
    if (pMid >= 6.7 && pMid < 7.0) return {0.996, 0.999};   // 6.7-7.0 GeV/c
    return {0.96, 1.03}; // Default range
}

// Function to get manual pion fit range based on momentum bin
pair<double, double> getPionFitRange(double pMid) {
    if (pMid >= 4.0 && pMid < 4.3) return {0.997, 1.015};    // 4.0-4.3 GeV/c
    if (pMid >= 4.3 && pMid < 4.6) return {0.997, 1.02};    // 4.3-4.6 GeV/c
    if (pMid >= 4.6 && pMid < 4.9) return {0.997, 1.02};    // 4.6-4.9 GeV/c
    if (pMid >= 4.9 && pMid < 5.2) return {0.998, 1.02};    // 4.9-5.2 GeV/c
    if (pMid >= 5.2 && pMid < 5.5) return {0.998, 1.02};    // 5.2-5.5 GeV/c
    if (pMid >= 5.5 && pMid < 5.8) return {0.998, 1.02};    // 5.5-5.8 GeV/c
    if (pMid >= 5.8 && pMid < 6.1) return {0.999, 1.02};    // 5.8-6.1 GeV/c
    if (pMid >= 6.1 && pMid < 6.4) return {0.999, 1.02};    // 6.1-6.4 GeV/c
    if (pMid >= 6.4 && pMid < 6.7) return {0.999, 1.02};    // 6.4-6.7 GeV/c
    if (pMid >= 6.7 && pMid < 7.0) return {0.999, 1.02};    // 6.7-7.0 GeV/c
    return {0.96, 1.03}; // Default range
}

double findIntersection(TF1* f1, TF1* f2, double xMin, double xMax, double betaPion, double betaKaon, double tolerance = 1e-6) {
    double midpoint = (betaPion + betaKaon) / 2.0;
    vector<double> intersections;
    double xLeft = xMin, xRight = xMax;
    double stepSize = (xMax - xMin) / 1000.0; // Increased precision

    for (double x = xMin; x <= xMax; x += stepSize) {
        double diff = f1->Eval(x) - f2->Eval(x);
        if (diff * (f1->Eval(x + stepSize) - f2->Eval(x + stepSize)) < 0) {
            double localXLeft = max(xMin, x - stepSize);
            double localXRight = min(xMax, x + stepSize);
            while (localXRight - localXLeft > tolerance) {
                double xMid = (localXLeft + localXRight) / 2.0;
                double f1Val = f1->Eval(xMid) - f2->Eval(xMid);
                if (fabs(f1Val) < tolerance) {
                    intersections.push_back(xMid);
                    break;
                }
                if (f1Val * (f1->Eval(localXLeft) - f2->Eval(localXLeft)) < 0) {
                    localXRight = xMid;
                } else {
                    localXLeft = xMid;
                }
            }
        }
    }

    if (intersections.empty()) {
        cout << "Warning: No intersection found in range [" << xMin << ", " << xMax << "]. Using midpoint.\n";
        return midpoint;
    }

    double closestIntersection = intersections[0];
    double minDistance = fabs(intersections[0] - midpoint);
    for (size_t j = 1; j < intersections.size(); j++) {
        double distance = fabs(intersections[j] - midpoint);
        if (distance < minDistance) {
            minDistance = distance;
            closestIntersection = intersections[j];
        }
    }
    return closestIntersection;
}

// Double Crystal Ball function
Double_t doubleCrystalBall(Double_t* x, Double_t* par) {
    Double_t t = (x[0] - par[1]) / par[2];
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

int main() {
    auto start = chrono::high_resolution_clock::now();

    ROOT::EnableImplicitMT(16);
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

    vector<double> pBins;
    double pMin = 4.0, pMax = 7.0;
    int nBins = static_cast<int>((pMax - pMin) / 0.3);
    for (int i = 0; i <= nBins; i++) {
        pBins.push_back(pMin + i * 0.3);
    }

    ofstream fitOutput("output/fit_parameters_4.0_to_7.0.txt");
    fitOutput << "Bin\tpLow\tpHigh\tPion_Mean\tPion_Sigma\tPion_Amplitude\tPion_Chi2NDF\tKaon_Mean\tKaon_Sigma\tKaon_Amplitude\tKaon_Chi2NDF\n";
    ofstream contamOutput("output/contamination_4.0_to_7.0.csv");
    contamOutput << "Momentum Bin (GeV/c),Intersection Point,Kaon Integral,Pion Integral,Kaon-to-Pion Ratio (%)\n";

    TCanvas* canvas = new TCanvas("canvas", "Beta Distribution per Momentum Bin (4.0 to 7.0 GeV/c)", 1200, 800);
    canvas->Print("output/beta_distributions_4.0_to_7.0.pdf[");

    const double pionMass = 0.1396;
    const double kaonMass = 0.4937;

    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];
        double pMid = (pLow + pHigh) / 2.0;

        auto filtered_df_pions = df_pions.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"})
                                        .Filter([](float p, float chi2) {
                                            double chi2Min = getChi2CutNeg(p);
                                            double chi2Max = getChi2CutPos(p);
                                            return chi2 > chi2Min && chi2 < chi2Max;
                                        }, {"p", "recomputed_chi2pid"});
        auto filtered_df_kaons = df_kaons.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"})
                                        .Filter([](float p, float chi2) {
                                            double chi2Min = getChi2CutNeg(p);
                                            double chi2Max = getChi2CutPos(p);
                                            return chi2 > chi2Min && chi2 < chi2Max;
                                        }, {"p", "recomputed_chi2pid"});

        ROOT::RDF::TH1DModel modelBetaPions(
            TString::Format("beta_pions_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #beta; Counts", pLow, pHigh),
            100, 0.96, 1.03
        );
        ROOT::RDF::TH1DModel modelBetaKaons(
            TString::Format("beta_kaons_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #beta; Counts", pLow, pHigh),
            100, 0.96, 1.03
        );

        auto histo_pions = filtered_df_pions.Histo1D(modelBetaPions, "beta");
        histo_pions->SetStats(0);
        auto histo_kaons = filtered_df_kaons.Histo1D(modelBetaKaons, "beta");
        histo_kaons->SetStats(0);

        double betaPion = getTheoreticalBeta(pMid, pionMass);
        double betaKaon = getTheoreticalBeta(pMid, kaonMass);

        auto pionRange = getPionFitRange(pMid);
        auto kaonRange = getKaonFitRange(pMid);

        double kaonFitLow = kaonRange.first;
        double kaonFitHigh = kaonRange.second;
        double pionFitLow = pionRange.first;
        double pionFitHigh = pionRange.second;

        TF1* fitPion = new TF1(TString::Format("fit_pion_%d", i), doubleCrystalBall, pionFitLow, pionFitHigh, 8);
        fitPion->SetParameter(0, histo_pions->GetMaximum() * 0.7);
        fitPion->FixParameter(1, betaPion);
        fitPion->SetParameter(2, 0.005);
        fitPion->SetParameter(3, 1.5);
        fitPion->SetParameter(4, 2.0);
        fitPion->SetParameter(5, histo_pions->GetMaximum() * 0.3);
        fitPion->SetParameter(6, 1.5);
        fitPion->SetParameter(7, 2.0);
        fitPion->SetParLimits(0, 0, histo_pions->GetMaximum() * 1.2);
        fitPion->SetParLimits(2, 0.001, 0.01);
        fitPion->SetParLimits(3, 0.5, 5.0);
        fitPion->SetParLimits(4, 1.0, 30.0);
        fitPion->SetParLimits(5, 0, histo_pions->GetMaximum() * 0.6);
        fitPion->SetParLimits(6, 0.5, 5.0);
        fitPion->SetParLimits(7, 1.0, 30.0);
        fitPion->SetLineColor(kBlue);
        fitPion->SetLineWidth(2);
        histo_pions->Fit(fitPion, "RQ");
        fitPion->SetRange(0.96, 1.03);

        TF1* fitKaon = new TF1(TString::Format("fit_kaon_%d", i), "gaus", kaonFitLow, kaonFitHigh);
        int binMax = histo_kaons->GetMaximumBin();
        fitKaon->SetParameter(0, histo_kaons->GetBinContent(binMax));
        fitKaon->FixParameter(1, betaKaon);
        fitKaon->SetParameter(2, histo_kaons->GetRMS()); // Use RMS for initial sigma
        fitKaon->SetParLimits(2, 0.001, 0.01); // Reasonable sigma limits
        fitKaon->SetLineColor(kGreen);
        fitKaon->SetLineWidth(2);
        histo_kaons->Fit(fitKaon, "RME+", "", kaonFitLow, kaonFitHigh);
        fitKaon->SetRange(0.96, 1.03);

        fitOutput << i << "\t" << pLow << "\t" << pHigh << "\t"
                  << fitPion->GetParameter(1) << "\t" << fitPion->GetParameter(2) << "\t" << (fitPion->GetParameter(0) + fitPion->GetParameter(5)) << "\t" << fitPion->GetChisquare() / fitPion->GetNDF() << "\t"
                  << fitKaon->GetParameter(1) << "\t" << fitKaon->GetParameter(2) << "\t" << fitKaon->GetParameter(0) << "\t" << fitKaon->GetChisquare() / fitKaon->GetNDF() << "\n";

        double intersectionPoint = -1.0;
        double contamination = -1.0;
        if (histo_pions->GetEntries() >= 10 && histo_kaons->GetEntries() >= 10) {
            intersectionPoint = findIntersection(fitKaon, fitPion, kaonFitLow, pionFitHigh, betaPion, betaKaon);
            double kaonIntegral = fitKaon->Integral(intersectionPoint, 1.03);
            double pionIntegral = fitPion->Integral(0.96, 1.03);
            if (pionIntegral > 0) {
                contamination = (kaonIntegral / pionIntegral) * 100.0;
            } else {
                cout << "Error: Pion integral is zero or negative for bin " << i << "\n";
            }
            cout << "Bin [" << pLow << "-" << pHigh << "): Intersection at beta = " << intersectionPoint
                 << ", Kaon Integral = " << kaonIntegral << ", Pion Integral = " << pionIntegral
                 << ", Contamination = " << contamination << "%\n";
            contamOutput << pLow << "-" << pHigh << "," << intersectionPoint << "," << kaonIntegral << "," << pionIntegral << "," << contamination << "\n";
        } else {
            cout << "Bin [" << pLow << "-" << pHigh << "): Not enough entries to calculate contamination\n";
            contamOutput << pLow << "-" << pHigh << ",N/A,N/A,N/A,N/A\n";
        }

        canvas->Clear();
        histo_pions->SetMarkerColor(kBlue);
        histo_pions->SetMarkerStyle(20);
        histo_pions->SetMarkerSize(1.2);
        histo_pions->Draw("P");
        histo_kaons->SetMarkerColor(kGreen);
        histo_kaons->SetMarkerStyle(20);
        histo_kaons->SetMarkerSize(1.2);
        histo_kaons->Draw("P SAME");
        fitPion->Draw("SAME");
        fitKaon->Draw("SAME");

        double maxHeight = max(histo_pions->GetMaximum(), histo_kaons->GetMaximum()) * 1.1;
        TLine* linePion = new TLine(betaPion, 0, betaPion, maxHeight);
        linePion->SetLineColor(kBlue);
        linePion->SetLineStyle(2);
        linePion->SetLineWidth(2);
        linePion->SetLineColorAlpha(kBlue, 0.3);
        linePion->Draw();
        TLine* lineKaon = new TLine(betaKaon, 0, betaKaon, maxHeight);
        lineKaon->SetLineColor(kGreen);
        lineKaon->SetLineStyle(2);
        lineKaon->SetLineWidth(2);
        lineKaon->SetLineColorAlpha(kGreen, 0.3);
        lineKaon->Draw();

        TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg->AddEntry(histo_pions.GetPtr(), "Pions", "p");
        leg->AddEntry(fitPion, "Pion Fit", "l");
        leg->AddEntry(histo_kaons.GetPtr(), "Kaons", "p");
        leg->AddEntry(fitKaon, "Kaon Fit", "l");
        leg->AddEntry(linePion, "Pion #beta_{theory}", "l");
        leg->AddEntry(lineKaon, "Kaon #beta_{theory}", "l");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);
        leg->Draw();

        canvas->Update();
        canvas->Print(TString::Format("output/beta_distributions_4.0_to_7.0.pdf"));
        delete leg;
        delete linePion;
        delete lineKaon;
        delete fitPion;
        delete fitKaon;
    }

    canvas->Print("output/beta_distributions_4.0_to_7.0.pdf]");
    fitOutput.close();
    contamOutput.close();

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