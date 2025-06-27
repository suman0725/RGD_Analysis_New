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

// Function to get manual Proton fit range based on momentum bin
pair<double, double> getProtonFitRange(double pMid) {
    if (pMid >= 4.0 && pMid < 4.3) return {0.98, 0.986};   // 4.0-4.3 GeV/c
    if (pMid >= 4.3 && pMid < 4.6) return {0.982, 0.987};   // 4.3-4.6 GeV/c
    if (pMid >= 4.6 && pMid < 4.9) return {0.983, 0.989};   // 4.6-4.9 GeV/c
    if (pMid >= 4.9 && pMid < 5.2) return {0.984, 0.9883};   // 4.9-5.2 GeV/c
    if (pMid >= 5.2 && pMid < 5.5) return {0.984, 0.99};   // 5.2-5.5 GeV/c
    if (pMid >= 5.5 && pMid < 5.8) return {0.985, 0.99};   // 5.5-5.8 GeV/c
    if (pMid >= 5.8 && pMid < 6.1) return {0.986, 0.992};   // 5.8-6.1 GeV/c
    if (pMid >= 6.1 && pMid < 6.4) return {0.986, 0.992};   // 6.1-6.4 GeV/c
    if (pMid >= 6.4 && pMid < 6.7) return {0.9865, 0.993};  // 6.4-6.7 GeV/c
    if (pMid >= 6.7 && pMid < 7.0) return {0.988, 0.994};   // 6.7-7.0 GeV/c
    // Added ranges for 1.0–4.0 GeV/c (placeholders, update later)
    if (pMid >= 1.0 && pMid < 1.3) return {0.96, 0.962};
    if (pMid >= 1.3 && pMid < 1.6) return {0.966, 0.974};
    if (pMid >= 1.6 && pMid < 1.9) return {0.976, 0.982};
    if (pMid >= 1.9 && pMid < 2.2) return {0.98, 0.987};
    if (pMid >= 2.2 && pMid < 2.5) return {0.98, 0.99};
    if (pMid >= 2.5 && pMid < 2.8) return {0.96, 0.965};
    if (pMid >= 2.8 && pMid < 3.1) return {0.983, 0.9915};
    if (pMid >= 3.1 && pMid < 3.4) return {0.974, 0.976};
    if (pMid >= 3.4 && pMid < 3.7) return {0.976, 0.978};
    if (pMid >= 3.7 && pMid < 4.0) return {0.978, 0.982};
    return {0.96, 1.03}; // Default range
}

// Function to get manual pion fit range based on momentum bin
pair<double, double> getPionFitRange(double pMid) {
    if (pMid >= 4.0 && pMid < 4.3) return {0.988, 1.0103};    // 4.0-4.3 GeV/c
    if (pMid >= 4.3 && pMid < 4.6) return {0.997, 1.015};    // 4.3-4.6 GeV/c
    if (pMid >= 4.6 && pMid < 4.9) return {0.997, 1.015};    // 4.6-4.9 GeV/c
    if (pMid >= 4.9 && pMid < 5.2) return {0.999, 1.015};    // 4.9-5.2 GeV/c
    if (pMid >= 5.2 && pMid < 5.5) return {0.998, 1.015};    // 5.2-5.5 GeV/c
    if (pMid >= 5.5 && pMid < 5.8) return {0.998, 1.02};    // 5.5-5.8 GeV/c
    if (pMid >= 5.8 && pMid < 6.1) return {0.999, 1.02};    // 5.8-6.1 GeV/c
    if (pMid >= 6.1 && pMid < 6.4) return {0.999, 1.02};    // 6.1-6.4 GeV/c
    if (pMid >= 6.4 && pMid < 6.7) return {0.999, 1.02};    // 6.4-6.7 GeV/c
    if (pMid >= 6.7 && pMid < 7.0) return {0.999, 1.02};    // 6.7-7.0 GeV/c
    // Placeholder ranges for 1.0–4.0 GeV/c (update later)
    if (pMid >= 1.0 && pMid < 1.3) return {0.98, 1.015};
    if (pMid >= 1.3 && pMid < 1.6) return {0.982, 1.015};
    if (pMid >= 1.6 && pMid < 1.9) return {0.98, 1.015};
    if (pMid >= 1.9 && pMid < 2.2) return {0.984, 1.01};
    if (pMid >= 2.2 && pMid < 2.5) return {0.987, 1.01};
    if (pMid >= 2.5 && pMid < 2.8) return {0.99, 1.014};
    if (pMid >= 2.8 && pMid < 3.1) return {0.995, 1.015};
    if (pMid >= 3.1 && pMid < 3.4) return {0.99, 1.015};
    if (pMid >= 3.4 && pMid < 3.7) return {0.99, 1.015};
    if (pMid >= 3.7 && pMid < 4.0) return {0.988, 1.0103};
    return {0.96, 1.03}; // Default range
}

// Function to get intersection search range
pair<double, double> getIntersectionRange(double pMid) {
    if (pMid >= 4.0 && pMid < 4.3) return {0.983, 0.985};   // Adjusted to match visual crossover
    if (pMid >= 4.3 && pMid < 4.6) return {0.9855, 0.986};   // Adjusted to match visual crossover
    if (pMid >= 4.6 && pMid < 4.9) return {0.987, 0.989};   // Adjust based on visual inspection
    if (pMid >= 4.9 && pMid < 5.2) return {0.9882, 0.99};   // Adjust based on visual inspection
    if (pMid >= 5.2 && pMid < 5.5) return {0.99, 0.993}; 
      // Adjust based on visual inspection
    if (pMid >= 5.5 && pMid < 5.8) return {0.9918, 0.992};   // Adjust based on visual inspection
    if (pMid >= 5.8 && pMid < 6.1) return {0.9927, 0.9933};   // Adjust based on visual inspection
    if (pMid >= 6.1 && pMid < 6.4) return {0.9927, 0.9932};
       // Adjust based on visual inspection
    if (pMid >= 6.4 && pMid < 6.7) return {0.993, 0.994};   // Adjust based on visual inspection
    if (pMid >= 6.7 && pMid < 7.0) return {0.994, 0.995};   // Adjust based on visual inspection // Adjust based on visual inspection
    // Placeholder ranges for 1.0–4.0 GeV/c (update later)
    if (pMid >= 1.0 && pMid < 1.3) return {0.968, 0.972};
    if (pMid >= 1.3 && pMid < 1.6) return {0.973, 0.975};
    if (pMid >= 1.6 && pMid < 1.9) return {0.976, 0.98};
    if (pMid >= 1.9 && pMid < 2.2) return {0.983, 0.985};
    if (pMid >= 2.2 && pMid < 2.5) return {0.987, 0.988};
    if (pMid >= 2.5 && pMid < 2.8) return {0.989, 0.991};
    if (pMid >= 2.8 && pMid < 3.1) return {0.9918, 0.9919};
    if (pMid >= 3.1 && pMid < 3.4) return {0.9927, 0.993};
    if (pMid >= 3.4 && pMid < 3.7) return {0.9939, 0.9941};
    if (pMid >= 3.7 && pMid < 4.0) return {0.994, 0.996};
    return {0.995, 0.999}; // Default range
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

// Find intersection point
double findIntersection(TF1* f1, TF1* f2, double pMid, double tolerance = 1e-6) {
    auto intersectRange = getIntersectionRange(pMid);
    double xMinIntersect = intersectRange.first;
    double xMaxIntersect = intersectRange.second;

    vector<double> intersections;
    double stepSize = (xMaxIntersect - xMinIntersect) / 100.0;

    for (double x = xMinIntersect; x <= xMaxIntersect; x += stepSize) {
        double diff = f1->Eval(x) - f2->Eval(x);
        if (diff * (f1->Eval(x + stepSize) - f2->Eval(x + stepSize)) < 0) {
            double localXLeft = max(xMinIntersect, x - stepSize);
            double localXRight = min(xMaxIntersect, x + stepSize);
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
        cout << "Warning: No intersection found in range [" << xMinIntersect << ", " << xMaxIntersect << "]. Using midpoint.\n";
        return (xMinIntersect + xMaxIntersect) / 2.0;
    }

    return intersections[0]; // Return the first intersection found
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
    ROOT::RDataFrame df_protons("EB_pid_protons", file);

    gSystem->mkdir("output", kTRUE);

    vector<double> pBins;
    double pMin = 1.0, pMax = 7.0;
    int nBins = static_cast<int>((pMax - pMin) / 0.3);
    for (int i = 0; i <= nBins; i++) {
        pBins.push_back(pMin + i * 0.3);
    }

    ofstream fitOutput("output/fit_parameters_1.0_to_7.0.txt");
    fitOutput << "Bin\tpLow\tpHigh\tPion_Mean\tPion_Sigma\tPion_Amplitude\tPion_Chi2NDF\tProton_Mean\tProton_Sigma\tProton_Amplitude\tProton_Chi2NDF\n";
    ofstream contamOutput("output/contamination_1.0_to_7.0.csv");
    contamOutput << "Momentum Bin (GeV/c),Intersection Point,Proton-to-Pion Ratio (%)\n";

    TCanvas* canvas = new TCanvas("canvas", "Beta Distribution per Momentum Bin (1.0 to 7.0 GeV/c)", 1200, 800);
    canvas->Print("output/beta_distributions_1.0_to_7.0.pdf[");

    const double pionMass = 0.1396;
    const double ProtonMass = 0.93827;

    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];
        double pMid = (pLow + pHigh) / 2.0;

        auto filtered_df_pions = df_pions.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"})
                                        .Filter([](float p, float chi2) {
                                            double chi2Min = getChi2CutNeg(p);
                                            double chi2Max = getChi2CutPos(p);
                                            return chi2 > chi2Min && chi2 < chi2Max;
                                        }, {"p", "recomputed_chi2pid"});
        auto filtered_df_protons = df_protons.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"})
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
        ROOT::RDF::TH1DModel modelBetaprotons(
            TString::Format("beta_protons_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #beta; Counts", pLow, pHigh),
            100, 0.96, 1.03
        );

        auto histo_pions = filtered_df_pions.Histo1D(modelBetaPions, "beta");
        histo_pions->SetStats(0);
        histo_pions->Sumw2();
        auto histo_protons = filtered_df_protons.Histo1D(modelBetaprotons, "beta");
        histo_protons->SetStats(0);
        histo_protons->Sumw2();

        double betaPion = getTheoreticalBeta(pMid, pionMass);
        double betaProton = getTheoreticalBeta(pMid, ProtonMass);

        auto pionRange = getPionFitRange(pMid);
        auto ProtonRange = getProtonFitRange(pMid);

        double ProtonFitLow = ProtonRange.first;
        double ProtonFitHigh = ProtonRange.second;
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
        fitPion->SetNpx(1000);
        histo_pions->Fit(fitPion, "RQ");
        fitPion->SetRange(0.96, 1.03);

        TF1* fitProton = new TF1(TString::Format("fit_Proton_%d", i), "gaus", ProtonFitLow, ProtonFitHigh);
        fitProton->SetParameter(0, histo_protons->GetMaximum());
        if (i == 3) {
            fitProton->FixParameter(1, 0.9961);
            fitProton->SetParameter(2, 0.003);
        } else if (i == 0) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.003);
        } else if (i == 1) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.003);
        } else if (i == 2) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.003);
        } else if (i == 3) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.003);
        } else if (i == 4) {
            fitProton->FixParameter(1, 0.9966);
            fitProton->SetParameter(2, 0.003);
        } else if (i == 5) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.005);
        } else if (i == 6) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.005);
        } else if (i == 7) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.005);
        } else if (i == 8) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.005);
        } else if (i == 9) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.005);
        } 
        else if (i == 10) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.005);
        } else if (i == 11) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.005);
        } else if (i == 12) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.005);
        } else if (i == 13) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.005);
        } else if (i == 14) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.005);
        } else if (i == 15) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.005);
        } else if (i == 16) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.005);
        } else if (i == 17) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.005);
        } else if (i == 18) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.005);
        } else if (i == 19) {
            fitProton->FixParameter(1, betaProton);
            fitProton->SetParameter(2, 0.005);
        }
        
        else {
            fitProton->SetParameter(1, betaProton);
            fitProton->SetParameter(2, 0.003);
        }
        fitProton->SetLineColor(kRed);
        fitProton->SetLineWidth(2);
        histo_protons->Fit(fitProton, "RME+", "", ProtonFitLow, ProtonFitHigh);
        fitProton->SetRange(0.96, 1.03);

        fitOutput << i << "\t" << pLow << "\t" << pHigh << "\t"
                  << fitPion->GetParameter(1) << "\t" << fitPion->GetParameter(2) << "\t" << (fitPion->GetParameter(0) + fitPion->GetParameter(5)) << "\t" << fitPion->GetChisquare() / fitPion->GetNDF() << "\t"
                  << fitProton->GetParameter(1) << "\t" << fitProton->GetParameter(2) << "\t" << fitProton->GetParameter(0) << "\t" << fitProton->GetChisquare() / fitProton->GetNDF() << "\n";

        double intersectionPoint = -1.0;
        double contamination = -1.0;
        if (histo_pions->GetEntries() >= 10 && histo_protons->GetEntries() >= 10) {
            intersectionPoint = findIntersection(fitProton, fitPion, pMid);
            double ProtonIntegral = fitProton->Integral(intersectionPoint, 1.03);
            double pionIntegral = fitPion->Integral(0.96, 1.03);
            if (pionIntegral > 0) {
                contamination = (ProtonIntegral / pionIntegral) * 100.0;
            } else {
                cout << "Error: Pion integral is zero or negative for bin " << i << "\n";
            }
            cout << "Bin [" << pLow << "-" << pHigh << "): Intersection at beta = " << intersectionPoint
                 << ", Proton Integral = " << ProtonIntegral << ", Pion Integral = " << pionIntegral
                 << ", Contamination = " << contamination << "%\n";
            contamOutput << pLow << "-" << pHigh << "," << intersectionPoint << "," << contamination << "\n";
        } else {
            cout << "Bin [" << pLow << "-" << pHigh << "): Not enough entries to calculate contamination\n";
            contamOutput << pLow << "-" << pHigh << ",N/A,N/A\n";
        }

        canvas->Clear();
        histo_pions->SetMarkerColor(kBlue);
        histo_pions->SetMarkerStyle(20);
        histo_pions->SetMarkerSize(1.2);
        histo_pions->Draw("P");
        histo_protons->SetMarkerColor(kRed);
        histo_protons->SetMarkerStyle(20);
        histo_protons->SetMarkerSize(1.2);
        histo_protons->Draw("P SAME");
        fitPion->Draw("SAME");
        fitProton->Draw("SAME");

        double maxHeight = max(histo_pions->GetMaximum(), histo_protons->GetMaximum()) * 1.1;
        TLine* linePion = new TLine(betaPion, 0, betaPion, maxHeight);
        linePion->SetLineColor(kBlue);
        linePion->SetLineStyle(2);
        linePion->SetLineWidth(2);
        linePion->SetLineColorAlpha(kBlue, 0.3);
        linePion->Draw();
        TLine* lineProton = new TLine(betaProton, 0, betaProton, maxHeight);
        lineProton->SetLineColor(kRed);
        lineProton->SetLineStyle(2);
        lineProton->SetLineWidth(2);
        lineProton->SetLineColorAlpha(kRed, 0.3);
        lineProton->Draw();

        TLine* lineIntersection = new TLine(intersectionPoint, 0, intersectionPoint, maxHeight);
        lineIntersection->SetLineColor(kRed);
        lineIntersection->SetLineStyle(2);
        lineIntersection->SetLineWidth(2);
        //lineIntersection->Draw();

        TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg->AddEntry(histo_pions.GetPtr(), "Pions", "p");
        leg->AddEntry(fitPion, "Pion Fit", "l");
        leg->AddEntry(histo_protons.GetPtr(), "protons", "p");
        leg->AddEntry(fitProton, "Proton Fit", "l");
        leg->AddEntry(linePion, "Pion #beta_{theory}", "l");
        leg->AddEntry(lineProton, "Proton #beta_{theory}", "l");
        //leg->AddEntry(lineIntersection, "Intersection", "l");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);
        leg->Draw();

        canvas->Update();
        canvas->Print(TString::Format("output/beta_distributions_1.0_to_7.0.pdf"));
        delete leg;
        delete linePion;
        delete lineProton;
        delete lineIntersection;
        delete fitPion;
        delete fitProton;
    }

    canvas->Print("output/beta_distributions_1.0_to_7.0.pdf]");
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