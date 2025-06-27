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
    if (pMid >= 1.0 && pMid < 1.3) return {0.96, 0.94};    // 1.0-1.3 GeV/c
    if (pMid >= 1.3 && pMid < 1.6) return {0.96, 0.975};    // 1.3-1.6 GeV/c
    if (pMid >= 1.6 && pMid < 1.9) return {0.978, 0.982};    // 1.6-1.9 GeV/c
    if (pMid >= 1.9 && pMid < 2.2) return {0.98, 0.987};    // 1.9-2.2 GeV/c
    if (pMid >= 2.2 && pMid < 2.5) return {0.98, 0.99};     // 2.2-2.5 GeV/c
    if (pMid >= 2.5 && pMid < 2.8) return {0.98, 0.988};    // 2.5-2.8 GeV/c
    if (pMid >= 2.8 && pMid < 3.1) return {0.982, 0.989};   // 2.8-3.1 GeV/c
    if (pMid >= 3.1 && pMid < 3.4) return {0.984, 0.991};   // 3.1-3.4 GeV/c
    if (pMid >= 3.4 && pMid < 3.7) return {0.984, 0.993};   // 3.4-3.7 GeV/c
    if (pMid >= 3.7 && pMid < 4.0) return {0.984, 0.994};   // 3.7-4.0 GeV/c
    return {0.96, 1.03}; // Default range (not used here)
}

// Function to get manual pion fit range based on momentum bin
pair<double, double> getPionFitRange(double pMid) {
    if (pMid >= 1.0 && pMid < 1.3) return {0.988, 0.996};    // 1.0-1.3 GeV/c
    if (pMid >= 1.3 && pMid < 1.6) return {0.992, 1};       // 1.3-1.6 GeV/c
    if (pMid >= 1.6 && pMid < 1.9) return {0.993, 1.2};     // 1.6-1.9 GeV/c
    if (pMid >= 1.9 && pMid < 2.2) return {0.994, 1.2};     // 1.9-2.2 GeV/c
    if (pMid >= 2.2 && pMid < 2.5) return {0.995, 1.3};     // 2.2-2.5 GeV/c
    if (pMid >= 2.5 && pMid < 2.8) return {0.995, 1.3};     // 2.5-2.8 GeV/c
    if (pMid >= 2.8 && pMid < 3.1) return {0.996, 1.3};     // 2.8-3.1 GeV/c
    if (pMid >= 3.1 && pMid < 3.4) return {0.996, 1.4};     // 3.1-3.4 GeV/c
    if (pMid >= 3.4 && pMid < 3.7) return {0.996, 1.4};     // 3.4-3.7 GeV/c
    if (pMid >= 3.7 && pMid < 4.0) return {0.997, 1.4};     // 3.7-4.0 GeV/c
    return {0.96, 1.03}; // Default range (not used here)
}

double findIntersection(TF1* f1, TF1* f2, double xMin, double xMax, double betaPion, double betaKaon, double tolerance = 1e-6) {
    double midpoint = (betaPion + betaKaon) / 2.0;
    vector<double> intersections;
    double xLeft = xMin, xRight = xMax;
    double stepSize = (xMax - xMin) / 100.0; // Coarse scan to find potential intersections

    // Coarse scan to identify potential intersection regions
    for (double x = xMin; x <= xMax; x += stepSize) {
        double diff = f1->Eval(x) - f2->Eval(x);
        if (diff * (f1->Eval(x + stepSize) - f2->Eval(x + stepSize)) < 0) {
            // Refine using binary search in this region
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

    // Find the intersection closest to the midpoint
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

int main() {
    // Start total program timer
    auto start = chrono::high_resolution_clock::now();

    // Step 1: Enable multi-threading
    ROOT::EnableImplicitMT(16); // Adjust based on your CPU
    cout << "Multi-threading enabled for RDataFrame processing." << endl;

    // Step 2: Open the ROOT file and create RDataFrames
    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/Skim/pkptreeCxC_9_test_modified.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file pkptreeCxC_9_test_modified.root\n";
        ROOT::DisableImplicitMT();
        return 1;
    }
    ROOT::RDataFrame df_pions("EB_pid_pions", file); // Pion tree
    ROOT::RDataFrame df_kaons("EB_pid_kaons", file); // Kaon tree

    // Step 3: Create output directories
    gSystem->mkdir("output", kTRUE);
    
    // Step 4: Set up momentum bins (1 to 4 GeV, width 0.3)
    vector<double> pBins;
    double pMin = 1.0, pMax = 4.0;
    int nBins = static_cast<int>((pMax - pMin) / 0.3); // 10 bins (1.0-1.3, 1.3-1.6, ..., 3.7-4.0)
    for (int i = 0; i <= nBins; i++) {
        pBins.push_back(pMin + i * 0.3);
    }

    // Step 5: Open files to save fit parameters and contamination
    ofstream fitOutput("output/fit_parameters_upto_4.0.txt");
    fitOutput << "Bin\tpLow\tpHigh\tPion_Mean\tPion_Sigma\tPion_Amplitude\tPion_Chi2NDF\tKaon_Mean\tKaon_Sigma\tKaon_Amplitude\tKaon_Chi2NDF\n";
    ofstream contamOutput("output/contamination_upto_4.0.csv");
    contamOutput << "Momentum Bin (GeV/c),Intersection Point,Kaon-to-Pion Ratio (%)\n";

    // Step 6: Initialize canvases and print start
    TCanvas* canvas = new TCanvas("canvas", "Beta Distribution per Momentum Bin (up to 4.0 GeV/c)", 1200, 800);
    canvas->Print("output/beta_distributions_upto_4.0.pdf[");

    // Particle masses (GeV/c^2)
    const double pionMass = 0.1396; // Pion mass
    const double kaonMass = 0.4937; // Kaon mass

    // Step 7: Process each momentum bin
    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];
        double pMid = (pLow + pHigh) / 2.0; // Midpoint momentum for theoretical beta

        // Filter for momentum range and apply chi2pid cut
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

        // Create histogram models for beta
        ROOT::RDF::TH1DModel modelBetaPions(
            TString::Format("beta_pions_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #beta; Counts", pLow, pHigh),
            100, 0.96, 1.03 // Beta range set to 0.96-1.03
        );
        ROOT::RDF::TH1DModel modelBetaKaons(
            TString::Format("beta_kaons_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #beta; Counts", pLow, pHigh),
            100, 0.96, 1.03 // Beta range set to 0.96-1.03
        );

        // Fill histograms
        auto histo_pions = filtered_df_pions.Histo1D(modelBetaPions, "beta");
        auto histo_kaons = filtered_df_kaons.Histo1D(modelBetaKaons, "beta");

        // Calculate theoretical beta
        double betaPion = getTheoreticalBeta(pMid, pionMass);
        double betaKaon = getTheoreticalBeta(pMid, kaonMass);

        // Get manual fit ranges
        auto pionRange = getPionFitRange(pMid);
        auto kaonRange = getKaonFitRange(pMid);

        // Use exact manual fit ranges
        double kaonFitLow = kaonRange.first;
        double kaonFitHigh = kaonRange.second;
        double pionFitLow = max(0.96, pionRange.first);
        double pionFitHigh = min(1.03, pionRange.second);

        // Define Gaussian fit functions for pions
        TF1* fitPion = new TF1(TString::Format("fit_pion_%d", i), "gaus", pionFitLow, pionFitHigh);
        fitPion->SetParameter(0, histo_pions->GetMaximum()); // Amplitude
        fitPion->SetParameter(1, betaPion);                  // Initial mean
        fitPion->SetParameter(2, 0.005);                     // Initial sigma
        fitPion->SetLineColor(kBlue);
        fitPion->SetLineWidth(2);
        histo_pions->Fit(fitPion, "RQ");
        fitPion->SetRange(0.96, 1.03); // Display over full range

        // Define Gaussian fit functions for kaons
        TF1* fitKaon = new TF1(TString::Format("fit_kaon_%d", i), "gaus", kaonFitLow, kaonFitHigh);
        // Initialize amplitude and mean from maximum within the fit range
        int binMax = histo_kaons->GetMaximumBin();
        double xMax = histo_kaons->GetBinCenter(binMax);
        if (xMax >= kaonFitLow && xMax <= kaonFitHigh) {
            fitKaon->SetParameter(0, histo_kaons->GetBinContent(binMax)); // Amplitude from max within range
            fitKaon->SetParameter(1, xMax); // Initial mean from max within range
        } else {
            int bin_K = histo_kaons->FindBin(betaKaon);
            fitKaon->SetParameter(0, histo_kaons->GetBinContent(bin_K)); // Fall back to beta_K bin
            fitKaon->SetParameter(1, betaKaon); // Initial mean set to theoretical beta
        }
        fitKaon->SetParameter(2, 0.003);    // Initial sigma
        fitKaon->SetLineColor(kGreen);
        fitKaon->SetLineWidth(2);
        histo_kaons->Fit(fitKaon, "RME+", "", kaonFitLow, kaonFitHigh); // Use RME+ for better optimization
        fitKaon->SetRange(0.96, 1.03); // Display over full range

        // Save fit parameters
        fitOutput << i << "\t" << pLow << "\t" << pHigh << "\t"
                  << fitPion->GetParameter(1) << "\t" << fitPion->GetParameter(2) << "\t" << fitPion->GetParameter(0) << "\t" << fitPion->GetChisquare() / fitPion->GetNDF() << "\t"
                  << fitKaon->GetParameter(1) << "\t" << fitKaon->GetParameter(2) << "\t" << fitKaon->GetParameter(0) << "\t" << fitKaon->GetChisquare() / fitKaon->GetNDF() << "\n";

        // Calculate intersection and contamination
double intersectionPoint = -1.0;
double contamination = -1.0;
if (histo_pions->GetEntries() >= 10 && histo_kaons->GetEntries() >= 10) {
    // Find intersection between kaon and pion Gaussian fits
    intersectionPoint = findIntersection(fitKaon, fitPion, kaonFitLow, pionFitHigh, betaPion, betaKaon);

    // Integrate kaon Gaussian from intersection to upper limit (1.03)
    double kaonIntegral = fitKaon->Integral(intersectionPoint, 1.03);

    // Integrate pion Gaussian over its fit range
    double pionIntegral = fitPion->Integral(0.96, 1.03);

    // Compute contamination ratio
    if (pionIntegral > 0) {
        contamination = (kaonIntegral / pionIntegral) * 100.0;
    } else {
        cout << "Error: Pion integral is zero or negative for bin " << i << "\n";
    }

    cout << "Bin [" << pLow << "-" << pHigh << "): Intersection at beta = " << intersectionPoint
         << ", Kaon Integral = " << kaonIntegral << ", Pion Integral = " << pionIntegral
         << ", Contamination = " << contamination << "%\n";

    // Save to CSV
    contamOutput << pLow << "-" << pHigh << "," << intersectionPoint << "," << contamination << "\n";
} else {
    cout << "Bin [" << pLow << "-" << pHigh << "): Not enough entries to calculate contamination\n";
    contamOutput << pLow << "-" << pHigh << ",N/A,N/A\n";
}

        // Plot on the same canvas with point style
        canvas->Clear();
        histo_pions->SetMarkerColor(kBlue); // Pions in blue
        histo_pions->SetMarkerStyle(20);    // Circle markers
        histo_pions->SetMarkerSize(1.2);    // Slightly smaller marker size
        histo_pions->Draw("P");             // Point style
        histo_kaons->SetMarkerColor(kGreen); // Kaons in green
        histo_kaons->SetMarkerStyle(20);     // Circle markers (same as pions)
        histo_kaons->SetMarkerSize(1.2);     // Slightly smaller marker size
        histo_kaons->Draw("P SAME");         // Point style, same canvas

        // Draw fits
        fitPion->Draw("SAME");
        fitKaon->Draw("SAME");

        // Get maximum height for vertical lines
        double maxHeight = max(histo_pions->GetMaximum(), histo_kaons->GetMaximum()) * 1.1;

        // Add theoretical vertical lines
        TLine* linePion = new TLine(betaPion, 0, betaPion, maxHeight);
        linePion->SetLineColor(kBlue);
        linePion->SetLineStyle(2); // Dashed
        linePion->SetLineWidth(2);
        linePion->SetLineColorAlpha(kBlue, 0.3); // 30% opacity for fade
        linePion->Draw();

        TLine* lineKaon = new TLine(betaKaon, 0, betaKaon, maxHeight);
        lineKaon->SetLineColor(kGreen);
        lineKaon->SetLineStyle(2); // Dashed
        lineKaon->SetLineWidth(2);
        lineKaon->SetLineColorAlpha(kGreen, 0.3); // 30% opacity for fade
        lineKaon->Draw();

        // Add legend
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

        // Update and save
        canvas->Update();
        canvas->Print(TString::Format("output/beta_distributions_upto_4.0.pdf"));
        delete leg;
        delete linePion;
        delete lineKaon;
        delete fitPion;
        delete fitKaon;
    }

    // Step 8: Close PDF and output files
    canvas->Print("output/beta_distributions_upto_4.0.pdf]");
    fitOutput.close();
    contamOutput.close();

    // Step 9: Clean up
    delete canvas;
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