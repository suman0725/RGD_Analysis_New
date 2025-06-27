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

int main() {
    // Start total program timer
    auto start = chrono::high_resolution_clock::now();

    // Step 1: Enable multi-threading
    ROOT::EnableImplicitMT(15); // Adjust based on your CPU
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

    // Step 3: Create output directory
    gSystem->mkdir("output", kTRUE);
    
    // Step 4: Set up momentum bins (4 to 7 GeV, width 0.3)
    vector<double> pBins;
    double pMin = 4.0, pMax = 7.0;
    int nBins = static_cast<int>((pMax - pMin) / 0.3); // 10 bins (4.0-4.3, 4.3-4.6, ..., 6.7-7.0)
    for (int i = 0; i <= nBins; i++) {
        pBins.push_back(pMin + i * 0.3);
    }

    // Step 5: Open file to save fit parameters
    ofstream fitOutput("output/fit_parameters.txt");
    fitOutput << "Bin\tpLow\tpHigh\tPion_Mean\tPion_Sigma\tPion_Amplitude\tPion_Chi2NDF\tKaon_Mean\tKaon_Sigma\tKaon_Amplitude\tKaon_Chi2NDF\n";

    // Step 6: Initialize canvas and print start
    TCanvas* canvas = new TCanvas("canvas", "Beta Distribution per Momentum Bin", 1200, 800);
    canvas->Print("output/beta_distributions_chi2cut_fivesigma_linear.pdf[");

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

        // Define Gaussian fit for pions
        TF1* fitPion = new TF1(TString::Format("fit_pion_%d", i), "gaus", 0.96, 1.03);
        fitPion->SetParameter(0, histo_pions->GetMaximum()); // Amplitude
        fitPion->SetParameter(1, 0.999);                     // Initial mean
        fitPion->SetParameter(2, 0.005);                     // Initial sigma
        fitPion->SetLineColor(kBlue);
        fitPion->SetLineWidth(2);
        histo_pions->Fit(fitPion, "RQ");

        // Define Gaussian fit for kaons
        TF1* fitKaon = new TF1(TString::Format("fit_kaon_%d", i), "gaus", 0.96, 1.03);
        if (i == 0) { // 4.0-4.3 GeV/c
            fitKaon->SetParameters(histo_kaons->GetMaximum(), 0.986, 0.993);
        } else if (i == 1) { // 4.3-4.6 GeV/c
            fitKaon->SetParameters(histo_kaons->GetMaximum(), 0.994, 0.01);
        } else if (i == 2) { // 4.6-4.9 GeV/c
            fitKaon->SetParameters(histo_kaons->GetMaximum(), 0.994, 0.01);
        } else if (i == 3) { // 4.9-5.2 GeV/c
            fitKaon->SetParameters(histo_kaons->GetMaximum(), 0.996, 0.012);
        } else if (i == 4) { // 5.2-5.5 GeV/c
            fitKaon->SetParameters(histo_kaons->GetMaximum(), 0.994, 0.012);
        } else if (i == 5) { // 5.5-5.8 GeV/c
            fitKaon->SetParameters(histo_kaons->GetMaximum(), 0.994, 0.012);
        } else if (i == 6) { // 5.8-6.1 GeV/c
            fitKaon->SetParameters(histo_kaons->GetMaximum(), 0.994, 0.01);
        } else if (i == 7) { // 6.1-6.4 GeV/c
            fitKaon->SetParameters(histo_kaons->GetMaximum(), 0.994, 0.01);
        } else if (i == 8) { // 6.4-6.7 GeV/c
            fitKaon->SetParameters(histo_kaons->GetMaximum(), 0.994, 0.01);
        } else if (i == 9) { // 6.7-7.0 GeV/c
            fitKaon->SetParameters(histo_kaons->GetMaximum(), 0.994, 0.01);
        }
        fitKaon->SetLineColor(kGreen);
        fitKaon->SetLineWidth(2);
        histo_kaons->Fit(fitKaon, "RQ");

        // Save fit parameters
        fitOutput << i << "\t" << pLow << "\t" << pHigh << "\t"
                  << fitPion->GetParameter(1) << "\t" << fitPion->GetParameter(2) << "\t" << fitPion->GetParameter(0) << "\t" << fitPion->GetChisquare() / fitPion->GetNDF() << "\t"
                  << fitKaon->GetParameter(1) << "\t" << fitKaon->GetParameter(2) << "\t" << fitKaon->GetParameter(0) << "\t" << fitKaon->GetChisquare() / fitKaon->GetNDF() << "\n";

        // Plot on the same canvas with point style
        canvas->Clear();
        histo_pions->SetMarkerColor(kBlue); // Pions in blue
        histo_pions->SetMarkerStyle(20);    // Circle markers
        histo_pions->SetMarkerSize(1.2);    // Slightly smaller marker size
        histo_pions->Draw("P");             // Point style
        histo_kaons->SetMarkerColor(kGreen); // Kaons in green
        histo_kaons->SetMarkerStyle(20);     // Circle markers
        histo_kaons->SetMarkerSize(1.2);     // Slightly smaller marker size
        histo_kaons->Draw("P SAME");         // Point style

        // Draw fits
        fitPion->Draw("SAME");
        fitKaon->Draw("SAME");

        // Calculate theoretical beta
        double betaPion = getTheoreticalBeta(pMid, pionMass);
        double betaKaon = getTheoreticalBeta(pMid, kaonMass);

        // Get maximum height for vertical lines
        double maxHeight = max(histo_pions->GetMaximum(), histo_kaons->GetMaximum()) * 1.1;

        // Add theoretical vertical lines
        TLine* linePion = new TLine(betaPion, 0, betaPion, maxHeight);
        linePion->SetLineColor(kBlue);
        linePion->SetLineStyle(2); // Dashed
        linePion->SetLineWidth(2);
        linePion->SetLineColorAlpha(kBlue, 0.3); // 30% opacity
        linePion->Draw();

        TLine* lineKaon = new TLine(betaKaon, 0, betaKaon, maxHeight);
        lineKaon->SetLineColor(kGreen);
        lineKaon->SetLineStyle(2); // Dashed
        lineKaon->SetLineWidth(2);
        lineKaon->SetLineColorAlpha(kGreen, 0.3); // 30% opacity
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
        canvas->Print(TString::Format("output/beta_distributions_chi2cut_fivesigma_linear.pdf"));
        delete leg;
        delete linePion;
        delete lineKaon;
        delete fitPion;
        delete fitKaon;
    }

    // Step 8: Close PDF and output file
    canvas->Print("output/beta_distributions_chi2cut_fivesigma_linear.pdf]");
    fitOutput.close();

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