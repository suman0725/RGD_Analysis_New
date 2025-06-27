#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TF1.h>
#include <TLine.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <TStyle.h>

using namespace std;
using namespace ROOT;

int main() {
    auto start = chrono::high_resolution_clock::now();
    ROOT::EnableImplicitMT(16);
    cout << "Multi-threading enabled for RDataFrame processing." << endl;
    gStyle->SetOptStat(0);

    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/Skim/pkptreeCxC_9_test_modified.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file pkptreeCxC_9_test_modified.root\n";
        ROOT::DisableImplicitMT();
        return 1;
    }
    ROOT::RDataFrame df_all("EB_pos_pion_assumed", file);

    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/dt_plots", kTRUE);

    vector<double> pBins;
    double pMin = 0.1, pMax = 7.0;
    double binWidth = 0.3;
    int nBins = static_cast<int>((pMax - pMin) / binWidth);
    for (int i = 0; i <= nBins; i++) {
        pBins.push_back(pMin + i * binWidth);
    }

    TCanvas* canvas = new TCanvas("canvas", "Delta t Distribution per Momentum Bin", 1200, 800);
    canvas->Print("output/dt_plots/dt_distributions_linear.pdf[");

    ofstream outFile("output/dt_fit_parameters.txt");
    outFile << "p_min(GeV)\tmean(ns)\tsigma(ns)\tmean-3sigma(ns)\tmean+3sigma(ns)\n";

    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];
    
        auto filtered_df = df_all.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});
    
        ROOT::RDF::TH1DModel modelDt(
            TString::Format("dt_%d", i),
            TString::Format("p: [%.2f-%.2f) GeV/c; #Delta t (ns); Counts", pLow, pHigh),
            100, -2.0, 2.0
        );
    
        auto histo_dt = filtered_df.Histo1D(modelDt, "dt");
    
        canvas->Clear();
        histo_dt->SetMarkerStyle(20);
        histo_dt->SetMarkerSize(1.2);
        histo_dt->SetMarkerColor(kBlack);
        histo_dt->Draw("P");
    
        double mean = 0.0, sigma = 0.0;
        TLegend* leg = new TLegend(0.6, 0.6, 0.85, 0.85);
        leg->AddEntry(histo_dt.GetPtr(), Form("Entries: %.0f", histo_dt->GetEntries()), "p");
       
        TF1* gaussian = new TF1("gaussian", "gaus", -0.15, 0.15);
        gaussian->SetLineColor(kBlue);
        gaussian->SetParameter(1, 0.0);
        gaussian->SetParameter(2, 0.05);
        gaussian->SetNpx(1000);
        histo_dt->Fit(gaussian, "RQ");
        mean = gaussian->GetParameter(1);
        sigma = gaussian->GetParameter(2);
    
        gaussian->Draw("SAME");
        
        // Adding lines at mean - 3*sigma and mean + 3*sigma
        TLine* line_m3s = new TLine(mean - 3 * sigma, 0, mean - 3 * sigma, histo_dt->GetMaximum());
        TLine* line_p3s = new TLine(mean + 3 * sigma, 0, mean + 3 * sigma, histo_dt->GetMaximum());
        line_m3s->SetLineColor(kBlue);
        line_p3s->SetLineColor(kBlue);
        line_m3s->SetLineStyle(2);
        line_p3s->SetLineStyle(2);
        line_m3s->Draw("SAME");
        line_p3s->Draw("SAME");
    
        leg->AddEntry(gaussian, Form("Mean: %.3f ns", mean), "l");
        leg->AddEntry(gaussian, Form("Sigma: %.3f ns", sigma), "l");
        leg->AddEntry(line_m3s, "#mu #pm 3#sigma", "l");
        
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);
        leg->Draw();
    
        canvas->Update();
        canvas->Print(TString::Format("output/dt_plots/dt_distributions_linear.pdf"));
    
        // Delete objects after printing
        delete gaussian;
        delete line_m3s;
        delete line_p3s;
        delete leg;
    }

    canvas->Print("output/dt_plots/dt_distributions_linear.pdf]");
    outFile.close();
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