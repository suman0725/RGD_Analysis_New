#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>
#include <TF1.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <iostream>
#include <chrono>
#include <vector>
#include <sstream>
#include <iomanip>

using namespace std;
using namespace ROOT;

int main() {
    // Start total program timer
    auto start = chrono::high_resolution_clock::now();

    // Enable multi-threading
    ROOT::EnableImplicitMT(4);
    cout << "Multi-threading enabled for RDataFrame processing." << endl;

    // Set style to remove stat box and enable log scale on z-axis
    gStyle->SetOptStat(0);
    gStyle->SetOptLogz(1);

    // Create output directory
    gSystem->mkdir("output", kTRUE);

    // Open the ROOT file and create RDataFrame
    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/Skim/pkptreeCxC_9_test_modified.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file pkptreeCxC_9_test_modified.root\n";
        ROOT::DisableImplicitMT();
        return 1;
    }
    ROOT::RDataFrame df_pion("EB_all_pion_assumed", file);

    // Create 2D histogram model for dt vs p
    auto h2_dt_p = df_pion.Histo2D(
        {"h2_dt_p", "chi2pid vs Momentum; Momentum (GeV/c); chi2pid; Counts", 
         200, 0.0, 7.0,  
         200, -15.0, 15.0}, 
        "p", "recomputed_chi2pid");

    // Initialize canvas
    TCanvas* canvas = new TCanvas("canvas", "chi2pid vs Momentum 2D Plot", 800, 600);
    
    // Draw 2D histogram
    h2_dt_p->Draw("COLZ");

    // Updated vectors with data up to 3.7 GeV/c
   #include <vector>
using namespace std;

vector<double> p = {
    0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55,
    1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75,
    2.85, 2.95
};

vector<double> mean_values = {
    0.188361, 0.0987575, 0.0940236, 0.0548412, 0.0292177, 0.0162554, 0.00353505, -0.00739133,
    -0.0186999, -0.0271535, -0.0359679, -0.0425498, -0.0468213, -0.0511101, -0.0510728, -0.0552002,
    -0.0584904, -0.0582249, -0.0621789, -0.065081, -0.0672496, -0.0632703, -0.0684037, -0.0665262,
    -0.0702069, -0.0707249
};

vector<double> mean_m5sigma = {
    -9.48387, -7.0994, -5.73342, -5.30606, -5.14887, -5.03047, -4.94465, -4.86987,
    -4.82894, -4.80306, -4.79058, -4.79006, -4.75566, -4.75122, -4.72338, -4.73995,
    -4.69308, -4.69223, -4.67605, -4.66903, -4.66251, -4.63706, -4.61914, -4.60513,
    -4.63307, -4.57132
};

vector<double> mean_p5sigma = {
    9.86059, 7.29692, 5.92147, 5.41574, 5.20731, 5.06299, 4.95172, 4.85508,
    4.79154, 4.74875, 4.71865, 4.70496, 4.66202, 4.649, 4.62123, 4.62955,
    4.57609, 4.57578, 4.55169, 4.53887, 4.52801, 4.51052, 4.48234, 4.47208,
    4.49265, 4.42987
};


    // Create TGraphs for the bounds and mean
    TGraph* graph_m5s = new TGraph(p.size(), &p[0], &mean_m5sigma[0]);
    TGraph* graph_p5s = new TGraph(p.size(), &p[0], &mean_p5sigma[0]);
    TGraph* graph_mean = new TGraph(p.size(), &p[0], &mean_values[0]);

    // Define exponential decay functions with adjusted initial parameters
    TF1* fit_m5s = new TF1("fit_m5s", "[0] + [1]*exp(-[2]*x)", 0.4, 3.7);
    TF1* fit_p5s = new TF1("fit_p5s", "[0] + [1]*exp(-[2]*x)", 0.4, 3.7);
    fit_m5s->SetParameters(-4.9, -1,1); // Adjusted initial guess
    fit_p5s->SetParameters(4.8, 2, 1);
    fit_m5s->SetLineColor(kRed);
    fit_p5s->SetLineColor(kBlue);
    fit_m5s->SetLineWidth(1);
    fit_p5s->SetLineWidth(1);

    // Perform fits
    graph_m5s->Fit(fit_m5s, "RN");
    graph_p5s->Fit(fit_p5s, "RN");

    // Clone and draw
    TF1* fit_m5s_draw = (TF1*)fit_m5s->Clone("fit_m5s_draw");
    TF1* fit_p5s_draw = (TF1*)fit_p5s->Clone("fit_p5s_draw");
    fit_m5s_draw->Draw("SAME");
    fit_p5s_draw->Draw("SAME");

    // Draw data points
    graph_m5s->SetMarkerStyle(20); // Filled circles
    graph_m5s->SetMarkerColor(kRed);
    graph_m5s->SetMarkerSize(0.8);
    graph_m5s->Draw("P SAME");

    graph_p5s->SetMarkerStyle(20); // Filled circles
    graph_p5s->SetMarkerColor(kBlue);
    graph_p5s->SetMarkerSize(0.8);
    graph_p5s->Draw("P SAME");

    graph_mean->SetMarkerStyle(20); // Filled circles for mean
    graph_mean->SetMarkerColor(kBlack);
    graph_mean->SetMarkerSize(0.8);
    graph_mean->Draw("P SAME");

    // Add legend (moved slightly to avoid overlap)
    TLegend* legend = new TLegend(0.65, 0.65, 0.85, 0.85);
    legend->AddEntry(fit_m5s, "#mu - 5#sigma (fit)", "l");
    legend->AddEntry(fit_p5s, "#mu + 5#sigma (fit)", "l");
    legend->AddEntry(graph_m5s, "#mu - 5#sigma (data)", "p");
    legend->AddEntry(graph_p5s, "#mu + 5#sigma (data)", "p");
    legend->AddEntry(graph_mean, "#mu values", "p");
    legend->SetBorderSize(0);
    legend->SetTextSize(0.03); // Smaller text size
    legend->Draw();

    // Create formatted strings for fit equations
    auto formatEquation = [](TF1* fit, const string& prefix) {
        stringstream ss;
        ss << prefix << ": " << fixed << setprecision(3) << fit->GetParameter(0) << " + " 
           << showpos << setprecision(3) << fit->GetParameter(1) << "*exp(-" 
           << noshowpos << setprecision(3) << fit->GetParameter(2) << "*x)";
        return ss.str();
    };

    // Add fit equations to the plot
    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.025);
    tex->SetTextAlign(12);
    
    tex->DrawLatex(0.15, 0.85, formatEquation(fit_m5s, "#mu - 5#sigma").c_str());
    tex->DrawLatex(0.15, 0.82, formatEquation(fit_p5s, "#mu + 5#sigma").c_str());

    // Print fit parameters to console
    cout << "\nFit results for #mu - 5σ:" << endl;
    fit_m5s->Print();
    cout << "\nFit results for #mu + 5σ:" << endl;
    fit_p5s->Print();

    // Save to PDF
    canvas->Print("output/chi2pid_vs_p_2d_211.pdf");

    // Clean up
    delete graph_m5s;
    delete graph_p5s;
    delete graph_mean;
    delete fit_m5s;
    delete fit_p5s;
    delete legend;
    delete tex;
    delete canvas;
    delete fit_m5s_draw;
    delete fit_p5s_draw;
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