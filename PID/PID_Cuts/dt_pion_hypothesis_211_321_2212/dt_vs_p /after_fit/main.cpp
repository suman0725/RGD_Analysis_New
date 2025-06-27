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
    ROOT::RDataFrame df_all("EB_all_pion_assumed", file);

    // Create 2D histogram model for dt vs p
    auto h2_dt_p = df_all.Histo2D(
        {"h2_dt_p", "#Delta t vs Momentum; Momentum (GeV/c); #Delta t (ns); Counts", 
         200, 0.0, 7.0,  // p range: 0 to 7 GeV
         200, -2.0, 2.0}, // dt range: -2 to 2 ns
        "p", "dt");

    // Initialize canvas
    TCanvas* canvas = new TCanvas("canvas", "#Delta t vs Momentum 2D Plot", 800, 600);
    
    // Draw 2D histogram
    h2_dt_p->Draw("COLZ");

    // Data points for fitting (limited to 0.4 to 2.7 GeV/c)
    vector<double> p = {0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400, 1.500, 1.600, 1.700, 1.800, 1.900, 2.000, 2.100, 2.200, 2.300, 2.400, 2.500, 2.600, 2.700};
    vector<double> mean_m3sigma = {-0.485, -0.339, -0.273, -0.256, -0.250, -0.246, -0.243, -0.239, -0.238, -0.237, -0.236, -0.235, -0.234, -0.232, -0.231, -0.231, -0.228, -0.228, -0.227, -0.226, -0.225, -0.224, -0.224, -0.223};
    vector<double> mean_p3sigma = {0.515, 0.350, 0.284, 0.262, 0.253, 0.248, 0.243, 0.238, 0.235, 0.233, 0.230, 0.229, 0.227, 0.225, 0.224, 0.222, 0.219, 0.220, 0.217, 0.217, 0.215, 0.214, 0.214, 0.213};

    // Calculate mean values
    vector<double> mean_values(p.size());
    for (size_t i = 0; i < p.size(); ++i) {
        mean_values[i] = (mean_m3sigma[i] + mean_p3sigma[i]) / 2.0;
    }

   // Create TGraphs for the bounds and mean
TGraph* graph_m3s = new TGraph(p.size(), &p[0], &mean_m3sigma[0]);
TGraph* graph_p3s = new TGraph(p.size(), &p[0], &mean_p3sigma[0]);
TGraph* graph_mean = new TGraph(p.size(), &p[0], &mean_values[0]);

// Define exponential decay functions [0] + [1]*exp(-[2]*x)
TF1* fit_m3s = new TF1("fit_m3s", "[0] + [1]*exp(-[2]*x)", 0.4, 2.7);
TF1* fit_p3s = new TF1("fit_p3s", "[0] + [1]*exp(-[2]*x)", 0.4, 2.7);

// Set initial parameters based on data trend
fit_m3s->SetParameters(-0.225, -0.3, 2.0);
fit_p3s->SetParameters(0.225, 0.3, 2.0);

// FIRST set all visual properties BEFORE fitting
fit_m3s->SetLineColor(kRed);
fit_p3s->SetLineColor(kBlue);
fit_m3s->SetLineWidth(1);
fit_p3s->SetLineWidth(1);
fit_m3s->SetRange(0, 7);
fit_p3s->SetRange(0, 7);

// Now perform the fits with "N" option to prevent immediate drawing
graph_m3s->Fit(fit_m3s, "RNQ");  // R=range, N=no draw, Q=quiet
graph_p3s->Fit(fit_p3s, "RNQ");

// Create CLONES of the fit functions to ensure independent drawing
TF1* fit_m3s_draw = (TF1*)fit_m3s->Clone("fit_m3s_draw");
TF1* fit_p3s_draw = (TF1*)fit_p3s->Clone("fit_p3s_draw");

// Set properties again on the clones (extra safety)
fit_m3s_draw->SetLineColor(kRed);
fit_p3s_draw->SetLineColor(kBlue);
fit_m3s_draw->SetLineWidth(1);
fit_p3s_draw->SetLineWidth(1);

// Draw the CLONED fit functions
fit_m3s_draw->Draw("SAME");
fit_p3s_draw->Draw("SAME");

// Don't forget to clean up the clones later
// (Add these to your cleanup section at the end)
// delete fit_m3s_draw;
// delete fit_p3s_draw;

    // Draw data points
    graph_m3s->SetMarkerStyle(20); // Filled circles
    graph_m3s->SetMarkerColor(kRed);
    graph_m3s->SetMarkerSize(0.8);
    graph_m3s->Draw("P SAME");

    graph_p3s->SetMarkerStyle(20); // Filled circles
    graph_p3s->SetMarkerColor(kBlue);
    graph_p3s->SetMarkerSize(0.8);
    graph_p3s->Draw("P SAME");

    graph_mean->SetMarkerStyle(20); // Filled circles for mean
    graph_mean->SetMarkerColor(kBlack);
    graph_mean->SetMarkerSize(0.8);
    graph_mean->Draw("P SAME");

    // Add legend (moved slightly to avoid overlap)
    TLegend* legend = new TLegend(0.65, 0.65, 0.85, 0.85);
    legend->AddEntry(fit_m3s, "#mu - 3#sigma (fit)", "l");
    legend->AddEntry(fit_p3s, "#mu + 3#sigma (fit)", "l");
    legend->AddEntry(graph_m3s, "#mu - 3#sigma (data)", "p");
    legend->AddEntry(graph_p3s, "#mu + 3#sigma (data)", "p");
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
    
    tex->DrawLatex(0.15, 0.85, formatEquation(fit_m3s, "#mu - 3#sigma").c_str());
    tex->DrawLatex(0.15, 0.82, formatEquation(fit_p3s, "#mu + 3#sigma").c_str());

    // Print fit parameters to console
    cout << "\nFit results for #mu - 3σ:" << endl;
    fit_m3s->Print();
    cout << "\nFit results for #mu + 3σ:" << endl;
    fit_p3s->Print();

    // Save to PDF
    canvas->Print("output/dt_vs_p_2d.pdf");

    // Clean up
    delete graph_m3s;
    delete graph_p3s;
    delete graph_mean;
    delete fit_m3s;
    delete fit_p3s;
    delete legend;
    delete tex;
    delete canvas;
    delete fit_m3s_draw;
delete fit_p3s_draw;
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