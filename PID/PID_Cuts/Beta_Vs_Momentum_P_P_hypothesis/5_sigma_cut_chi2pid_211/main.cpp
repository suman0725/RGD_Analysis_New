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
    ROOT::RDataFrame df_pion("EB_pid_pions", file);

    // Create 2D histogram model for dt vs p
    auto h2_dt_p = df_pion.Histo2D(
        {"h2_dt_p", "chi2pid vs Momentum; Momentum (GeV/c); chi2pid; Counts", 
         200, 0.0, 7.0,  
         200, -15.0, 15.0}, 
        "p", "orig_chi2pid");

    // Initialize canvas
    TCanvas* canvas = new TCanvas("canvas", "chi2pid vs Momentum 2D Plot", 800, 600);
    
    // Draw 2D histogram
    h2_dt_p->Draw("COLZ");

    vector<double> p = {0.4, 0.7, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 3.4, 3.7};
vector<double> mean_values = {0.0922588, 0.0301706, -0.00706712, -0.0344082, -0.0490577, -0.0573817, -0.0651223, -0.0674492, -0.0732826, -0.0776369, -0.0864284, -0.0799175};
vector<double> mean_m3sigma = {-6.21599, -5.19239, -4.90447, -4.81546, -4.7656, -4.72811, -4.68587, -4.64245, -4.61782, -4.55872, -4.47899, -4.46576};
vector<double> mean_p3sigma = {6.4005, 5.25273, 4.89034, 4.74665, 4.66749, 4.61335, 4.55563, 4.50755, 4.47126, 4.40344, 4.30613, 4.30592};


// Create TGraphs for the bounds and mean
TGraph* graph_m3s = new TGraph(p.size(), &p[0], &mean_m3sigma[0]);
TGraph* graph_p3s = new TGraph(p.size(), &p[0], &mean_p3sigma[0]);
TGraph* graph_mean = new TGraph(p.size(), &p[0], &mean_values[0]);

// Define exponential decay functions with adjusted initial parameters
TF1* fit_m3s = new TF1("fit_m3s", "[0] + [1]*exp(-[2]*x)", 0.4, 3.7);
TF1* fit_p3s = new TF1("fit_p3s", "[0] + [1]*exp(-[2]*x)", 0.4, 3.7);
fit_m3s->SetParameters(-4.5, -2.0, 0.5); // Adjusted based on data trend
fit_p3s->SetParameters(4.5, 2.0, 0.5);
fit_m3s->SetLineColor(kRed);
fit_p3s->SetLineColor(kBlue);
fit_m3s->SetLineWidth(1);
fit_p3s->SetLineWidth(1);

// Perform fits and check quality
graph_m3s->Fit(fit_m3s, "RN");
graph_p3s->Fit(fit_p3s, "RN");

// Clone and draw
TF1* fit_m3s_draw = (TF1*)fit_m3s->Clone("fit_m3s_draw");
TF1* fit_p3s_draw = (TF1*)fit_p3s->Clone("fit_p3s_draw");
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
    legend->AddEntry(fit_m3s, "#mu - 5#sigma (fit)", "l");
    legend->AddEntry(fit_p3s, "#mu + 5#sigma (fit)", "l");
    legend->AddEntry(graph_m3s, "#mu - 5#sigma (data)", "p");
    legend->AddEntry(graph_p3s, "#mu + 5#sigma (data)", "p");
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
    
    tex->DrawLatex(0.15, 0.85, formatEquation(fit_m3s, "#mu - 5#sigma").c_str());
    tex->DrawLatex(0.15, 0.82, formatEquation(fit_p3s, "#mu + 5#sigma").c_str());

    // Print fit parameters to console
    cout << "\nFit results for #mu - 5σ:" << endl;
    fit_m3s->Print();
    cout << "\nFit results for #mu + 5σ:" << endl;
    fit_p3s->Print();

    // Save to PDF
    canvas->Print("output/chi2pid_vs_p_2d_211.pdf");

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