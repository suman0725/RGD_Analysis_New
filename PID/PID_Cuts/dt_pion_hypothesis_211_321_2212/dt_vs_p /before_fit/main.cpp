#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>
#include <iostream>
#include <chrono>

using namespace std;
using namespace ROOT;

int main() {
    // Start total program timer
    auto start = chrono::high_resolution_clock::now();

    // Enable multi-threading
    ROOT::EnableImplicitMT(4);
    cout << "Multi-threading enabled for RDataFrame processing." << endl;

    gStyle->SetOptLogz(1);
    // Set style to remove stat box
    gStyle->SetOptStat(0);

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
        {"h2_dt_p", "Delta t vs Momentum; Momentum (GeV/c); #Delta t (ns); Counts", 
         200, 0.0, 7.0,  // p range: 0 to 7 GeV
         00, -2.0, 2.0}, // dt range: -2 to 2 ns
        "p", "dt");

    // Initialize canvas
    TCanvas* canvas = new TCanvas("canvas", "Delta t vs Momentum 2D Plot", 800, 600);
    
    // Draw 2D histogram
    h2_dt_p->Draw("COLZ"); // Color map option

    // Save to PDF (no return value check due to void return in some ROOT versions)
    canvas->Print("output/dt_vs_p_2d.pdf"); // Note: Success depends on directory permissions and path

    // Clean up
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