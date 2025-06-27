#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>
#include <iostream>
#include <chrono>

using namespace std;
using namespace ROOT;

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
    ROOT::RDataFrame df_pions("EB_pos_pion_assumed", file);
    gSystem->mkdir("output", kTRUE);

    // Define 2D histogram model
    ROOT::RDF::TH2DModel modelDtVsP(
        "dt_vs_p",
        "Delta t vs Momentum; Momentum (GeV/c); #Delta t (ns)",
        200, 0.0, 7.0,  // Momentum range: 1.0 to 7.0 GeV/c
        200, -2.0, 2.0  // Delta t range: -2.0 to 2.0 ns
    );

    // Create 2D histogram with correct column specification
    auto histo_dt_vs_p = df_pions.Histo2D(modelDtVsP, {"p"}, {"dt"});
    histo_dt_vs_p->SetStats(0);
    histo_dt_vs_p->Sumw2();

    // Create canvas and draw
    TCanvas* canvas = new TCanvas("canvas", "Delta t vs Momentum", 800, 600);
    canvas->SetLogz(); 
    histo_dt_vs_p->Draw("COLZ"); // Draw as color map
    // Save to PDF
    canvas->Print("output/dt_vs_p_2d.pdf");

    // Clean up
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