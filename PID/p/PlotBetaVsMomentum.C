#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLatex.h"

void PlotBetaVsMomentum(const char* filename = "charged_particles.root", 
                        const char* treename = "charged_particle", 
                        const char* output = "beta_vs_p_posi_parti_cpp_ftof1a.png") {
    // Open the ROOT file
    TFile *file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file!" << std::endl;
        return;
    }

    // Access the TTree
    TTree *tree = (TTree*)file->Get(treename);
    if (!tree) {
        std::cerr << "Error: Cannot find tree!" << std::endl;
        file->Close();
        return;
    }

    // Variables to hold branch data
    float px, py, pz, beta, startTime;
    int charge, pid, detector, layer;
    Short_t status; // Corrected type for status

    // Set branch a

    // Set branch addresses
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("beta", &beta);
    tree->SetBranchAddress("charge", &charge);
    tree->SetBranchAddress("status", &status);
    tree->SetBranchAddress("startTime",&startTime);
    tree->SetBranchAddress("pid", &pid);
    tree->SetBranchAddress("detector",&detector);
    tree->SetBranchAddress("layer",&layer);


    // Create a 2D histogram for beta vs p
    TH2F *h_beta_vs_p = new TH2F("h_beta_vs_p", "Beta vs Momentum;Momentum (GeV);Beta",
                                 100, 0, 10, 200, 0.5, 1.1);

    // Loop over the tree entries
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        if (i >= 1000) break; 
        tree->GetEntry(i);
         if (std::abs(status / 2000) == 1 && startTime >= 0 && pid > 11 ) { 
            // Include only positively charged particles
            if (charge > 0) {
                // Calculate momentum magnitude (p)
                float p = std::sqrt(px * px + py * py + pz * pz);

                // Fill the histogram
                h_beta_vs_p->Fill(p, beta);
            }
        }
}

    // Draw the histogram
    TCanvas *canvas = new TCanvas("canvas", "Beta vs Momentum", 800, 600);
    const Int_t NRGBs = 7; // Number of colors in the gradient
    const Int_t NCont = 99; // Number of color levels
    Double_t stops[NRGBs] = {0.00, 0.17, 0.34, 0.51, 0.68, 0.85, 1.00}; // Gradient stops
    Double_t red[NRGBs]   = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0}; // Red component
    Double_t green[NRGBs] = {0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0}; // Green component
    Double_t blue[NRGBs]  = {0.5, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0}; // Blue component
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    h_beta_vs_p->SetStats(0); // Hide stats box
    h_beta_vs_p->SetOption("COLZ"); // Color map
    h_beta_vs_p->SetContour(99); // Smooth color gradient
    canvas->SetLogz();
    h_beta_vs_p->Draw("COLZ");

    // Add the total entries text
    TLatex latex;
    latex.SetNDC(); // Use normalized device coordinates
    latex.SetTextSize(0.04);
    latex.SetTextColor(kBlack);
    TString entriesText = Form("Total Entries: %lld", (Long64_t)h_beta_vs_p->GetEntries());
    latex.DrawLatex(0.15, 0.85, entriesText); // Position (x, y) in normalized coordinates

    // Save the plot
    canvas->SaveAs(output);

    // Clean up
    delete canvas;
    delete h_beta_vs_p;
    file->Close();
    delete file;

    std::cout << "Plot saved as " << output << std::endl;
}
