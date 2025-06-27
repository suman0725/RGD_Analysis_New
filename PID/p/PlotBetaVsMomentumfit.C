#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLatex.h"
#include <TLegend.h>
#include "TF1.h"

void PlotBetaVsMomentumfit(const char* filename = "charged_particles.root", 
                            const char* treename = "charged_particle", 
                            const char* output = "beta_vs_p_posi_parti_cpp_fit.png", 
                            const char* rootfile = "beta_vs_p_posi_parti_cpp_fit.root") {
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
    Short_t status;

    // Set branch addresses
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("beta", &beta);
    tree->SetBranchAddress("charge", &charge);
    tree->SetBranchAddress("status", &status);
    tree->SetBranchAddress("startTime", &startTime);
    tree->SetBranchAddress("pid", &pid);
    tree->SetBranchAddress("detector", &detector);
    tree->SetBranchAddress("layer", &layer);

    // Create a 2D histogram for beta vs p
    TH2F *h_beta_vs_p = new TH2F("h_beta_vs_p", "Beta vs Momentum;Momentum (GeV);Beta",
                                100, 0, 10, 200, 0.5, 1.1);

    // Loop over the tree entries
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        
        if (std::abs(status / 2000) == 1 && startTime >= 0 && pid > 11 ) { 
            // Include only positively charged particles
            if (charge > 0) {
                // Calculate momentum magnitude (p)
                float p = std::sqrt(px * px + py * py + pz * pz);
                
                // Fill the histogram for all particles (p, beta)
                h_beta_vs_p->Fill(p, beta);
            }
        }
    }

    // Create functions for beta(p)
    double m_pi = 0.13957;  // Pion mass in GeV
    double m_k = 0.49368;   // Kaon mass in GeV
    double m_p = 0.93827;   // Proton mass in GeV

    TF1 *f_pi = new TF1("f_pi", "x / sqrt(x*x + [0]*[0])", 0, 10);  // Beta function for pion
    f_pi->SetParameter(0, m_pi);  // Set pion mass parameter

    TF1 *f_k = new TF1("f_k", "x / sqrt(x*x + [0]*[0])", 0, 10);  // Beta function for kaon
    f_k->SetParameter(0, m_k);  // Set kaon mass parameter

    TF1 *f_p = new TF1("f_p", "x / sqrt(x*x + [0]*[0])", 0, 10);  // Beta function for proton
    f_p->SetParameter(0, m_p);  // Set proton mass parameter

    // Open a new ROOT file to save the histogram and fit functions
    TFile *outputFile = new TFile(rootfile, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Cannot open output ROOT file!" << std::endl;
        file->Close();
        return;
    }

    // Save the histogram and the fit functions to the ROOT file
    h_beta_vs_p->Write();  // Write histogram to ROOT file
    f_pi->Write();         // Write pion beta function to ROOT file
    f_k->Write();          // Write kaon beta function to ROOT file
    f_p->Write();          // Write proton beta function to ROOT file

    // Draw the histogram
    TCanvas *canvas = new TCanvas("canvas", "Beta vs Momentum", 800, 600);
    h_beta_vs_p->SetStats(0);
    h_beta_vs_p->SetOption("COLZ");
    h_beta_vs_p->SetContour(99);
    canvas->SetLogz();
    h_beta_vs_p->Draw("COLZ");

    // Add the total entries text
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextColor(kBlack);
    TString entriesText = Form("Total Entries: %lld", (Long64_t)h_beta_vs_p->GetEntries());
    latex.DrawLatex(0.15, 0.85, entriesText);

    // Save the plot as a PNG image
    canvas->SaveAs(output);

  
}

