#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <iostream>

using namespace std;
void plotBeta() {
    // Open the ROOT file with the tree
    TFile* file = TFile::Open("charged_particles_LTCC_HTCC.root", "READ");
    TTree* tree = (TTree*)file->Get("charged_particle");

    // Define variables to read from the tree
    int pid, charge, detector;
    short status;
    float px, py, pz, beta;
    
    // Set branches
    tree->SetBranchAddress("pid", &pid);
    tree->SetBranchAddress("charge", &charge);
    tree->SetBranchAddress("status", &status);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("beta", &beta);
    tree->SetBranchAddress("detector", &detector);

    // Create a PDF to store the histograms
    TCanvas* canvas = new TCanvas("canvas", "Beta Distributions", 1200, 800);
    canvas->Print("beta_histograms_LTCC_HTCC.pdf[");

    // Define momentum ranges and corresponding histogram bins
    float pLow, pHigh;
    const float pStep = 0.25;
    const float pMax = 8.5;
    
    // Loop over different ranges of momentum (p)
    for (pLow = 2.0; pLow < pMax; pLow += pStep) {
        pHigh = pLow + pStep;
        
        // Create a histogram for the current momentum range
        TH1F* hist = new TH1F(Form("hist_beta_%.2f_%.2f", pLow, pHigh), 
                              Form("Beta Distribution for %.2f < p < %.2f", pLow, pHigh),
                              100, 0.9, 1.02);  // 100 bins for beta from 0 to 1

        // Loop through the tree entries
        int nEntries = tree->GetEntries();
        for (int i = 0; i < nEntries; i++) {
            tree->GetEntry(i);
            //if (i >= 1000) break; 

            // Apply cuts: status check (abs(status) / 2000 == 1) and charge > 0
            if (abs(status) / 2000 == 1 && charge > 0 &&  detector ==16 ) {
                // Calculate the momentum
                float p = sqrt(px * px + py * py + pz * pz);
                
                // Check if the momentum is in the desired range
                if (p >= pLow && p < pHigh) {
                    hist->Fill(beta);
                }
            }
        }

       

        // Plot the histogram for this momentum range
        hist->Draw();
        canvas->Print("beta_histograms_LTCC_HTCC.pdf");

        // Clean up the histogram after plotting
        delete hist;
    }

    // Close the PDF file and the ROOT file
    canvas->Print("beta_histograms_LTCC_HTCC.pdf]");
    file->Close();

    std::cout << "Histograms saved to beta_histograms.pdf" << std::endl;
}
