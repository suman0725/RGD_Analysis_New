#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <iostream>
#include <TStopwatch.h>

using namespace std;
void plotBeta() {

    TStopwatch timer; // Start timing
    timer.Start();

    // Open the ROOT file with the tree
    TFile* file = TFile::Open("charged_particles.root", "READ");
    TTree* tree = (TTree*)file->Get("charged_particle");

    // Define variables to read from the tree
    int pid, charge;
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

    // Create a PDF to store the histograms
    TCanvas* canvas = new TCanvas("canvas", "Beta Distributions", 1200, 800);
    canvas->Print("beta_histograms.pdf[");

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
                              200, 0.88, 1.02);  // 100 bins for beta from 0 to 1

        // Loop through the tree entries
        int nEntries = tree->GetEntries();
        for (int i = 0; i < nEntries; i++) {
            tree->GetEntry(i);
            //if (i >= 1000) break; 

            // Apply cuts: status check (abs(status) / 2000 == 1) and charge > 0
            if (abs(status) / 2000 == 1 && charge > 0) {
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
        canvas->Print("beta_histograms.pdf");

        // Clean up the histogram after plotting
        delete hist;
    }

    // Close the PDF file and the ROOT file
    canvas->Print("beta_histograms.pdf]");
    file->Close();

    // Stop the timer and display the elapsed time
    timer.Stop();
    cout << "Histograms saved to beta_histograms.pdf" << endl;
    cout << "Total time for process: " << timer.RealTime() << " seconds (real time), "
         << timer.CpuTime() << " seconds (CPU time)" << endl;
}
