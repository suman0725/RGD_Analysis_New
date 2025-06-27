#include <iostream>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

void PlotChi2PidVsMomentum(const char* filename = "charged_particles.root", 
                           const char* treename = "charged_particle", 
                           const char* output1 ="chi2pid_posit_pions_cpp.png") {
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
    float chi2pid, startTime;
    int charge, pid;
    Short_t status;

    // Set branch addresses
    tree->SetBranchAddress("chi2pid", &chi2pid);
    tree->SetBranchAddress("charge", &charge);
    tree->SetBranchAddress("status", &status);
    tree->SetBranchAddress("startTime", &startTime);
    tree->SetBranchAddress("pid", &pid);

    // Create histogram
    TH1F *h_chi2pid = new TH1F("h_chi2pid", "Chi2Pid Distribution for π⁺;Chi2Pid;Counts", 200, -15, 15);

    // Loop over tree entries (limit to 10,000 events)
    Long64_t nEntries = tree->GetEntries();
    Long64_t maxEvents = std::min(nEntries, (Long64_t)1000000);

    for (Long64_t i = 0; i < maxEvents; i++) {
        tree->GetEntry(i);
        if (charge <= 0) continue;  // Select only positively charged particles
        if (pid != 211) continue;
        if (abs(status/2000) != 1) continue;   // Select only π⁺ (positive pions)
        h_chi2pid->Fill(chi2pid);
    }

    // Save histogram to ROOT file
    TFile *outputFile = new TFile("particles_data_chi2pid_positive_pion.root", "RECREATE");
    h_chi2pid->Write();

    // Draw and save histogram
    TCanvas *canvas1 = new TCanvas("canvas1", "Chi2Pid Histogram for π⁺", 800, 600);
    h_chi2pid->SetStats(0);
    h_chi2pid->Draw();
    canvas1->SaveAs(output1);

    std::cout << "Plot saved as " << output1 << std::endl;

    // Cleanup
    delete outputFile;
    delete file;
}
