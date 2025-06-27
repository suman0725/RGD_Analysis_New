#include <TFile.h>
#include <TTree.h>
#include <iostream>

void checkValues() {
    // Open the ROOT file with the tree
    TFile* file = TFile::Open("charged_particles.root", "READ");
    if (!file || file->IsOpen() == false) {
        std::cerr << "Failed to open file!" << std::endl;
        return;
    }
    
    TTree* tree = (TTree*)file->Get("charged_particle");
    if (!tree) {
        std::cerr << "Failed to get tree!" << std::endl;
        return;
    }

    // Define variables to read from the tree
    int charge;
    short status;
    float px, py, pz, beta;

    // Set branches
    tree->SetBranchAddress("charge", &charge);
    tree->SetBranchAddress("status", &status);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("beta", &beta);

    // Loop through the tree and check values for p and beta
    int nEntries = tree->GetEntries();
    for (int i = 0; i < nEntries && i < 10; ++i) {  // Loop over first 10 entries
        tree->GetEntry(i);
        if (i>=1000) break; 
        
        // Calculate momentum p
        float p = sqrt(px * px + py * py + pz * pz);
        
        // Print out values of p and beta
        std::cout << "Entry " << i << ": p = " << p << ", beta = " << beta << std::endl;
        
        // Check if p or beta has non-zero values
        if (p > 0) {
            std::cout << "Non-zero momentum found: " << p << std::endl;
        }
        if (beta > 0) {
            std::cout << "Non-zero beta found: " << beta << std::endl;
        }
    }

    // Close the ROOT file
    file->Close();
}
