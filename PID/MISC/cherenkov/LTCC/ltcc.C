#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include <cmath>

void ltcc() {
    // Open the ROOT file
    TFile* inFile = new TFile("charged_particles_LTCC_HTCC.root", "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Could not open the file!" << std::endl;
        return;
    }

    // Get the TTree
    TTree* tree = (TTree*)inFile->Get("charged_particle");
    if (!tree) {
        std::cerr << "Error: Could not find the tree 'charged_particle' in the file!" << std::endl;
        inFile->Close();
        return;
    }

    // Variables to read from the tree
    int detector;
    short status; 
    float px, py, pz, path, time, beta;

    // Set branch addresses
    tree->SetBranchAddress("detector", &detector);
    tree->SetBranchAddress("status", &status);
    tree->SetBranchAddress("path", &path);  // Path in cm
    tree->SetBranchAddress("time", &time);  // Time in ns

    
    const float c = 3.0e10;  // Speed of light in cm/s

    // Loop through the tree
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        // Check if the detector is LTCC
        if (detector == 16) {
            // Check if the particle is in the forward detector
            if (std::abs(status / 2000) == 1) {
                // Calculate beta using path and time (path in cm, time in ns)
                beta = path / (time * 1e-9 * c);  // Convert time from ns to seconds

                // Print the information
                std::cout << "Entry " << i << ": "
                          << "Beta = " << beta << "," << "path = " << path<< "," << "time= " << time << std::endl;
            }
        }
    }

    // Close the file
    inFile->Close();
    std::cout << "Processing complete!" << std::endl;
}
