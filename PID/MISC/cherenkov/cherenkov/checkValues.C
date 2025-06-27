#include <TFile.h>
#include <TTree.h>
#include <iostream>
using namespace std; 
void checkValues() {
    // Open the ROOT file with the tree
    TFile* file = TFile::Open("../LTCC/charged_particles_LTCC_HTCC.root", "READ");
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
    int charge, detector;
    short status;
    float px, py,nphe, pz, beta;
    

    // Set branches
    tree->SetBranchAddress("charge", &charge);
    tree->SetBranchAddress("status", &status);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("beta", &beta);
    tree->SetBranchAddress("nphe", &nphe);
    tree->SetBranchAddress("detector", &detector);

    // Loop through the tree and check values for p and beta
    int n = 0; 
    int nEntries = tree->GetEntries();
    for (int i = 0; i < nEntries; ++i) {  // Loop over first 10 entries
        tree->GetEntry(i);
        if (i>=100) break; 
        if (charge > 0  ) {
        //if (charge > 0  && detector==15 || detector==16) {
            cout << "nphe = " << nphe << endl;
            n++; 
        }
        
    }
    cout << "Total number of nphe values: " << n << endl;
    // Close the ROOT file
    file->Close();
}
