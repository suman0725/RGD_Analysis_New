#include <ROOT/RDataFrame.hxx>
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>

    void processFile() {

    const std::string& filename = "charged_particles.root";
    // Open the ROOT file
    TFile file(filename.c_str(), "READ");
    // Access the TTree
    TTree* tree = (TTree*)file.Get("charged_particle");

    // Create a 2D histogram for chi2pid vs p
    TH2F* h_chi2pid_vs_p = new TH2F("h_chi2pid_vs_p", "Chi2Pid vs Momentum;Momentum (GeV/c);Chi2Pid", 
                                    100, 0, 10, 100, -5, 5);

    // Set up the branches (assuming you know the names of the branches)
    Float_t px, py, pz, chi2pid;
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("chi2pid", &chi2pid);

    // Loop over the tree entries
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        // Calculate momentum magnitude (p)
        Float_t p = TMath::Sqrt(px*px + py*py + pz*pz);

        // Fill the histogram
        h_chi2pid_vs_p->Fill(p, chi2pid);
    }

    // Draw the histogram
    TCanvas* canvas = new TCanvas("canvas", "Chi2Pid vs Momentum", 800, 600);
    h_chi2pid_vs_p->Draw("COLZ");

    // Save the plot
    canvas->SaveAs("chi2pid_vs_p_posicharpart_cpp.png");

    // Close the ROOT file
    file.Close();
}


