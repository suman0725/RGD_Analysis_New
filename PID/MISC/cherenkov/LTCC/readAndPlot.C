#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <cmath>
#include <chrono>  // For the timer

void readAndPlot() {
    // Start the timer
    auto start_time = std::chrono::high_resolution_clock::now();

    // Open the ROOT file and get the tree
    TFile* file = TFile::Open("charged_particles_LTCC_HTCC.root");
    if (!file || !file->IsOpen()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("charged_particle");
    if (!tree) {
        std::cerr << "Error retrieving tree!" << std::endl;
        return;
    }

    // Branch variables
    int pid, charge, detector, sector;
    short status;  // status as short
    float px, py, pz, nphe, beta;

    // Set up branches
    tree->SetBranchAddress("pid", &pid);
    tree->SetBranchAddress("charge", &charge);
    tree->SetBranchAddress("detector", &detector);
    tree->SetBranchAddress("sector", &sector);
    tree->SetBranchAddress("status", &status);  // status as short
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("nphe", &nphe);
    tree->SetBranchAddress("beta", &beta);

    // Create histograms
    TH1F* h_momentum_all = new TH1F("h_momentum_all", "Momentum Distribution for All +ve Charged Particles", 100, 0, 8);
    TH1F* h_momentum_LTCC = new TH1F("h_momentum_LTCC", "Momentum Distribution for +ve Charged Particles with LTCC Signal", 100, 0, 8);
    TH1F* h_nphe = new TH1F("h_nphe", "Number of Photoelectrons (nphe) for Particles with LTCC Signal", 100, 0, 30);
    TH2F* h_nphe_vs_momentum = new TH2F("h_nphe_vs_momentum", "nphe vs. Momentum for Particles with LTCC Signal", 100, 0, 8, 100, 0, 30);

    // Define momentum ranges for beta histograms
    double momentumRanges[][2] = {
        {2.5, 3.0}, {3.0, 3.5}, {3.5, 4.0}, {4.0, 4.5},
        {4.5, 5.0}, {5.0, 5.5}, {5.5, 6.0}, {6.0, 6.5},
        {6.5, 7.0}, {7.0, 7.5}, {7.5, 8.0}, {8.0, 8.5}};
    const int numRanges = sizeof(momentumRanges) / sizeof(momentumRanges[0]);

    TH1F* h_beta[numRanges];  // Array to store beta histograms
    for (int i = 0; i < numRanges; i++) {
        TString name = Form("h_beta_p%.1f_%.1f", momentumRanges[i][0], momentumRanges[i][1]);
        TString title = Form("Beta Distribution for %.1f < p < %.1f", momentumRanges[i][0], momentumRanges[i][1]);
        h_beta[i] = new TH1F(name, title, 100, 0.8, 1.2);
    }

    // Loop through the tree
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        // Calculate momentum
        float momentum = sqrt(px * px + py * py + pz * pz);

        // Fill histogram for all positive-charged particles
        if (charge > 0 && (abs(static_cast<int>(status)) / 2000) == 1) {
            h_momentum_all->Fill(momentum);

            // Fill histogram only for LTCC signal (detector == 16)
            if (detector == 16) {
                h_momentum_LTCC->Fill(momentum);
                h_nphe->Fill(nphe);
                h_nphe_vs_momentum->Fill(momentum, nphe);

                // Fill beta histograms based on momentum ranges
                for (int j = 0; j < numRanges; j++) {
                    if (momentum >= momentumRanges[j][0] && momentum < momentumRanges[j][1]) {
                        h_beta[j]->Fill(beta);
                        break;
                    }
                }
            }
        }
    }

    // Canvas 1: Momentum comparison
    TCanvas* c1 = new TCanvas("c1", "Momentum Distribution Before and After LTCC Signal", 800, 600);
    h_momentum_all->SetLineColor(kBlack);
    h_momentum_all->Draw();
    h_momentum_LTCC->SetLineColor(kRed);
    h_momentum_LTCC->Draw("SAME");
    TLegend* legend1 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend1->AddEntry(h_momentum_all, "All +ve Charged Particles (Black)", "l");
    legend1->AddEntry(h_momentum_LTCC, "With LTCC Signal (Red)", "l");
    legend1->Draw();
    c1->SaveAs("momentum_LTCC_comparison.png");

    // Canvas 2: nphe histogram
    TCanvas* c2 = new TCanvas("c2", "nphe Distribution for Particles with LTCC Signal", 800, 600);
    h_nphe->SetLineColor(kBlue);
    h_nphe->Draw();
    c2->SaveAs("nphe_distribution.png");

    // Canvas 3: nphe vs. momentum
    TCanvas* c3 = new TCanvas("c3", "nphe vs. Momentum for Particles with LTCC Signal", 800, 600);
    h_nphe_vs_momentum->SetMarkerStyle(20);
    h_nphe_vs_momentum->Draw("COLZ");
    c3->SaveAs("nphe_vs_momentum.png");

    // Canvas 4: Beta histograms
    TCanvas* c4 = new TCanvas("c4", "Beta Distributions for Different Momentum Ranges", 1200, 800);
    c4->Divide(4, 3);  // Divide canvas into subplots

    for (int i = 0; i < numRanges; i++) {
        c4->cd(i + 1);
        h_beta[i]->SetLineColor(i + 1);  // Assign unique colors
        h_beta[i]->Draw();
    }
    c4->SaveAs("beta_distributions.png");

    // Stop the timer
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = end_time - start_time;

    // Print elapsed time
    std::cout << "Script execution time: " << elapsed_time.count() << " seconds." << std::endl;

    // Close the file
    file->Close();
}
