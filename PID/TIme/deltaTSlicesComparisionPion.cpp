#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <string>
#include <TF1.h>
#include "TPaveText.h"
#include <TLegend.h>

int main() {
    // Open the ROOT files
    TFile* file1 = TFile::Open("particles_data_pid.root");
    TFile* file2 = TFile::Open("particles_data.root");
    if (!file1 || file1->IsZombie() || !file2 || file2->IsZombie()) {
        std::cerr << "Error: Could not open one of the files." << std::endl;
        return 1;
    }

    // Access the trees
    TTree* tree1 = dynamic_cast<TTree*>(file1->Get("pion"));
    TTree* tree2 = dynamic_cast<TTree*>(file2->Get("pion"));
    if (!tree1 || !tree2) {
        std::cerr << "Error: Could not find tree in one of the files." << std::endl;
        file1->Close();
        file2->Close();
        delete file1;
        delete file2;
        return 1;
    }

    // Variables to hold branch data
    Float_t delta_t, momentum;

    // Set branch addresses for both files (with and without PID)
    tree1->SetBranchAddress("delta_t", &delta_t);
    tree1->SetBranchAddress("momentum", &momentum);
    
    tree2->SetBranchAddress("delta_t", &delta_t);
    tree2->SetBranchAddress("momentum", &momentum);

    // Define momentum slices
    const Float_t momentum_min = 0.2;
    const Float_t momentum_max = 7.0;
    const Float_t slice_size = 0.05;
    const int num_slices = static_cast<int>((momentum_max - momentum_min) / slice_size);

    // Create a PDF to save histograms
    TString pdf_file = "delta_t_comparison_legend_pion.pdf";
    TCanvas* canvas = new TCanvas("canvas", "Momentum Slices", 800, 600);

    // Loop over slices and process both files
    for (int i = 0; i < num_slices; i++) {
        Float_t slice_start = momentum_min + i * slice_size;
        Float_t slice_end = slice_start + slice_size;

        // Create histograms for both datasets
        TString hist_name1 = Form("DeltaT_slice_with_pid_%d", i);
        TString hist_title1 = Form("Momentum range [%.2f, %.2f] GeV (With PID)", slice_start, slice_end);
        TH1F* hist_with_pid = new TH1F(hist_name1, hist_title1, 100, -3, 3);
        hist_with_pid->GetXaxis()->SetTitle("#DeltaT (ns)");
        hist_with_pid->GetYaxis()->SetTitle("Counts");

        TString hist_name2 = Form("DeltaT_slice_without_pid_%d", i);
        TString hist_title2 = Form("Momentum range [%.2f, %.2f] GeV (Without PID)", slice_start, slice_end);
        TH1F* hist_without_pid = new TH1F(hist_name2, hist_title2, 100, -3, 3);
        hist_without_pid->GetXaxis()->SetTitle("#DeltaT (ns)");
        hist_without_pid->GetYaxis()->SetTitle("Counts");

        // Loop over tree1 (with PID) and fill the histogram
        Long64_t nentries1 = tree1->GetEntries();
        for (Long64_t j = 0; j < nentries1; j++) {
            tree1->GetEntry(j);
            if (momentum >= slice_start && momentum < slice_end) {
                hist_with_pid->Fill(delta_t);
            }
        }

        // Loop over tree2 (without PID) and fill the histogram
        Long64_t nentries2 = tree2->GetEntries();
        for (Long64_t j = 0; j < nentries2; j++) {
            tree2->GetEntry(j);
            if (momentum >= slice_start && momentum < slice_end) {
                hist_without_pid->Fill(delta_t);
            }
        }

        // Draw the histograms
        
        
        hist_with_pid->SetLineColor(kRed);
        hist_with_pid->SetFillColorAlpha(kGreen, 0.5); // Semi-transparent red
        hist_with_pid->SetFillStyle(3001);  // Solid fill

        hist_without_pid->SetLineColor(kBlack);
        hist_without_pid->SetFillColorAlpha(kBlack, 0.5); // Semi-transparent black
        hist_without_pid->SetFillStyle(3004);  // Hatch pattern
        hist_without_pid->SetStats(0); 

        // Draw the histograms
        hist_without_pid->Draw("HIST");   // Draw black histogram first
        hist_with_pid->Draw("HIST SAME"); // Draw red histogram on top

        /* // Fit the histogram with a Gaussian in the range [-0.5, 0.5]
        TF1* gaus_fit = new TF1("gaus_fit", "gaus", -0.2, 0.2);
        gaus_fit->SetParameters(hist_with_pid->GetMaximum(), 0, 0.1);
        gaus_fit->SetLineColor(kGreen);
        hist_with_pid->Fit(gaus_fit, "R"); */

        // Add a legend
        TLegend* legend = new TLegend(0.65, 0.75, 0.9, 0.9); // Adjust position if needed
        legend->SetBorderSize(0);  // No border
        legend->SetFillStyle(0);   // Transparent background
        legend->AddEntry(hist_without_pid, Form("Without PID (Entries: %d)", (int)hist_without_pid->GetEntries()), "f");
        legend->AddEntry(hist_with_pid, Form("With PID (Entries: %d)", (int)hist_with_pid->GetEntries()), "f");
        legend->Draw();


        // Save the histogram to the PDF
        if (i == 0) {
            canvas->Print(pdf_file + "(", "pdf");
        } else if (i == num_slices - 1) {
            canvas->Print(pdf_file + ")", "pdf");
        } else {
            canvas->Print(pdf_file, "pdf");
        }

        // Clean up
        delete hist_with_pid;
        delete hist_without_pid;
    }

    // Final cleanup
    canvas->Close();
    delete canvas;
    file1->Close();
    file2->Close();
    delete file1;
    delete file2;

    std::cout << "Histograms saved to " << pdf_file << std::endl;
    return 0;
}
