#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <string>
#include <TF1.h>
#include "TPaveText.h" 

int main() {
    // Open the ROOT file
    TFile* file = TFile::Open("particles_data_pid.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file." << std::endl;
        return 1;
    }

    // Access the tree
    TTree* tree_pion = dynamic_cast<TTree*>(file->Get("pion"));
    if (!tree_pion) {
        std::cerr << "Error: Could not find tree 'kaon' in the file." << std::endl;
        file->Close();
        delete file;
        return 1;
    }

    // Variables to hold branch data
    Float_t delta_t_pion, momentum_pion;
    tree_pion->SetBranchAddress("delta_t", &delta_t_pion);
    tree_pion->SetBranchAddress("momentum", &momentum_pion);

    // Define momentum slices
    //const Float_t momentum_min = 0.2;
    const Float_t momentum_min = 0.2;
    const Float_t momentum_max = 7;
    //const Float_t momentum_max = 0.6;
    //const Float_t slice_size = 0.01; // 50 MeV slices
    const Float_t slice_size = 0.05;
    const int num_slices = static_cast<int>((momentum_max - momentum_min) / slice_size);

    // Create a PDF to save histograms
    TString pdf_file = "delta_t_slices_pion_all.pdf";
    TCanvas* canvas = new TCanvas("canvas", "Momentum Slices", 800, 600);

    for (int i = 0; i < num_slices; i++) {
        Float_t slice_start = momentum_min + i * slice_size;
        Float_t slice_end = slice_start + slice_size;

        // Create a histogram for the current slice
        TString hist_name = Form("DeltaT_slice_%d", i);
        TString hist_title = Form("Momentum range [%.2f, %.2f] GeV", slice_start, slice_end);
        TH1F* hist = new TH1F(hist_name, hist_title, 100, -3, 3); // Adjust binning/range as needed
        hist->GetXaxis()->SetTitle("#DeltaT (ns)");
        hist->GetYaxis()->SetTitle("Counts");

        // Loop over the tree and fill the histogram
        Long64_t nentries = tree_pion->GetEntries();
        for (Long64_t j = 0; j < nentries; j++) {
            tree_pion->GetEntry(j);
            if (momentum_pion >= slice_start && momentum_pion < slice_end) {
                hist->Fill(delta_t_pion);
            }
        }


        // Create a canvas and draw the histogram
        canvas = new TCanvas(Form("canvas_%d", i), Form("Canvas %d", i), 800, 600);

         // Draw the histogram
        hist->Draw(); 

        // Fit the histogram with a Gaussian in the range [-0.5, 0.5]
        TF1* gaus_fit = new TF1("gaus_fit", "gaus", -0.2, 0.2);
        gaus_fit->SetParameters(hist->GetMaximum(), 0, 0.1); // Initial parameters: amplitude, mean, sigma
        gaus_fit->SetLineColor(kRed); // Set fit line color to red
        hist->Fit(gaus_fit, "R"); // "R" specifies the range

        // Get fit parameters
        Double_t mean = gaus_fit->GetParameter(1);  // Mean
        Double_t sigma = gaus_fit->GetParameter(2); // Standard deviation
        gaus_fit->Draw("same");

        // Add text box with mean and sigma in the top-left corner
        TPaveText* stats = new TPaveText(0.1, 0.8, 0.4, 0.9, "NDC"); // Top-left position in normalized coordinates
        stats->SetFillColor(0);
        stats->SetTextSize(0.04);
        stats->AddText(Form("Mean = %.3f", mean));
        stats->AddText(Form("SD = %.3f", sigma));
        stats->Draw();       

        // Save the histogram to the PDF
        if (i == 0) {
            canvas->Print(pdf_file + "(", "pdf"); // Open PDF
        } else if (i == num_slices - 1) {
            canvas->Print(pdf_file + ")", "pdf"); // Close PDF
        } else {
            canvas->Print(pdf_file, "pdf"); // Append
        }

        // Clean up
        delete hist;
    }

    // Final cleanup
    canvas->Close();
    delete canvas;
    file->Close();
    delete file;

    std::cout << "Histograms saved to " << pdf_file << std::endl;
    return 0;
}