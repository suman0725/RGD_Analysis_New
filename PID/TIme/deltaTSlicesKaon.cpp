#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <string>

int main() {
    // Open the ROOT file
    TFile* file = TFile::Open("particles_data_pid.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file." << std::endl;
        return 1;
    }

    // Access the tree
    TTree* tree_kaon = dynamic_cast<TTree*>(file->Get("kaon"));
    if (!tree_kaon) {
        std::cerr << "Error: Could not find tree 'kaon' in the file." << std::endl;
        file->Close();
        delete file;
        return 1;
    }

    // Variables to hold branch data
    Float_t delta_t_kaon, momentum_kaon;
    tree_kaon->SetBranchAddress("delta_t", &delta_t_kaon);
    tree_kaon->SetBranchAddress("momentum", &momentum_kaon);

    // Define momentum slices
    //const Float_t momentum_min = 0.2;
    const Float_t momentum_min = 1.3;
    const Float_t momentum_max = 7;
    //const Float_t momentum_max = 0.6;
    //const Float_t slice_size = 0.01; // 50 MeV slices
    const Float_t slice_size = 0.05;
    const int num_slices = static_cast<int>((momentum_max - momentum_min) / slice_size);

    // Create a PDF to save histograms
    TString pdf_file = "delta_t_slices_kaon.pdf";
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
        Long64_t nentries = tree_kaon->GetEntries();
        for (Long64_t j = 0; j < nentries; j++) {
            tree_kaon->GetEntry(j);
            if (momentum_kaon >= slice_start && momentum_kaon < slice_end) {
                hist->Fill(delta_t_kaon);
            }
        }

        // Draw the histogram
        hist->Draw();

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