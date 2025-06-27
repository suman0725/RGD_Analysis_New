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
    TTree* tree_proton = dynamic_cast<TTree*>(file->Get("proton"));
    if (!tree_proton) {
        std::cerr << "Error: Could not find tree 'proton' in the file." << std::endl;
        file->Close();
        delete file;
        return 1;
    }

    // Variables to hold branch data
    Float_t delta_t_proton, momentum_proton;
    tree_proton->SetBranchAddress("delta_t", &delta_t_proton);
    tree_proton->SetBranchAddress("momentum", &momentum_proton);

    // Define momentum slices
    const Float_t momentum_min = 0.2;
    const Float_t momentum_max = 0.6;
    const Float_t slice_size = 0.05; // 50 MeV slices
    const int num_slices = static_cast<int>((momentum_max - momentum_min) / slice_size);

    // Create a PDF to save histograms
    TString pdf_file = "delta_t_slices_proton_new.pdf";
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
        Long64_t nentries = tree_proton->GetEntries();
        for (Long64_t j = 0; j < nentries; j++) {
            tree_proton->GetEntry(j);
            if (momentum_proton >= slice_start && momentum_proton < slice_end) {
                hist->Fill(delta_t_proton);
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
