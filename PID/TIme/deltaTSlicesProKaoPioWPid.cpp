#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>

int main() {
    TFile* file = TFile::Open("particles_data_pid.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file." << std::endl;
        return 1;
    }

    TTree* tree_pion = dynamic_cast<TTree*>(file->Get("pion"));
    TTree* tree_kaon = dynamic_cast<TTree*>(file->Get("kaon"));
    TTree* tree_proton = dynamic_cast<TTree*>(file->Get("proton"));

    if (!tree_pion || !tree_kaon || !tree_proton) {
        std::cerr << "Error: Could not find one of the trees." << std::endl;
        file->Close();
        return 1;
    }

    Float_t delta_t_pion, momentum_pion;
    Float_t delta_t_kaon, momentum_kaon;
    Float_t delta_t_proton, momentum_proton;

    tree_pion->SetBranchAddress("delta_t", &delta_t_pion);
    tree_pion->SetBranchAddress("momentum", &momentum_pion);
    tree_kaon->SetBranchAddress("delta_t", &delta_t_kaon);
    tree_kaon->SetBranchAddress("momentum", &momentum_kaon);
    tree_proton->SetBranchAddress("delta_t", &delta_t_proton);
    tree_proton->SetBranchAddress("momentum", &momentum_proton);

    const Float_t momentum_min = 0.2;
    const Float_t momentum_max = 7;
    const Float_t slice_size = 0.05;
    const int num_slices = static_cast<int>((momentum_max - momentum_min) / slice_size);

    TString pdf_file = "delta_t_slices_all_particles_pid.pdf";
    TCanvas* canvas = new TCanvas("canvas", "Momentum Slices", 800, 600);

    for (int i = 0; i < num_slices; i++) {
        Float_t slice_start = momentum_min + i * slice_size;
        Float_t slice_end = slice_start + slice_size;

        TString hist_title = Form("Momentum range [%.2f, %.2f] GeV", slice_start, slice_end);
        TH1F* hist_pion = new TH1F("hist_pion", hist_title, 100, -3, 3);
        TH1F* hist_kaon = new TH1F("hist_kaon", hist_title, 100, -3, 3);
        TH1F* hist_proton = new TH1F("hist_proton", hist_title, 100, -3, 3);

        hist_pion->GetXaxis()->SetTitle("#DeltaT (ns)");
        hist_pion->GetYaxis()->SetTitle("Counts");
        hist_pion->SetLineColor(kRed);
        hist_pion->SetFillColorAlpha(kRed, 0.35);

        hist_kaon->GetXaxis()->SetTitle("#DeltaT (ns)");
        hist_kaon->GetYaxis()->SetTitle("Counts");
        hist_kaon->SetLineColor(kBlue);
        hist_kaon->SetFillColorAlpha(kBlue, 0.35);

        hist_proton->GetXaxis()->SetTitle("#DeltaT (ns)");
        hist_proton->GetYaxis()->SetTitle("Counts");
        hist_proton->SetLineColor(kGreen);
        hist_proton->SetFillColorAlpha(kGreen, 0.35);

        Long64_t nentries_pion = tree_pion->GetEntries();
        Long64_t nentries_kaon = tree_kaon->GetEntries();
        Long64_t nentries_proton = tree_proton->GetEntries(); 


        // Process only the first 10,000 entries for quick checking
        /* Long64_t nentries_pion = std::min(tree_pion->GetEntries(), 10000LL);
        Long64_t nentries_kaon = std::min(tree_kaon->GetEntries(), 10000LL);
        Long64_t nentries_proton = std::min(tree_proton->GetEntries(), 10000LL); */

        for (Long64_t j = 0; j < nentries_pion; j++) {
            tree_pion->GetEntry(j);
            if (momentum_pion >= slice_start && momentum_pion < slice_end) {
                hist_pion->Fill(delta_t_pion);
            }
        }

        for (Long64_t j = 0; j < nentries_kaon; j++) {
            tree_kaon->GetEntry(j);
            if (momentum_kaon >= slice_start && momentum_kaon < slice_end) {
                hist_kaon->Fill(delta_t_kaon);
            }
        }

        for (Long64_t j = 0; j < nentries_proton; j++) {
            tree_proton->GetEntry(j);
            if (momentum_proton >= slice_start && momentum_proton < slice_end) {
                hist_proton->Fill(delta_t_proton);
            }
        }

        // Find the maximum bin content across all histograms
        float max_pion = hist_pion->GetMaximum();
        float max_kaon = hist_kaon->GetMaximum();
        float max_proton = hist_proton->GetMaximum();
        float max_all = std::max({max_pion, max_kaon, max_proton});

        // Set the Y-axis range for all histograms
        hist_pion->SetMaximum(max_all * 1.2); // Add 20% padding
        hist_kaon->SetMaximum(max_all * 1.2);
        hist_proton->SetMaximum(max_all * 1.2);

        canvas->Clear(); // Clear the canvas before drawing new histograms

        hist_pion->Draw("HIST");
        hist_kaon->Draw("HIST SAME");
        hist_proton->Draw("HIST SAME");

        TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(hist_pion, "Pion", "f");
        legend->AddEntry(hist_kaon, "Kaon", "f");
        legend->AddEntry(hist_proton, "Proton", "f");
        legend->Draw();

        gPad->Update(); // Force update

        if (i == 0) {
            canvas->Print(pdf_file + "(", "pdf");
        } else if (i == num_slices - 1) {
            canvas->Print(pdf_file + ")", "pdf");
        } else {
            canvas->Print(pdf_file, "pdf");
        }

        delete hist_pion;
        delete hist_kaon;
        delete hist_proton;
    }

    canvas->Close();
    delete canvas;
    file->Close();
    delete file;

    std::cout << "Histograms saved to " << pdf_file << std::endl;
    return 0;
}