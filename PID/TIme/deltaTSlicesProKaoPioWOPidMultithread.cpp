#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <iostream>

struct ParticleData {
    TTree* tree;
    const char* name;
    int color;
};

int main() {
    TFile* file = TFile::Open("particles_data.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file." << std::endl;
        return 1;
    }

    ParticleData particles[] = {
        {dynamic_cast<TTree*>(file->Get("pion")), "Pion", kRed},
        {dynamic_cast<TTree*>(file->Get("kaon")), "Kaon", kBlue},
        {dynamic_cast<TTree*>(file->Get("proton")), "Proton", kGreen}
    };

    // Validate trees
    for (auto& p : particles) {
        if (!p.tree) {
            std::cerr << "Error: Could not find tree for " << p.name << std::endl;
            file->Close();
            return 1;
        }
        p.tree->SetBranchStatus("*", 0); // Disable all branches
        p.tree->SetBranchStatus("delta_t", 1); // Enable only needed branches
        p.tree->SetBranchStatus("momentum", 1);
    }

    const Float_t momentum_min = 0.2;
    const Float_t momentum_max = 7;
    const Float_t slice_size = 0.05;
    const int num_slices = static_cast<int>((momentum_max - momentum_min) / slice_size);

    TString pdf_file = "delta_t_slices_all_particles_MT.pdf";
    TCanvas* canvas = new TCanvas("canvas", "Momentum Slices", 800, 600);

    for (int i = 0; i < num_slices; i++) {
        const Float_t slice_start = momentum_min + i * slice_size;
        const Float_t slice_end = slice_start + slice_size;

        TString hist_title = Form("Momentum range [%.2f, %.2f] GeV", slice_start, slice_end);
        TH1F* hist_pion = new TH1F("hist_pion", hist_title, 100, -3, 3);
        TH1F* hist_kaon = new TH1F("hist_kaon", hist_title, 100, -3, 3);
        TH1F* hist_proton = new TH1F("hist_proton", hist_title, 100, -3, 3);

        // Process PION tree
        {
            TTreeReader reader(particles[0].tree);
            TTreeReaderValue<Float_t> delta_t(reader, "delta_t");
            TTreeReaderValue<Float_t> momentum(reader, "momentum");
            while (reader.Next()) {
                if (*momentum >= slice_start && *momentum < slice_end) {
                    hist_pion->Fill(*delta_t);
                }
            }
        }

        // Process KAON tree
        {
            TTreeReader reader(particles[1].tree);
            TTreeReaderValue<Float_t> delta_t(reader, "delta_t");
            TTreeReaderValue<Float_t> momentum(reader, "momentum");
            while (reader.Next()) {
                if (*momentum >= slice_start && *momentum < slice_end) {
                    hist_kaon->Fill(*delta_t);
                }
            }
        }

        // Process PROTON tree
        {
            TTreeReader reader(particles[2].tree);
            TTreeReaderValue<Float_t> delta_t(reader, "delta_t");
            TTreeReaderValue<Float_t> momentum(reader, "momentum");
            while (reader.Next()) {
                if (*momentum >= slice_start && *momentum < slice_end) {
                    hist_proton->Fill(*delta_t);
                }
            }
        }

        // Draw histograms
        hist_pion->SetLineColor(kRed);
        hist_kaon->SetLineColor(kBlue);
        hist_proton->SetLineColor(kGreen);

        hist_pion->Draw("HIST");
        hist_kaon->Draw("HIST SAME");
        hist_proton->Draw("HIST SAME");

        TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(hist_pion, "Pion", "l");
        legend->AddEntry(hist_kaon, "Kaon", "l");
        legend->AddEntry(hist_proton, "Proton", "l");
        legend->Draw();

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