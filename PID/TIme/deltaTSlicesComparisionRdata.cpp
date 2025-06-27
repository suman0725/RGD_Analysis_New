#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <ROOT/TThreadExecutor.hxx>
#include <iostream>
#include <future>  // Required for std::future
#include <vector>  // Required for std::vector
#include <TLegend.h>

const Float_t momentum_min = 0.2;
const Float_t momentum_max = 7.0;
const Float_t slice_size = 0.05;
const int num_slices = static_cast<int>((momentum_max - momentum_min) / slice_size);

struct HistogramPair {
    TH1F* hist_with_pid;
    TH1F* hist_without_pid;
};

// Function to process momentum slices
HistogramPair processSlice(TTree* tree1, TTree* tree2, int slice_index) {
    Float_t slice_start = momentum_min + slice_index * slice_size;
    Float_t slice_end = slice_start + slice_size;

    // Histograms for this slice
    TString hist_name1 = Form("DeltaT_slice_with_pid_%d", slice_index);
    TH1F* hist_with_pid = new TH1F(hist_name1, "", 100, -3, 3);

    TString hist_name2 = Form("DeltaT_slice_without_pid_%d", slice_index);
    TH1F* hist_without_pid = new TH1F(hist_name2, "", 100, -3, 3);

    // TTreeReader for better performance
    TTreeReader reader1(tree1);
    TTreeReaderValue<Float_t> delta_t1(reader1, "delta_t");
    TTreeReaderValue<Float_t> momentum1(reader1, "momentum");

    TTreeReader reader2(tree2);
    TTreeReaderValue<Float_t> delta_t2(reader2, "delta_t");
    TTreeReaderValue<Float_t> momentum2(reader2, "momentum");

    // Fill histograms for tree1 (with PID)
    while (reader1.Next()) {
        if (*momentum1 >= slice_start && *momentum1 < slice_end) {
            hist_with_pid->Fill(*delta_t1);
        }
    }

    // Fill histograms for tree2 (without PID)
    while (reader2.Next()) {
        if (*momentum2 >= slice_start && *momentum2 < slice_end) {
            hist_without_pid->Fill(*delta_t2);
        }
    }

    return {hist_with_pid, hist_without_pid};
}

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
        return 1;
    }

    // Use multithreading for parallel processing
    ROOT::TThreadExecutor executor;
    std::vector<std::future<HistogramPair>> results;

    for (int i = 0; i < num_slices; i++) {
         results.push_back(std::async(std::launch::async, processSlice, tree1, tree2, i));
    }

    // Create PDF for histograms
    TString pdf_file = "delta_t_comparison_rdata.pdf";
    TCanvas* canvas = new TCanvas("canvas", "Momentum Slices", 800, 600);

    // Retrieve results and plot histograms
    for (int i = 0; i < num_slices; i++) {
        auto pair = results[i].get();  // Retrieve histograms

        pair.hist_with_pid->SetLineColor(kRed);
        pair.hist_with_pid->SetFillColorAlpha(kRed, 0.5);
        pair.hist_with_pid->SetFillStyle(3001);

        pair.hist_without_pid->SetLineColor(kBlack);
        pair.hist_without_pid->SetFillColorAlpha(kBlack, 0.5);
        pair.hist_without_pid->SetFillStyle(3004);

        pair.hist_without_pid->Draw("HIST");
        pair.hist_with_pid->Draw("HIST SAME");

        TLegend* legend = new TLegend(0.65, 0.75, 0.9, 0.9); // Adjust position if needed
        legend->SetBorderSize(0);  // No border
        legend->SetFillStyle(0);   // Transparent background
        legend->AddEntry(pair.hist_without_pid, Form("Without PID (Entries: %d)", (int)pair.hist_without_pid->GetEntries()), "f");
        legend->AddEntry(pair.hist_with_pid, Form("With PID (Entries: %d)", (int)pair.hist_with_pid->GetEntries()), "f");
        legend->Draw();

        if (i == 0)
            canvas->Print(pdf_file + "(", "pdf");
        else if (i == num_slices - 1)
            canvas->Print(pdf_file + ")", "pdf");
        else
            canvas->Print(pdf_file, "pdf");

        delete pair.hist_with_pid;
        delete pair.hist_without_pid;
    }

    canvas->Close();
    delete canvas;
    file1->Close();
    file2->Close();
    delete file1;
    delete file2;

    std::cout << "Histograms saved to " << pdf_file << std::endl;
    return 0;
}