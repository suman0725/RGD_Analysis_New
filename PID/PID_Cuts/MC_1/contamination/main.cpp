#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <cmath>

using namespace std;
using namespace ROOT;

// Timing cut functions
double getdtCutNeg(double p) {
    return -0.232 + (-6.405) * exp(-8.123 * p);
}
double getdtCutPos(double p) {
    return 0.225 + 5.245 * exp(-7.298 * p);
}

int main() {
    auto start = chrono::high_resolution_clock::now();

    // Input file
    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/new_RGD_Analysis/PID/PID_Cuts/MC_1/mc_pion_candidates.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file" << endl;
        return 1;
    }

    ROOT::RDataFrame df_all("EB_pos_pion_assumed_mc", file);

    // Count total unmatched tracks for diagnostics
    auto n_unmatched_total = df_all.Filter("is_matched == 0").Count();
    cout << "Total unmatched tracks: " << *n_unmatched_total << endl;

    // Define momentum bins (0.4 to 8.2 GeV/c, 0.3 GeV/c width)
    const double pMin = 0.4, pMax = 8.2, binWidth = 0.3;
    const int nBins = static_cast<int>((pMax - pMin) / binWidth);
    vector<double> pLowEdges(nBins), pHighEdges(nBins);
    for (int i = 0; i < nBins; ++i) {
        pLowEdges[i] = pMin + i * binWidth;
        pHighEdges[i] = pLowEdges[i] + binWidth;
    }

    // Output file
    ofstream outFile("contamination_results_all_bins.txt");
    if (!outFile.is_open()) {
        cerr << "Error: Cannot open contamination_results_all_bins.txt" << endl;
        file->Close();
        delete file;
        return 1;
    }
    outFile << fixed << setprecision(4);
    outFile << "Contamination Results for MC Positive Tracks\n";
    outFile << "----------------------------------------\n";
    outFile << "Bin | p_low | p_high | p_mid | N_pion_before | N_pion_after | N_kaon | N_proton | N_other | Total | Eff (%) | Kaon Cont (%) | Proton Cont (%) | Other Cont (%) | N_unmatched\n";
    outFile << "----------------------------------------\n";

    // Process each bin
    for (int i = 0; i < nBins; ++i) {
        double pLow = pLowEdges[i], pHigh = pHighEdges[i];

        // Filter tracks in momentum bin and ensure truth-matching
        auto df_bin = df_all.Filter([pLow, pHigh](float p, int is_matched) {
            return p >= pLow && p < pHigh && is_matched == 1;
        }, {"p", "is_matched"});

        // Count unmatched tracks in this bin for diagnostics
        auto df_unmatched = df_all.Filter([pLow, pHigh](float p, int is_matched) {
            return p >= pLow && p < pHigh && is_matched == 0;
        }, {"p", "is_matched"});
        auto n_unmatched = df_unmatched.Count();

        // Calculate mean momentum
        auto p_mean_action = df_bin.Mean("p");
        double p_mid = *p_mean_action;

        // Count true pions before dt cut
        auto n_pion_before = df_bin.Filter("true_pid == 211").Count();

        // Apply dt cut
        auto df_after_dt = df_bin.Filter([pLow](float p, float dt) {
            return dt > getdtCutNeg(pLow) && dt < getdtCutPos(pLow);
        }, {"p", "dt"});

        // Count tracks after dt cut
        auto n_total = df_after_dt.Count();
        auto n_pion_after = df_after_dt.Filter("true_pid == 211").Count();
        auto n_kaon = df_after_dt.Filter("true_pid == 321").Count();
        auto n_proton = df_after_dt.Filter("true_pid == 2212").Count();
        auto n_other = df_after_dt.Filter("true_pid != 211 && true_pid != 321 && true_pid != 2212").Count();

        // Wait for counts to resolve
        double N_pion_before = *n_pion_before;
        double N_pion_after = *n_pion_after;
        double N_kaon = *n_kaon;
        double N_proton = *n_proton;
        double N_other = *n_other;
        double N_total = *n_total;
        double N_unmatched = *n_unmatched;

        // Calculate efficiency and contamination
        double efficiency = N_pion_before > 0 ? N_pion_after / N_pion_before : 0.0;
        double kaon_contamination = N_total > 0 ? N_kaon / N_total : 0.0;
        double proton_contamination = N_total > 0 ? N_proton / N_total : 0.0;
        double other_contamination = N_total > 0 ? N_other / N_total : 0.0;

        // Poisson errors for counts
        double sigma_pion_before = sqrt(N_pion_before);
        double sigma_pion_after = sqrt(N_pion_after);
        double sigma_kaon = sqrt(N_kaon);
        double sigma_proton = sqrt(N_proton);
        double sigma_other = sqrt(N_other);
        double sigma_total = sqrt(N_total);

        // Error propagation for ratios (assuming zero covariance)
        double eff_error = (N_pion_before > 0 && N_pion_after > 0) ?
            efficiency * sqrt(1.0 / N_pion_after + 1.0 / N_pion_before) : 0.0;
        double kaon_cont_error = (N_total > 0 && N_kaon > 0) ?
            kaon_contamination * sqrt(1.0 / N_kaon + 1.0 / N_total) : 0.0;
        double proton_cont_error = (N_total > 0 && N_proton > 0) ?
            proton_contamination * sqrt(1.0 / N_proton + 1.0 / N_total) : 0.0;
        double other_cont_error = (N_total > 0 && N_other > 0) ?
            other_contamination * sqrt(1.0 / N_other + 1.0 / N_total) : 0.0;

        // Output results with ± format
        outFile << i + 1 << " | "
                << pLow << " | " << pHigh << " | " << p_mid << " | "
                << N_pion_before << " ± " << sigma_pion_before << " | "
                << N_pion_after << " ± " << sigma_pion_after << " | "
                << N_kaon << " ± " << sigma_kaon << " | "
                << N_proton << " ± " << sigma_proton << " | "
                << N_other << " ± " << sigma_other << " | "
                << N_total << " ± " << sigma_total << " | "
                << efficiency * 100 << " ± " << eff_error * 100 << " | "
                << kaon_contamination * 100 << " ± " << kaon_cont_error * 100 << " | "
                << proton_contamination * 100 << " ± " << proton_cont_error * 100 << " | "
                << other_contamination * 100 << " ± " << other_cont_error * 100 << " | "
                << N_unmatched << "\n";
    }

    outFile << "----------------------------------------\n";
    outFile.close();

    file->Close();
    delete file;

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Total program time: " << duration.count() / 1000.0 << " seconds" << endl;

    return 0;
}