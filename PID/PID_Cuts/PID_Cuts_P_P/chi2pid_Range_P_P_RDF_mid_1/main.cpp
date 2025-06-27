#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <algorithm>

using namespace std;

int main() {
    auto start = chrono::high_resolution_clock::now();

    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/Skim/pkptreeCxC_9_test_modified.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file\n";
        return 1;
    }

    TTree* tree = dynamic_cast<TTree*>(file->Get("EB_all_pion_assumed"));
    if (!tree) {
        cerr << "Error: Cannot find tree\n";
        file->Close();
        return 1;
    }

    float p, beta, recomputed_chi2pid;
    tree->SetBranchAddress("p", &p);
    tree->SetBranchAddress("beta", &beta);
    tree->SetBranchAddress("recomputed_chi2pid", &recomputed_chi2pid);

    vector<double> momenta = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0};
    vector<double> uncappedThresholds, cappedThresholds, pCenters, uncappedErrors, cappedErrors;

    Long64_t nEntries = tree->GetEntries();
    const double pTolerance = 0.05;
    const double betaTolerance = 0.0001;

    for (double targetP : momenta) {
        vector<double> chi2pid_values;
        for (Long64_t j = 0; j < nEntries; j++) {
            tree->GetEntry(j);
            if (p >= targetP - pTolerance && p <= targetP + pTolerance && beta >= 0.0 && beta <= 1.1) {
                float mPi = 0.1396f;
                float mK = 0.4937f;
                float betaPi = targetP / sqrtf(targetP * targetP + mPi * mPi);
                float betaK = targetP / sqrtf(targetP * targetP + mK * mK);
                float betaMid = (betaPi + betaK) / 2.0f;
                if (fabs(beta - betaMid) < betaTolerance) {
                    double chi2pid = recomputed_chi2pid;
                    if (chi2pid > -100.0 && chi2pid < 1000.0) {
                        chi2pid_values.push_back(chi2pid);
                    }
                }
            }
        }

        double uncappedThreshold = 0.0, cappedThreshold = 0.0;
        double uncappedError = 0.0, cappedError = 0.0;
        if (!chi2pid_values.empty()) {
            uncappedThreshold = *max_element(chi2pid_values.begin(), chi2pid_values.end()); // Real maximum chi2pid
            cappedThreshold = (uncappedThreshold > 5.0) ? 5.0 : uncappedThreshold; // Cap at 5 if above
            if (chi2pid_values.size() > 1) {
                double sum = 0.0, sum_sq = 0.0;
                for (double val : chi2pid_values) {
                    sum += val;
                    sum_sq += val * val;
                }
                double mean = sum / chi2pid_values.size();
                double variance = sum_sq / chi2pid_values.size() - mean * mean;
                uncappedError = sqrt(variance); // Standard deviation for uncapped
                cappedError = uncappedError; // Keep original error for capped file
                if (uncappedError == 0.0) {
                    uncappedError = 0.1;
                    cappedError = 0.1;
                }
            } else {
                uncappedError = 0.1;
                cappedError = 0.1;
            }
        }
        if (uncappedThreshold > 0.0) {
            uncappedThresholds.push_back(uncappedThreshold);
            cappedThresholds.push_back(cappedThreshold);
            pCenters.push_back(targetP);
            uncappedErrors.push_back(uncappedError);
            cappedErrors.push_back(cappedError);
            cout << "p = " << targetP << ": Uncapped Threshold = " << uncappedThreshold 
                 << ", Capped Threshold = " << cappedThreshold 
                 << ", Uncapped Error = " << uncappedError 
                 << ", Capped Error = " << cappedError << endl;
        }
    }

    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/PID_Cuts_P_P", kTRUE);

    // Save uncapped thresholds and errors
    ofstream uncappedFile("output/PID_Cuts_P_P/chi2pid_thresholds_uncapped.txt");
    if (uncappedFile.is_open()) {
        uncappedFile << "Momentum (GeV/c), Threshold, Error\n";
        for (size_t i = 0; i < pCenters.size(); i++) {
            uncappedFile << pCenters[i] << ", " << uncappedThresholds[i] << ", " << uncappedErrors[i] << "\n";
        }
        uncappedFile.close();
        cout << "Uncapped thresholds saved to chi2pid_thresholds_uncapped.txt" << endl;
    } else {
        cerr << "Error: Unable to open chi2pid_thresholds_uncapped.txt for writing" << endl;
    }

    // Save capped thresholds with original errors
    ofstream cappedFile("output/PID_Cuts_P_P/chi2pid_thresholds_capped.txt");
    if (cappedFile.is_open()) {
        cappedFile << "Momentum (GeV/c), Threshold, Error\n";
        for (size_t i = 0; i < pCenters.size(); i++) {
            cappedFile << pCenters[i] << ", " << cappedThresholds[i] << ", " << cappedErrors[i] << "\n";
        }
        cappedFile.close();
        cout << "Capped thresholds saved to chi2pid_thresholds_capped.txt" << endl;
    } else {
        cerr << "Error: Unable to open chi2pid_thresholds_capped.txt for writing" << endl;
    }

    file->Close();
    delete file;

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Total program time: " << duration.count() / 1000.0 << " seconds" << endl;

    cout << "Program completed successfully." << endl;
    return 0;
}