#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <filesystem>
#include "reader.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include <TStopwatch.h>

using namespace std;
namespace fs = std::filesystem;

int main() {
    TStopwatch timer;
    timer.Start();

    // Create output directory
    fs::create_directory("output");

    // Read directories from directories.txt
    ifstream inputFile("directories.txt");
    if (!inputFile.is_open()) {
        cerr << "Error: Could not open directories.txt" << endl;
        return 1;
    }
    vector<string> directories;
    string dir;
    while (getline(inputFile, dir)) {
        if (!dir.empty()) directories.push_back(dir);
    }
    inputFile.close();
    if (directories.empty()) {
        cerr << "Error: No directories found" << endl;
        return 1;
    }

    const int maxEvents = 100000000;
    //const int maxEvents = 1000;
    int total_events_processed = 0;
    int event_count = 0;
    int trigger_electrons_events = 0;

    // Beta vs. p histograms
    TH2F* hBetaVsPAllPositive = new TH2F("hBetaVsPAllPositive", "Beta vs. Momentum (All Positive Particles, No Positrons);Momentum p (GeV);Beta", 50, 0, 10, 50, 0, 1.2);
    TH2F* hBetaVsPPID211 = new TH2F("hBetaVsPPID211", "Beta vs. Momentum (PID 211);Momentum p (GeV);Beta", 50, 0, 10, 50, 0, 1.2);
    TH2F* hBetaVsPPID321 = new TH2F("hBetaVsPPID321", "Beta vs. Momentum (PID 321);Momentum p (GeV);Beta", 50, 0, 10, 50, 0, 1.2);
    TH2F* hBetaVsPPID2212 = new TH2F("hBetaVsPPID2212", "Beta vs. Momentum (PID 2212);Momentum p (GeV);Beta", 50, 0, 10, 50, 0, 1.2);

    // Momentum bins from 1 to 8 GeV with 0.35 GeV width for chi2pid analysis
    const double p_min = 1.0;
    const double p_max = 8.0;
    const double p_bin_width = 0.35;
    const int n_p_bins = static_cast<int>((p_max - p_min) / p_bin_width);
    vector<TH1F*> hChi2PIDPionvsPBins(n_p_bins);  // For pions (PID 211)
    vector<TH1F*> hChi2PIDKaonvsPBins(n_p_bins);  // For kaons (PID 321)
    vector<TH1F*> hChi2PIDProtonvsPBins(n_p_bins); // For protons (PID 2212)

    // Initialize chi2pid histograms for each momentum bin
    for (int i = 0; i < n_p_bins; ++i) {
        double p_low = p_min + i * p_bin_width;
        double p_high = p_low + p_bin_width;
        string hist_name_pion = "hChi2PIDPion_p_" + to_string(p_low) + "_to_" + to_string(p_high);
        string hist_name_kaon = "hChi2PIDKaon_p_" + to_string(p_low) + "_to_" + to_string(p_high);
        string hist_name_proton = "hChi2PIDProton_p_" + to_string(p_low) + "_to_" + to_string(p_high);
        string hist_title_pion = "Chi2PID for Positive Pions (p = " + to_string(p_low) + " to " + to_string(p_high) + " GeV);Chi2PID;Counts";
        string hist_title_kaon = "Chi2PID for Positive Kaons (p = " + to_string(p_low) + " to " + to_string(p_high) + " GeV);Chi2PID;Counts";
        string hist_title_proton = "Chi2PID for Positive Protons (p = " + to_string(p_low) + " to " + to_string(p_high) + " GeV);Chi2PID;Counts";
        hChi2PIDPionvsPBins[i] = new TH1F(hist_name_pion.c_str(), hist_title_pion.c_str(), 100, -5, 5);
        hChi2PIDKaonvsPBins[i] = new TH1F(hist_name_kaon.c_str(), hist_title_kaon.c_str(), 100, -5, 5);
        hChi2PIDProtonvsPBins[i] = new TH1F(hist_name_proton.c_str(), hist_title_proton.c_str(), 100, -5, 5);
    }

    // Setup TTree to store REC::Particle variables
    TFile* outputFile = new TFile("output/particle_data.root", "RECREATE");
    TTree* tree = new TTree("ParticleTree", "Tree containing REC::Particle data");

    // Variables for the TTree (event-by-event, vector for all particles in the event)
    vector<int> pid;
    vector<float> px, py, pz;
    vector<float> vx, vy, vz, vt;
    vector<int> charge;
    vector<float> beta;
    vector<float> chi2pid;
    vector<int> status;

    // Branch definitions
    tree->Branch("pid", &pid);
    tree->Branch("px", &px);
    tree->Branch("py", &py);
    tree->Branch("pz", &pz);
    tree->Branch("vx", &vx);
    tree->Branch("vy", &vy);
    tree->Branch("vz", &vz);
    tree->Branch("vt", &vt);
    tree->Branch("charge", &charge);
    tree->Branch("beta", &beta);
    tree->Branch("chi2pid", &chi2pid);
    tree->Branch("status", &status);

    for (const auto& directory : directories) {
        if (event_count >= maxEvents) break;
        cout << "Processing directory: " << directory << endl;

        for (const auto& entry : fs::directory_iterator(directory)) {
            if (event_count >= maxEvents) break;

            string file_path = entry.path().string();
            if (file_path.substr(file_path.find_last_of(".") + 1) != "hipo") continue;

            cout << "Processing HIPO file: " << file_path << endl;

            hipo::reader reader;
            reader.open(file_path.c_str());
            if (!reader.is_open()) {
                cerr << "Error: Failed to open HIPO file: " << file_path << endl;
                continue;
            }

            hipo::dictionary dict;
            reader.readDictionary(dict);
            if (!dict.hasSchema("REC::Particle")) {
                cerr << "Error: REC::Particle schema not found in dictionary for file: " << file_path << endl;
                continue;
            }

            hipo::bank PART(dict.getSchema("REC::Particle"));
            hipo::event event;

            while (reader.next() && event_count < maxEvents) {
                reader.read(event);
                event.getStructure(PART);
                if (PART.getRows() == 0) continue;

                total_events_processed++;
                event_count++;

                bool hasTriggerElectron = false;

                // Check for trigger electron
                for (int i = 0; i < PART.getRows(); ++i) {
                    int pid_i = PART.getInt("pid", i);
                    int status = PART.getShort("status", i);
                    if (status == 0) continue;

                    if (pid_i == 11) { // Electron
                        bool is_trigger = (status < 0);
                        if (is_trigger) {
                            float chi2pid_tele = PART.getFloat("chi2pid", i);
                            if (chi2pid_tele > -5 && chi2pid_tele < 5) {
                                float vz_tele = PART.getFloat("vz", i);
                                if (vz_tele >= -20 && vz_tele <= 5) {
                                    if (abs(status) / 1000 == 2) {
                                        hasTriggerElectron = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }

                if (hasTriggerElectron) {
                    trigger_electrons_events++;

                    // Clear vectors for the new event
                    pid.clear();
                    px.clear();
                    py.clear();
                    pz.clear();
                    vx.clear();
                    vy.clear();
                    vz.clear();
                    vt.clear();
                    charge.clear();
                    beta.clear();
                    chi2pid.clear();
                    status.clear();

                    // Process all particles in the event
                    for (int i = 0; i < PART.getRows(); ++i) {
                        int pid_i = PART.getInt("pid", i);
                        int charge_i = PART.getByte("charge", i);
                        float chi2pid_i = PART.getFloat("chi2pid", i);
                        float px_i = PART.getFloat("px", i);
                        float py_i = PART.getFloat("py", i);
                        float pz_i = PART.getFloat("pz", i);
                        float p = std::sqrt(px_i * px_i + py_i * py_i + pz_i * pz_i);
                        float beta_i = PART.getFloat("beta", i);
                        int status_i = PART.getShort("status", i);
                        if (status_i == 0) continue;
                        if (charge_i <= 0) continue; 
                        // Fill TTree variables
                        pid.push_back(pid_i);
                        px.push_back(px_i);
                        py.push_back(py_i);
                        pz.push_back(pz_i);
                        vx.push_back(PART.getFloat("vx", i));
                        vy.push_back(PART.getFloat("vy", i));
                        vz.push_back(PART.getFloat("vz", i));
                        vt.push_back(PART.getFloat("vt", i));
                        charge.push_back(charge_i);
                        beta.push_back(beta_i);
                        chi2pid.push_back(chi2pid_i);
                        status.push_back(status_i);

                        // Fill histograms for positive charge particles
                        if (charge_i > 0 && beta_i > 0 && beta_i <= 1.2 && p > 0) {
                            if (pid_i != -11) { // Exclude positrons
                                hBetaVsPAllPositive->Fill(p, beta_i);
                            }
                            if (pid_i == 211) {
                                hBetaVsPPID211->Fill(p, beta_i);
                                if (chi2pid_i != 9999.0 && p >= p_min && p < p_max) {
                                    int bin_index = static_cast<int>((p - p_min) / p_bin_width);
                                    if (bin_index >= 0 && bin_index < n_p_bins) {
                                        hChi2PIDPionvsPBins[bin_index]->Fill(chi2pid_i);
                                    }
                                }
                            }
                            if (pid_i == 321) {
                                hBetaVsPPID321->Fill(p, beta_i);
                                if (chi2pid_i != 9999.0 && p >= p_min && p < p_max) {
                                    int bin_index = static_cast<int>((p - p_min) / p_bin_width);
                                    if (bin_index >= 0 && bin_index < n_p_bins) {
                                        hChi2PIDKaonvsPBins[bin_index]->Fill(chi2pid_i);
                                    }
                                }
                            }
                            if (pid_i == 2212) {
                                hBetaVsPPID2212->Fill(p, beta_i);
                                if (chi2pid_i != 9999.0 && p >= p_min && p < p_max) {
                                    int bin_index = static_cast<int>((p - p_min) / p_bin_width);
                                    if (bin_index >= 0 && bin_index < n_p_bins) {
                                        hChi2PIDProtonvsPBins[bin_index]->Fill(chi2pid_i);
                                    }
                                }
                            }
                        }
                    }

                    // Fill the TTree for this event
                    tree->Fill();
                }

                if (event_count % 1000 == 0) {
                    cout << "Processed " << event_count << " events" << endl;
                }
            }
        }
    }

    // Write histograms and TTree to the ROOT file
    hBetaVsPAllPositive->Write();
    hBetaVsPPID211->Write();
    hBetaVsPPID321->Write();
    hBetaVsPPID2212->Write();
    for (int i = 0; i < n_p_bins; ++i) {
        hChi2PIDPionvsPBins[i]->Write();
        hChi2PIDKaonvsPBins[i]->Write();
        hChi2PIDProtonvsPBins[i]->Write();
    }
    tree->Write();
    outputFile->Close();
    delete outputFile;

    timer.Stop();
    double execution_time = timer.RealTime();

    // Write basic stats to CSV
    ofstream csvFile("output/stats_basic.csv");
    if (!csvFile.is_open()) {
        cerr << "Error: Could not open output/stats_basic.csv for writing!" << endl;
        return 1;
    }

    csvFile << "Statistic,Value\n";
    csvFile << "Total events processed," << total_events_processed << "\n";
    csvFile << "Events with trigger electrons," << trigger_electrons_events << "\n";
    csvFile << "Execution time (seconds)," << execution_time << "\n";
    csvFile.close();
    cout << "Basic statistics saved to output/stats_basic.csv" << endl;

    // Clean up
    delete hBetaVsPAllPositive;
    delete hBetaVsPPID211;
    delete hBetaVsPPID321;
    delete hBetaVsPPID2212;
    for (int i = 0; i < n_p_bins; ++i) {
        delete hChi2PIDPionvsPBins[i];
        delete hChi2PIDKaonvsPBins[i];
        delete hChi2PIDProtonvsPBins[i];
    }

    return 0;
}