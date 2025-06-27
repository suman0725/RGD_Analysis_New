#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <filesystem>
#include <map>
#include "reader.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"

using namespace std;
namespace fs = std::filesystem;

int main() {
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

    const int maxEvents = 100000;
    int event_count = 0;
    int total_events_processed = 0;
    int pre_trigger_electron_count = 0;
    int trigger_elec_events = 0;

    // Pre-trigger counters
    int total_particles = 0;
    int total_positive_charge = 0;
    int pre_trigger_total_negative_charge = 0;
    int total_neutral_particles = 0;

    // PID and chi2pid counters
    map<int, int> pid_counts;
    map<int, int> chi2pid_valid_counts;
    map<int, int> chi2pid_invalid_counts;
    map<int, int> positive_pid_counts;
    map<int, int> postive_chi2pid_valid_counts;
    map<int, int> pre_trigger_positive_chi2pid_invalid_counts;

    // List of events with a trigger electron
    vector<int> trigger_electron_event_indices;

    // Electron cut statistics
    int trigger_elec_electrons_pid_11 = 0;
    int trigger_elec_status_negative = 0;
    int trigger_elec_vz_range = 0;
    int trigger_elec_final = 0;

    // Histogram for electron momentum (before final trigger selection)
    TH1F* hElectronMomentum = new TH1F("hElectronMomentum", "Electron Momentum (Before Final Trigger Selection);Momentum (GeV);Counts", 200, 0, 10);
    hElectronMomentum->SetLineColor(kBlue);
    hElectronMomentum->SetFillColor(kBlue);
    hElectronMomentum->SetFillStyle(3004);

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

                // Process particles for pre-trigger statistics
                for (int i = 0; i < PART.getRows(); ++i) {
                    total_particles++;

                    int pid = PART.getInt("pid", i);
                    pid_counts[pid]++;

                    float chi2pid = PART.getFloat("chi2pid", i);
                    if (chi2pid != 9999.0) {
                        chi2pid_valid_counts[pid]++;
                    } else {
                        chi2pid_invalid_counts[pid]++;
                    }

                    int charge = PART.getByte("charge", i);
                    if (charge == 0) {
                        total_neutral_particles++;
                    } else if (charge > 0) {
                        total_positive_charge++;
                        positive_pid_counts[pid]++;
                        if (chi2pid != 9999.0) {
                            postive_chi2pid_valid_counts[pid]++;
                        } else {
                            pre_trigger_positive_chi2pid_invalid_counts[pid]++;
                        }
                    } else if (charge < 0) {
                        pre_trigger_total_negative_charge++;
                    }

                    int status = PART.getShort("status", i);
                    if (status == 0) continue;

                    if (pid == 11) { // Electron
                        trigger_elec_electrons_pid_11++;
                        pre_trigger_electron_count++;

                        float px = PART.getFloat("px", i);
                        float py = PART.getFloat("py", i);
                        float pz = PART.getFloat("pz", i);
                        float p_e = std::sqrt(px * px + py * py + pz * pz);
                        hElectronMomentum->Fill(p_e);

                        bool is_trigger = (status < 0);
                        if (is_trigger) {
                            trigger_elec_status_negative++;

                            float vz_tele = PART.getFloat("vz", i);
                            if (vz_tele >= -20 && vz_tele <= 5) {
                                trigger_elec_vz_range++;

                                if (abs(status) / 1000 == 2) {
                                    trigger_elec_final++;
                                    hasTriggerElectron = true;
                                }
                            }
                        }
                    }
                }

                if (hasTriggerElectron) {
                    trigger_elec_events++;
                    trigger_electron_event_indices.push_back(event_count - 1);
                }

                if (event_count % 1000 == 0) {
                    cout << "Processed " << event_count << " events" << endl;
                }
            }
        }
    }

    // Validation checks
    int sum_pid_counts = 0;
    for (const auto& pid_count : pid_counts) {
        sum_pid_counts += pid_count.second;
    }
    bool pid_count_matches = (sum_pid_counts == total_particles);
    string pid_validation_message = pid_count_matches ? "Yes" : "No (Possible data inconsistency)";

    int sum_valid_chi2pid = 0, sum_invalid_chi2pid = 0;
    for (const auto& pid_count : pid_counts) {
        int pid = pid_count.first;
        sum_valid_chi2pid += chi2pid_valid_counts[pid];
        sum_invalid_chi2pid += chi2pid_invalid_counts[pid];
    }
    bool chi2pid_sum_matches = (sum_valid_chi2pid + sum_invalid_chi2pid == sum_pid_counts);
    string chi2pid_validation_message = chi2pid_sum_matches ? "Yes" : "No (Possible data inconsistency)";

    // Write pre-trigger statistics to CSV
    ofstream csvFile("output/stats.csv");
    if (!csvFile.is_open()) {
        cerr << "Error: Could not open output/stats.csv for writing!" << endl;
        return 1;
    }

    csvFile << "Statistic,Value\n";
    csvFile << "Total events processed," << total_events_processed << "\n";
    csvFile << "Pre-trigger total particles," << total_particles << "\n";
    csvFile << "Pre-trigger total neutral particles (charge == 0)," << total_neutral_particles << "\n";
    csvFile << "Pre-trigger total positive charged particles (charge > 0)," << total_positive_charge << "\n";
    csvFile << "Pre-trigger total negative charged particles (charge < 0)," << pre_trigger_total_negative_charge << "\n";
    csvFile << "Pre-trigger total electrons (pid == 11)," << pre_trigger_electron_count << "\n";
    csvFile << "Trigger electron events," << trigger_elec_events << "\n";

    csvFile << "\nPre-Trigger PID Counts (REC::Particle),\n";
    csvFile << "PID,Count\n";
    for (const auto& pid_count : pid_counts) {
        csvFile << pid_count.first << "," << pid_count.second << "\n";
    }

    csvFile << "\nPre-Trigger PID Count Validation,\n";
    csvFile << "Sum of PID counts," << sum_pid_counts << "\n";
    csvFile << "Matches pre-trigger total particles?," << pid_validation_message << "\n";

    csvFile << "\nTrigger Electron Cut Statistics,\n";
    csvFile << "Cut,Count\n";
    csvFile << "Electrons with pid == 11," << trigger_elec_electrons_pid_11 << "\n";
    csvFile << "Electrons with pid == 11 and status < 0," << trigger_elec_status_negative << "\n";
    csvFile << "Electrons with pid == 11, status < 0, and vz_tele in [-20, 5]," << trigger_elec_vz_range << "\n";
    csvFile << "Electrons with pid == 11, status < 0, vz_tele in [-20, 5], and abs(status) / 1000 == 2," << trigger_elec_final << "\n";

    csvFile << "\nPre-Trigger PID Counts with Valid Chi2PID (REC::Particle),\n";
    csvFile << "PID,Count\n";
    for (const auto& pid_count : pid_counts) {
        int pid = pid_count.first;
        int valid_count = chi2pid_valid_counts[pid];
        csvFile << pid << "," << valid_count << "\n";
    }

    csvFile << "\nPre-Trigger PID Counts with Invalid Chi2PID (REC::Particle),\n";
    csvFile << "PID,Count\n";
    for (const auto& pid_count : pid_counts) {
        int pid = pid_count.first;
        int invalid_count = chi2pid_invalid_counts[pid];
        csvFile << pid << "," << invalid_count << "\n";
    }

    csvFile << "\nPre-Trigger Chi2PID Count Validation,\n";
    csvFile << "Sum of particles with valid Chi2PID," << sum_valid_chi2pid << "\n";
    csvFile << "Sum of particles with invalid Chi2PID," << sum_invalid_chi2pid << "\n";
    csvFile << "Total (valid + invalid)," << (sum_valid_chi2pid + sum_invalid_chi2pid) << "\n";
    csvFile << "Matches sum of pre-trigger PID counts?," << chi2pid_validation_message << "\n";

    csvFile.close();
    cout << "Pre-trigger statistics saved to output/stats.csv" << endl;

    // Save event indices with trigger electrons
    ofstream eventFile("output/trigger_electron_events.txt");
    if (!eventFile.is_open()) {
        cerr << "Error: Could not open output/trigger_electron_events.txt for writing!" << endl;
        return 1;
    }
    for (int idx : trigger_electron_event_indices) {
        eventFile << idx << "\n";
    }
    eventFile.close();
    cout << "Trigger electron event indices saved to output/trigger_electron_events.txt" << endl;

    // Plot electron momentum (before final trigger selection)
    TCanvas* canvas = new TCanvas("canvas", "Electron Momentum (Before Final Trigger Selection)", 800, 600);
    hElectronMomentum->Draw("HIST");
    canvas->Print("output/electron_momentum.pdf");
    cout << "Electron momentum histogram saved to output/electron_momentum.pdf" << endl;

    // Clean up
    delete hElectronMomentum;
    delete canvas;

    return 0;
}