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
#include "TF1.h"
#include "TMath.h"
#include <TStopwatch.h>
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"

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

    const int maxEvents = 100000;
    int event_count = 0;
    int total_events_processed = 0;
    int pre_trigger_electron_count = 0;
    int trigger_elec_events = 0;

    // Counters for particle statistics (before trigger electron selection)
    int pre_trigger_total_particles = 0;
    int pre_trigger_total_positive_charge = 0;
    int pre_trigger_total_negative_charge = 0;
    int pre_trigger_total_neutral_particles = 0;
    int post_trigger_total_positive_charge = 0;
    int post_trigger_total_negative_charge = 0;

    // Counters for beta and p validation (after trigger electron selection)
    int post_trigger_positive_beta_invalid = 0;
    int post_trigger_positive_beta_1_to_1p2 = 0;
    int post_trigger_positive_beta_p_valid = 0;

    // Counter for particles within the pion mass range (after trigger electron selection)
    int post_trigger_pions_in_mass_range = 0;
    int post_trigger_pions_in_mass_range_updated = 0;
    int post_trigger_pions_pid_211_in_mass_range = 0;
    int post_trigger_pions_pid_211_in_mass_range_updated = 0;
    int post_trigger_pions_pid_211_outside_mass_range = 0;

    // Region-separated counters for particles in the pion mass range (after trigger electron selection)
    int post_trigger_pions_in_mass_range_updated_forward = 0;
    int post_trigger_pions_pid_211_in_mass_range_updated_forward = 0;
    int post_trigger_pions_in_mass_range_updated_central = 0;
    int post_trigger_pions_pid_211_in_mass_range_updated_central = 0;

    // Diagnostic counters for chi2pid asymmetry (after trigger electron selection)
    int post_trigger_pid_211_chi2pid_positive = 0; // pid == 211 with chi2pid > 0
    int post_trigger_pid_211_chi2pid_negative = 0; // pid == 211 with chi2pid < 0
    int post_trigger_pid_211_chi2pid_zero = 0;     // pid == 211 with chi2pid == 0

    // Vectors to store masses, momenta, betas, PIDs, chi2pids, and regions of positive hadrons (after trigger electron selection)
    vector<float> post_trigger_positive_hadron_masses;
    vector<int> post_trigger_positive_hadron_pids;
    vector<float> post_trigger_positive_hadron_chi2pids;
    vector<bool> post_trigger_positive_hadron_is_forward;
    vector<bool> post_trigger_positive_hadron_is_central;

    // Map to count particles by PID (before trigger electron selection)
    map<int, int> pre_trigger_pid_counts;
    map<int, int> pre_trigger_chi2pid_valid_counts;
    map<int, int> pre_trigger_chi2pid_invalid_counts;

    // Counters for positive charged particles by PID and chi2pid validity (before trigger electron selection)
    map<int, int> pre_trigger_positive_pid_counts;           // Total positive charged particles by PID
    map<int, int> pre_trigger_positive_chi2pid_valid_counts; // Positive charged particles with valid chi2pid by PID
    map<int, int> pre_trigger_positive_chi2pid_invalid_counts; // Positive charged particles with invalid chi2pid by PID

    // Counters for positive charged particles by PID, chi2pid validity, and region (after trigger electron selection)
    map<int, int> post_trigger_positive_pid_counts_forward;           // Forward region
    map<int, int> post_trigger_positive_chi2pid_valid_counts_forward;
    map<int, int> post_trigger_positive_chi2pid_invalid_counts_forward;
    map<int, int> post_trigger_positive_pid_counts_central;           // Central region
    map<int, int> post_trigger_positive_chi2pid_valid_counts_central;
    map<int, int> post_trigger_positive_chi2pid_invalid_counts_central;

    // List of PIDs to visualize
    vector<int> pids_to_visualize = {0, 2212, 321, 211, -11, 45};

    // Counters for electron cut statistics (during trigger electron selection)
    int trigger_elec_electrons_pid_11 = 0;
    int trigger_elec_status_negative = 0;
    int trigger_elec_vz_range = 0;
    int trigger_elec_final = 0;

    // Counters for positive pions in events with a trigger electron (after trigger electron selection)
    int post_trigger_positive_pions_all_regions = 0;
    int post_trigger_positive_pions_forward = 0;
    int post_trigger_positive_pions_central = 0;

    // Histogram for electron momentum
    TH1F* hElectronMomentum = new TH1F("hElectronMomentum", "Electron Momentum (REC::PID);Momentum (GeV);Counts", 200, 0, 10);
    hElectronMomentum->SetLineColor(kBlue);
    hElectronMomentum->SetFillColor(kBlue);
    hElectronMomentum->SetFillStyle(3004);

    // Histogram for mass of all positive hadrons (after trigger electron selection)
    TH1F* hMassPositiveAll = new TH1F("hMassPositiveAll", "Mass of All Positive Hadrons (Post-Trigger);Mass (GeV/c^2);Counts", 100, 0, 2.5);

    // Histograms for chi2pid of positive pions (PID 211) in Forward and Central regions (after trigger electron selection)
    TH1F* hChi2PIDPositivePionsForward = new TH1F("hChi2PIDPositivePionsForward", "Chi2PID for Positive Pions (PID 211, Forward, Post-Trigger);Chi2PID;Counts", 100, -5, 5);
    TH1F* hChi2PIDPositivePionsCentral = new TH1F("hChi2PIDPositivePionsCentral", "Chi2PID for Positive Pions (PID 211, Central, Post-Trigger);Chi2PID;Counts", 100, -5, 5);

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
                std::vector<std::string> schemas = dict.getSchemaList();
                cout << "Available schemas in dictionary:" << endl;
                for (const auto& schema : schemas) {
                    cout << "  " << schema << endl;
                }
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

                // First loop: Identify the trigger electron
                for (int i = 0; i < PART.getRows(); ++i) {
                    pre_trigger_total_particles++;

                    int pid = PART.getInt("pid", i);
                    pre_trigger_pid_counts[pid]++;

                    float chi2pid = PART.getFloat("chi2pid", i);
                    if (chi2pid != 9999.0) {
                        pre_trigger_chi2pid_valid_counts[pid]++;
                    } else {
                        pre_trigger_chi2pid_invalid_counts[pid]++;
                    }

                    int charge = PART.getByte("charge", i);
                    if (charge == 0) {
                        pre_trigger_total_neutral_particles++;
                    } else if (charge > 0) {
                        pre_trigger_total_positive_charge++;
                        pre_trigger_positive_pid_counts[pid]++;
                        if (chi2pid != 9999.0) {
                            pre_trigger_positive_chi2pid_valid_counts[pid]++;
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

                    // Second loop: Process positive hadrons only if a trigger electron is found
                    for (int i = 0; i < PART.getRows(); ++i) {
                        int pid = PART.getInt("pid", i);
                        int charge = PART.getByte("charge", i);
                        float chi2pid = PART.getFloat("chi2pid", i);
                        float px = PART.getFloat("px", i);
                        float py = PART.getFloat("py", i);
                        float pz = PART.getFloat("pz", i);
                        float p = std::sqrt(px * px + py * py + pz * pz);
                        float beta = PART.getFloat("beta", i);
                        int status = PART.getShort("status", i);
                        if (status == 0) continue;

                        bool isForward = (abs(status) / 2000 == 1);
                        bool isCentral = (abs(status) / 4000 == 1);

                        if (charge > 0) {
                            post_trigger_total_positive_charge++;

                            if (beta <= 0 || beta > 1.2) {
                                post_trigger_positive_beta_invalid++;
                            }
                            if (beta > 0 && beta <= 1.2 && p > 0) {
                                post_trigger_positive_beta_p_valid++;
                            }
                            if (beta >= 1 && beta <= 1.2) {
                                post_trigger_positive_beta_1_to_1p2++;
                            }

                            // Fill chi2pid histograms for positive pions
                            if (pid == 211 && chi2pid != 9999.0) {
                                if (chi2pid > 0) {
                                    post_trigger_pid_211_chi2pid_positive++;
                                } else if (chi2pid < 0) {
                                    post_trigger_pid_211_chi2pid_negative++;
                                } else {
                                    post_trigger_pid_211_chi2pid_zero++;
                                }
                                if (isForward) {
                                    hChi2PIDPositivePionsForward->Fill(chi2pid);
                                }
                                if (isCentral) {
                                    hChi2PIDPositivePionsCentral->Fill(chi2pid);
                                }
                            }

                            if (beta > 0 && beta < 1) {
                                float mass = (p / beta) * std::sqrt(1 - beta * beta);
                                hMassPositiveAll->Fill(mass);
                                post_trigger_positive_hadron_masses.push_back(mass);
                                post_trigger_positive_hadron_pids.push_back(pid);
                                post_trigger_positive_hadron_chi2pids.push_back(chi2pid);
                                post_trigger_positive_hadron_is_forward.push_back(isForward);
                                post_trigger_positive_hadron_is_central.push_back(isCentral);
                            }

                            if (isForward) {
                                post_trigger_positive_pid_counts_forward[pid]++;
                                if (chi2pid != 9999.0) {
                                    post_trigger_positive_chi2pid_valid_counts_forward[pid]++;
                                } else {
                                    post_trigger_positive_chi2pid_invalid_counts_forward[pid]++;
                                }
                            }
                            if (isCentral) {
                                post_trigger_positive_pid_counts_central[pid]++;
                                if (chi2pid != 9999.0) {
                                    post_trigger_positive_chi2pid_valid_counts_central[pid]++;
                                } else {
                                    post_trigger_positive_chi2pid_invalid_counts_central[pid]++;
                                }
                            }
                        } else if (charge < 0) {
                            post_trigger_total_negative_charge++;
                        }

                        // Count positive pions in events with a trigger electron
                        if (pid == 211) {
                            post_trigger_positive_pions_all_regions++;
                            if (abs(status) / 2000 == 1) {
                                post_trigger_positive_pions_forward++;
                            }
                            if (abs(status) / 4000 == 1) {
                                post_trigger_positive_pions_central++;
                            }
                        }
                    }
                }

                if (event_count % 1000 == 0) {
                    cout << "Processed " << event_count << " events" << endl;
                }
            }
        }
    }

    // Fit the mass histogram with a simple Gaussian
    TCanvas* canvasMassPositiveAll = new TCanvas("canvasMassPositiveAll", "Mass of All Positive Hadrons (Post-Trigger)", 800, 600);
    hMassPositiveAll->Draw("HIST");

    TF1* gaussFit = new TF1("gaussFit", "gaus", 0.08, 0.22); // Wider range to capture the pion peak
    gaussFit->SetParameters(4000, 0.14, 0.02);
    gaussFit->SetParNames("Amplitude", "Mean", "Sigma");

    hMassPositiveAll->Fit(gaussFit, "R");

    double pion_mass_mean = gaussFit->GetParameter(1);
    double pion_mass_sigma = gaussFit->GetParameter(2);
    double pion_mass_low = pion_mass_mean - 2 * pion_mass_sigma;
    double pion_mass_high = pion_mass_mean + 2 * pion_mass_sigma;
    cout << "Pion peak fit (Gaussian) - Mean: " << pion_mass_mean << " GeV/c^2, Sigma: " << pion_mass_sigma << " GeV/c^2" << endl;
    cout << "Pion mass selection range (Gaussian): [" << pion_mass_low << ", " << pion_mass_high << "] GeV/c^2" << endl;

    // Fit with Gaussian + linear background
    TCanvas* canvasMassWithBkg = new TCanvas("canvasMassWithBkg", "Mass of Positive Hadrons with Background Fit (Post-Trigger)", 800, 600);
    hMassPositiveAll->Draw("HIST");

    TF1* gaussPlusBkg = new TF1("gaussPlusBkg", "gaus(0) + pol1(3)", 0.05, 0.4); // Wider range to better estimate the background
    gaussPlusBkg->SetParameters(4000, 0.14, 0.02, 10, -10); // Reduced initial background parameters
    gaussPlusBkg->SetParNames("Amplitude", "Mean", "Sigma", "Bkg_Const", "Bkg_Slope");
    gaussPlusBkg->SetParLimits(3, -50, 50); // Constrain Bkg_Const
    gaussPlusBkg->SetParLimits(4, -100, 100); // Constrain Bkg_Slope

    hMassPositiveAll->Fit(gaussPlusBkg, "R");

    double pion_mass_mean_bkg = gaussPlusBkg->GetParameter(1);
    double pion_mass_sigma_bkg = gaussPlusBkg->GetParameter(2);
    double amplitude_bkg = gaussPlusBkg->GetParameter(0);
    double bkg_const = gaussPlusBkg->GetParameter(3);
    double bkg_slope = gaussPlusBkg->GetParameter(4);
    double pion_mass_low_bkg = pion_mass_mean_bkg - 2 * pion_mass_sigma_bkg;
    double pion_mass_high_bkg = pion_mass_mean_bkg + 2 * pion_mass_sigma_bkg;
    cout << "Pion peak fit (with background) - Mean: " << pion_mass_mean_bkg << " GeV/c^2, Sigma: " << pion_mass_sigma_bkg << " GeV/c^2" << endl;
    cout << "Background parameters - Const: " << bkg_const << ", Slope: " << bkg_slope << endl;
    cout << "Updated pion mass selection range (with background): [" << pion_mass_low_bkg << ", " << pion_mass_high_bkg << "] GeV/c^2" << endl;

    // Count particles in both mass ranges (overall and by region)
    for (size_t i = 0; i < post_trigger_positive_hadron_masses.size(); ++i) {
        float mass = post_trigger_positive_hadron_masses[i];
        int pid = post_trigger_positive_hadron_pids[i];
        bool isForward = post_trigger_positive_hadron_is_forward[i];
        bool isCentral = post_trigger_positive_hadron_is_central[i];

        // Gaussian-only range
        if (mass >= pion_mass_low && mass <= pion_mass_high) {
            post_trigger_pions_in_mass_range++;
            if (pid == 211) {
                post_trigger_pions_pid_211_in_mass_range++;
            }
        }

        // Gaussian + background range
        if (mass >= pion_mass_low_bkg && mass <= pion_mass_high_bkg) {
            post_trigger_pions_in_mass_range_updated++;
            if (pid == 211) {
                post_trigger_pions_pid_211_in_mass_range_updated++;
            }
            if (isForward) {
                post_trigger_pions_in_mass_range_updated_forward++;
                if (pid == 211) {
                    post_trigger_pions_pid_211_in_mass_range_updated_forward++;
                }
            }
            if (isCentral) {
                post_trigger_pions_in_mass_range_updated_central++;
                if (pid == 211) {
                    post_trigger_pions_pid_211_in_mass_range_updated_central++;
                }
            }
        }

        // Count pid == 211 outside the updated mass range
        if (!(mass >= pion_mass_low_bkg && mass <= pion_mass_high_bkg)) {
            if (pid == 211) {
                post_trigger_pions_pid_211_outside_mass_range++;
            }
        }
    }

    // Estimate contamination for Gaussian + background fit (overall)
    double bkg_integral_bkg = (bkg_const * (pion_mass_high_bkg - pion_mass_low_bkg)) + 
                              (0.5 * bkg_slope * (pion_mass_high_bkg * pion_mass_high_bkg - pion_mass_low_bkg * pion_mass_low_bkg));
    double gauss_integral_bkg = amplitude_bkg * sqrt(2 * TMath::Pi()) * pion_mass_sigma_bkg * 
                                (TMath::Erf((pion_mass_high_bkg - pion_mass_mean_bkg) / (sqrt(2) * pion_mass_sigma_bkg)) - 
                                 TMath::Erf((pion_mass_low_bkg - pion_mass_mean_bkg) / (sqrt(2) * pion_mass_sigma_bkg))) / 2;
    double total_integral_bkg = bkg_integral_bkg + gauss_integral_bkg;
    double contamination_fraction_bkg = (total_integral_bkg > 0) ? (bkg_integral_bkg / total_integral_bkg) * 100.0 : 0.0;
    cout << "\nContamination estimate for Gaussian + background fit (overall, before chi2pid cut):" << endl;
    cout << "Background integral in updated mass range: " << bkg_integral_bkg << endl;
    cout << "Gaussian integral in updated mass range: " << gauss_integral_bkg << endl;
    cout << "Total integral (background + Gaussian) in updated range: " << total_integral_bkg << endl;
    cout << "Estimated contamination fraction in pion sample (with background): " << contamination_fraction_bkg << "%" << endl;

    // Estimate contamination by region (Forward and Central)
    double bkg_integral_bkg_forward = bkg_integral_bkg * (post_trigger_pions_in_mass_range_updated_forward / (double)post_trigger_pions_in_mass_range_updated);
    double bkg_integral_bkg_central = bkg_integral_bkg * (post_trigger_pions_in_mass_range_updated_central / (double)post_trigger_pions_in_mass_range_updated);
    double total_integral_bkg_forward = total_integral_bkg * (post_trigger_pions_in_mass_range_updated_forward / (double)post_trigger_pions_in_mass_range_updated);
    double total_integral_bkg_central = total_integral_bkg * (post_trigger_pions_in_mass_range_updated_central / (double)post_trigger_pions_in_mass_range_updated);
    double contamination_fraction_bkg_forward = (total_integral_bkg_forward > 0) ? (bkg_integral_bkg_forward / total_integral_bkg_forward) * 100.0 : 0.0;
    double contamination_fraction_bkg_central = (total_integral_bkg_central > 0) ? (bkg_integral_bkg_central / total_integral_bkg_central) * 100.0 : 0.0;
    cout << "\nContamination estimate for Forward region (before chi2pid cut):" << endl;
    cout << "Background integral in updated mass range (Forward): " << bkg_integral_bkg_forward << endl;
    cout << "Total integral in updated mass range (Forward): " << total_integral_bkg_forward << endl;
    cout << "Total particles in updated mass range (Forward): " << post_trigger_pions_in_mass_range_updated_forward << endl;
    cout << "Particles with PID 211 in updated mass range (Forward): " << post_trigger_pions_pid_211_in_mass_range_updated_forward << endl;
    cout << "Estimated contamination fraction (Forward): " << contamination_fraction_bkg_forward << "%" << endl;
    cout << "\nContamination estimate for Central region (before chi2pid cut):" << endl;
    cout << "Background integral in updated mass range (Central): " << bkg_integral_bkg_central << endl;
    cout << "Total integral in updated mass range (Central): " << total_integral_bkg_central << endl;
    cout << "Total particles in updated mass range (Central): " << post_trigger_pions_in_mass_range_updated_central << endl;
    cout << "Particles with PID 211 in updated mass range (Central): " << post_trigger_pions_pid_211_in_mass_range_updated_central << endl;
    cout << "Estimated contamination fraction (Central): " << contamination_fraction_bkg_central << "%" << endl;

    // Estimate contamination for Gaussian-only fit (overall)
    double bkg_integral_gauss = (bkg_const * (pion_mass_high - pion_mass_low)) + 
                                (0.5 * bkg_slope * (pion_mass_high * pion_mass_high - pion_mass_low * pion_mass_low));
    double gauss_integral_gauss = amplitude_bkg * sqrt(2 * TMath::Pi()) * pion_mass_sigma_bkg * 
                                  (TMath::Erf((pion_mass_high - pion_mass_mean_bkg) / (sqrt(2) * pion_mass_sigma_bkg)) - 
                                   TMath::Erf((pion_mass_low - pion_mass_mean_bkg) / (sqrt(2) * pion_mass_sigma_bkg))) / 2;
    double total_integral_gauss = bkg_integral_gauss + gauss_integral_gauss;
    double contamination_fraction_gauss = (total_integral_gauss > 0) ? (bkg_integral_gauss / total_integral_gauss) * 100.0 : 0.0;
    cout << "\nContamination estimate for Gaussian-only fit (overall):" << endl;
    cout << "Background integral in Gaussian-only mass range: " << bkg_integral_gauss << endl;
    cout << "Gaussian integral in Gaussian-only mass range: " << gauss_integral_gauss << endl;
    cout << "Total integral (background + Gaussian) in Gaussian-only range: " << total_integral_gauss << endl;
    cout << "Estimated contamination fraction in pion sample (Gaussian-only): " << contamination_fraction_gauss << "%" << endl;

    // Fit the chi2pid histogram for positive pions (Forward region)
    TCanvas* canvasChi2PIDForward = new TCanvas("canvasChi2PIDForward", "Chi2PID for Positive Pions (PID 211, Forward, Post-Trigger)", 800, 600);
    hChi2PIDPositivePionsForward->Draw("HIST");

    TF1* chi2pidFitForward = new TF1("chi2pidFitForward", "gaus", -2, 2);
    chi2pidFitForward->SetParameters(1000, 0, 1);
    chi2pidFitForward->SetParNames("Amplitude", "Mean", "Sigma");

    hChi2PIDPositivePionsForward->Fit(chi2pidFitForward, "R");

    double chi2pid_mean_forward = chi2pidFitForward->GetParameter(1);
    double chi2pid_sigma_forward = chi2pidFitForward->GetParameter(2);
    double chi2pid_low_cut_forward = chi2pid_mean_forward - 3.0;
    double chi2pid_high_cut_forward = chi2pid_mean_forward + 3.0;
    cout << "\nChi2PID fit for positive pions (PID 211, Forward, Post-Trigger):" << endl;
    cout << "Mean: " << chi2pid_mean_forward << endl;
    cout << "Sigma: " << chi2pid_sigma_forward << endl;
    cout << "Cut range (|chi2pid - mean| < 3): [" << chi2pid_low_cut_forward << ", " << chi2pid_high_cut_forward << "]" << endl;

    chi2pidFitForward->SetLineColor(kRed);
    chi2pidFitForward->Draw("SAME");
    canvasChi2PIDForward->Print("output/chi2pid_positive_pions_forward.pdf");
    cout << "Chi2PID histogram for positive pions (Forward, Post-Trigger) saved to output/chi2pid_positive_pions_forward.pdf" << endl;

    // Fit the chi2pid histogram for positive pions (Central region)
    TCanvas* canvasChi2PIDCentral = new TCanvas("canvasChi2PIDCentral", "Chi2PID for Positive Pions (PID 211, Central, Post-Trigger)", 800, 600);
    hChi2PIDPositivePionsCentral->Draw("HIST");

    TF1* chi2pidFitCentral = new TF1("chi2pidFitCentral", "gaus", -2, 2);
    chi2pidFitCentral->SetParameters(1000, 0, 1);
    chi2pidFitCentral->SetParNames("Amplitude", "Mean", "Sigma");

    hChi2PIDPositivePionsCentral->Fit(chi2pidFitCentral, "R");

    double chi2pid_mean_central = chi2pidFitCentral->GetParameter(1);
    double chi2pid_sigma_central = chi2pidFitCentral->GetParameter(2);
    double chi2pid_low_cut_central = chi2pid_mean_central - 3.0;
    double chi2pid_high_cut_central = chi2pid_mean_central + 3.0;
    cout << "\nChi2PID fit for positive pions (PID 211, Central, Post-Trigger):" << endl;
    cout << "Mean: " << chi2pid_mean_central << endl;
    cout << "Sigma: " << chi2pid_sigma_central << endl;
    cout << "Cut range (|chi2pid - mean| < 3): [" << chi2pid_low_cut_central << ", " << chi2pid_high_cut_central << "]" << endl;

    chi2pidFitCentral->SetLineColor(kRed);
    chi2pidFitCentral->Draw("SAME");
    canvasChi2PIDCentral->Print("output/chi2pid_positive_pions_central.pdf");
    cout << "Chi2PID histogram for positive pions (Central, Post-Trigger) saved to output/chi2pid_positive_pions_central.pdf" << endl;

    // Apply chi2pid cut and recalculate contamination (overall and by region)
    int post_trigger_pions_pid_211_in_mass_range_after_chi2pid = 0;
    int post_trigger_pions_pid_211_in_mass_range_after_chi2pid_forward = 0;
    int post_trigger_pions_pid_211_in_mass_range_after_chi2pid_central = 0;

    for (size_t i = 0; i < post_trigger_positive_hadron_masses.size(); ++i) {
        float mass = post_trigger_positive_hadron_masses[i];
        int pid = post_trigger_positive_hadron_pids[i];
        float chi2pid = post_trigger_positive_hadron_chi2pids[i];
        bool isForward = post_trigger_positive_hadron_is_forward[i];
        bool isCentral = post_trigger_positive_hadron_is_central[i];

        if (mass >= pion_mass_low_bkg && mass <= pion_mass_high_bkg && pid == 211 && chi2pid != 9999.0) {
            // Forward region
            if (isForward && chi2pid >= chi2pid_low_cut_forward && chi2pid <= chi2pid_high_cut_forward) {
                post_trigger_pions_pid_211_in_mass_range_after_chi2pid_forward++;
                post_trigger_pions_pid_211_in_mass_range_after_chi2pid++; // Increment overall counter
            }

            // Central region
            if (isCentral && chi2pid >= chi2pid_low_cut_central && chi2pid <= chi2pid_high_cut_central) {
                post_trigger_pions_pid_211_in_mass_range_after_chi2pid_central++;
                // Avoid double-counting if the particle is in both regions
                if (!isForward || !(chi2pid >= chi2pid_low_cut_forward && chi2pid <= chi2pid_high_cut_forward)) {
                    post_trigger_pions_pid_211_in_mass_range_after_chi2pid++;
                }
            }
        }
    }

    // Recalculate contamination after chi2pid cut (overall)
    double total_after_chi2pid = post_trigger_pions_in_mass_range_updated - (post_trigger_pions_pid_211_in_mass_range_updated - post_trigger_pions_pid_211_in_mass_range_after_chi2pid);
    double contamination_fraction_after_chi2pid = (total_after_chi2pid > 0) ? (bkg_integral_bkg / total_after_chi2pid) * 100.0 : 0.0;
    cout << "\nContamination estimate after chi2pid cut (overall):" << endl;
    cout << "Particles with PID 211 in updated mass range after chi2pid cut: " << post_trigger_pions_pid_211_in_mass_range_after_chi2pid << endl;
    cout << "Total particles in updated mass range after chi2pid cut (adjusted): " << total_after_chi2pid << endl;
    cout << "Estimated contamination fraction after chi2pid cut: " << contamination_fraction_after_chi2pid << "%" << endl;

    // Recalculate contamination after chi2pid cut (by region)
    double total_after_chi2pid_forward = post_trigger_pions_in_mass_range_updated_forward - (post_trigger_pions_pid_211_in_mass_range_updated_forward - post_trigger_pions_pid_211_in_mass_range_after_chi2pid_forward);
    double total_after_chi2pid_central = post_trigger_pions_in_mass_range_updated_central - (post_trigger_pions_pid_211_in_mass_range_updated_central - post_trigger_pions_pid_211_in_mass_range_after_chi2pid_central);
    double contamination_fraction_after_chi2pid_forward = (total_after_chi2pid_forward > 0) ? (bkg_integral_bkg_forward / total_after_chi2pid_forward) * 100.0 : 0.0;
    double contamination_fraction_after_chi2pid_central = (total_after_chi2pid_central > 0) ? (bkg_integral_bkg_central / total_after_chi2pid_central) * 100.0 : 0.0;
    cout << "\nContamination estimate after chi2pid cut (Forward):" << endl;
    cout << "Particles with PID 211 in updated mass range after chi2pid cut (Forward): " << post_trigger_pions_pid_211_in_mass_range_after_chi2pid_forward << endl;
    cout << "Total particles in updated mass range after chi2pid cut (Forward, adjusted): " << total_after_chi2pid_forward << endl;
    cout << "Estimated contamination fraction after chi2pid cut (Forward): " << contamination_fraction_after_chi2pid_forward << "%" << endl;
    cout << "\nContamination estimate after chi2pid cut (Central):" << endl;
    cout << "Particles with PID 211 in updated mass range after chi2pid cut (Central): " << post_trigger_pions_pid_211_in_mass_range_after_chi2pid_central << endl;
    cout << "Total particles in updated mass range after chi2pid cut (Central, adjusted): " << total_after_chi2pid_central << endl;
    cout << "Estimated contamination fraction after chi2pid cut (Central): " << contamination_fraction_after_chi2pid_central << "%" << endl;

    // Save the fits
    gaussFit->SetLineColor(kRed);
    gaussFit->Draw("SAME");
    canvasMassPositiveAll->Print("output/mass_positive_all_with_fit.pdf");
    cout << "Mass histogram with pion peak fit (Gaussian) saved to output/mass_positive_all_with_fit.pdf" << endl;

    gaussPlusBkg->SetLineColor(kRed);
    gaussPlusBkg->Draw("SAME");
    canvasMassWithBkg->Print("output/mass_positive_all_with_fit_and_bkg.pdf");
    cout << "Mass histogram with fit (including background) saved to output/mass_positive_all_with_fit_and_bkg.pdf" << endl;

    timer.Stop();
    double execution_time = timer.RealTime();

    // Validation checks
    int sum_pid_counts = 0;
    for (const auto& pid_count : pre_trigger_pid_counts) {
        sum_pid_counts += pid_count.second;
    }
    bool pid_count_matches = (sum_pid_counts == pre_trigger_total_particles);
    string pid_validation_message = pid_count_matches ? "Yes" : "No (Possible data inconsistency)";

    int sum_positive_pions_regions = post_trigger_positive_pions_forward + post_trigger_positive_pions_central;
    bool positive_pions_sum_matches = (sum_positive_pions_regions == post_trigger_positive_pions_all_regions);
    string positive_pions_validation_message = positive_pions_sum_matches ? "Yes" : "No (Possible data inconsistency)";

    int sum_valid_chi2pid = 0;
    int sum_invalid_chi2pid = 0;
    for (const auto& pid_count : pre_trigger_pid_counts) {
        int pid = pid_count.first;
        sum_valid_chi2pid += pre_trigger_chi2pid_valid_counts[pid];
        sum_invalid_chi2pid += pre_trigger_chi2pid_invalid_counts[pid];
    }
    bool chi2pid_sum_matches = (sum_valid_chi2pid + sum_invalid_chi2pid == sum_pid_counts);
    string chi2pid_validation_message = chi2pid_sum_matches ? "Yes" : "No (Possible data inconsistency)";

    int mass_histogram_entries = hMassPositiveAll->GetEntries();
    cout << "Entries in hMassPositiveAll: " << mass_histogram_entries << endl;
    cout << "Expected entries (post_trigger_positive_beta_p_valid - post_trigger_positive_beta_1_to_1p2): " 
         << post_trigger_positive_beta_p_valid - post_trigger_positive_beta_1_to_1p2 << endl;

    // Write to CSV
    ofstream csvFile("output/particle_stats.csv");
    if (!csvFile.is_open()) {
        cerr << "Error: Could not open output/particle_stats.csv for writing!" << endl;
        return 1;
    }

    csvFile << "Statistic,Value\n";
    csvFile << "Total events processed," << total_events_processed << "\n";
    csvFile << "Pre-trigger total particles," << pre_trigger_total_particles << "\n";
    csvFile << "Pre-trigger total neutral particles (charge == 0)," << pre_trigger_total_neutral_particles << "\n";
    csvFile << "Pre-trigger total positive charged particles (charge > 0)," << pre_trigger_total_positive_charge << "\n";
    csvFile << "Pre-trigger total negative charged particles (charge < 0)," << pre_trigger_total_negative_charge << "\n";
    csvFile << "Post-trigger total positive charged particles," << post_trigger_total_positive_charge << "\n";
    csvFile << "Post-trigger total negative charged particles," << post_trigger_total_negative_charge << "\n";
    csvFile << "Post-trigger positive particles with invalid beta (beta <= 0 or beta > 1.2)," << post_trigger_positive_beta_invalid << "\n";
    csvFile << "Post-trigger positive particles with valid beta and p (plotted in beta vs. p)," << post_trigger_positive_beta_p_valid << "\n";
    csvFile << "Post-trigger positive particles with 1 <= beta <= 1.2 (excluded from mass histogram)," << post_trigger_positive_beta_1_to_1p2 << "\n";
    csvFile << "Entries in mass histogram (hMassPositiveAll)," << mass_histogram_entries << "\n";
    csvFile << "Pion peak mean (Gaussian) (GeV/c^2)," << pion_mass_mean << "\n";
    csvFile << "Pion peak sigma (Gaussian) (GeV/c^2)," << pion_mass_sigma << "\n";
    csvFile << "Pion mass selection range low (Gaussian) (GeV/c^2)," << pion_mass_low << "\n";
    csvFile << "Pion mass selection range high (Gaussian) (GeV/c^2)," << pion_mass_high << "\n";
    csvFile << "Post-trigger number of particles in pion mass range (Gaussian)," << post_trigger_pions_in_mass_range << "\n";
    csvFile << "Post-trigger number of particles with pid == 211 in Gaussian-only mass range," << post_trigger_pions_pid_211_in_mass_range << "\n";
    csvFile << "Estimated contamination fraction in pion sample (Gaussian-only) (%)," << contamination_fraction_gauss << "\n";
    csvFile << "Pion peak mean (with background) (GeV/c^2)," << pion_mass_mean_bkg << "\n";
    csvFile << "Pion peak sigma (with background) (GeV/c^2)," << pion_mass_sigma_bkg << "\n";
    csvFile << "Pion mass selection range low (with background) (GeV/c^2)," << pion_mass_low_bkg << "\n";
    csvFile << "Pion mass selection range high (with background) (GeV/c^2)," << pion_mass_high_bkg << "\n";
    csvFile << "Post-trigger number of particles in updated pion mass range (with background)," << post_trigger_pions_in_mass_range_updated << "\n";
    csvFile << "Post-trigger number of particles with pid == 211 in updated mass range," << post_trigger_pions_pid_211_in_mass_range_updated << "\n";
    csvFile << "Estimated contamination fraction in pion sample (with background) (%)," << contamination_fraction_bkg << "\n";
    csvFile << "Post-trigger number of particles with pid == 211 outside the updated mass range," << post_trigger_pions_pid_211_outside_mass_range << "\n";
    csvFile << "Post-trigger number of particles in updated pion mass range (Forward)," << post_trigger_pions_in_mass_range_updated_forward << "\n";
    csvFile << "Post-trigger number of particles with pid == 211 in updated mass range (Forward)," << post_trigger_pions_pid_211_in_mass_range_updated_forward << "\n";
    csvFile << "Estimated contamination fraction in pion sample (Forward) (%)," << contamination_fraction_bkg_forward << "\n";
    csvFile << "Post-trigger number of particles in updated pion mass range (Central)," << post_trigger_pions_in_mass_range_updated_central << "\n";
    csvFile << "Post-trigger number of particles with pid == 211 in updated mass range (Central)," << post_trigger_pions_pid_211_in_mass_range_updated_central << "\n";
    csvFile << "Estimated contamination fraction in pion sample (Central) (%)," << contamination_fraction_bkg_central << "\n";
    csvFile << "Post-trigger number of particles with pid == 211 and chi2pid > 0," << post_trigger_pid_211_chi2pid_positive << "\n";
    csvFile << "Post-trigger number of particles with pid == 211 and chi2pid < 0," << post_trigger_pid_211_chi2pid_negative << "\n";
    csvFile << "Post-trigger number of particles with pid == 211 and chi2pid == 0," << post_trigger_pid_211_chi2pid_zero << "\n";
    csvFile << "Chi2PID mean for positive pions (Forward, Post-Trigger)," << chi2pid_mean_forward << "\n";
    csvFile << "Chi2PID sigma for positive pions (Forward, Post-Trigger)," << chi2pid_sigma_forward << "\n";
    csvFile << "Chi2PID cut range low (Forward)," << chi2pid_low_cut_forward << "\n";
    csvFile << "Chi2PID cut range high (Forward)," << chi2pid_high_cut_forward << "\n";
    csvFile << "Chi2PID mean for positive pions (Central, Post-Trigger)," << chi2pid_mean_central << "\n";
    csvFile << "Chi2PID sigma for positive pions (Central, Post-Trigger)," << chi2pid_sigma_central << "\n";
    csvFile << "Chi2PID cut range low (Central)," << chi2pid_low_cut_central << "\n";
    csvFile << "Chi2PID cut range high (Central)," << chi2pid_high_cut_central << "\n";
    csvFile << "Post-trigger number of particles with pid == 211 in updated mass range after chi2pid cut (overall)," << post_trigger_pions_pid_211_in_mass_range_after_chi2pid << "\n";
    csvFile << "Estimated contamination fraction after chi2pid cut (overall) (%)," << contamination_fraction_after_chi2pid << "\n";
    csvFile << "Post-trigger number of particles with pid == 211 in updated mass range after chi2pid cut (Forward)," << post_trigger_pions_pid_211_in_mass_range_after_chi2pid_forward << "\n";
    csvFile << "Estimated contamination fraction after chi2pid cut (Forward) (%)," << contamination_fraction_after_chi2pid_forward << "\n";
    csvFile << "Post-trigger number of particles with pid == 211 in updated mass range after chi2pid cut (Central)," << post_trigger_pions_pid_211_in_mass_range_after_chi2pid_central << "\n";
    csvFile << "Estimated contamination fraction after chi2pid cut (Central) (%)," << contamination_fraction_after_chi2pid_central << "\n";
    csvFile << "Pre-trigger total electrons (pid == 11)," << pre_trigger_electron_count << "\n";
    csvFile << "Trigger electron events," << trigger_elec_events << "\n";
    csvFile << "Execution time (seconds)," << execution_time << "\n";

    csvFile << "\nPre-Trigger PID Counts (REC::Particle),\n";
    csvFile << "PID,Count\n";
    for (const auto& pid_count : pre_trigger_pid_counts) {
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

    csvFile << "\nPost-Trigger Positive Pion Statistics (Events with Trigger Electron),\n";
    csvFile << "Statistic,Count\n";
    csvFile << "Pre-trigger total positive pions (pid == 211) across all events," << pre_trigger_pid_counts[211] << "\n";
    csvFile << "Post-trigger positive pions (pid == 211) in events with trigger electron (all regions)," << post_trigger_positive_pions_all_regions << "\n";
    csvFile << "Post-trigger positive pions in Forward region (abs(status) / 2000 == 1)," << post_trigger_positive_pions_forward << "\n";
    csvFile << "Post-trigger positive pions in Central region (abs(status) / 4000 == 1)," << post_trigger_positive_pions_central << "\n";
    csvFile << "Sum of Forward and Central regions," << sum_positive_pions_regions << "\n";
    csvFile << "Sum matches post-trigger total positive pions in events with trigger electron?," << positive_pions_validation_message << "\n";

    csvFile << "\nPre-Trigger PID Counts with Valid Chi2PID (REC::Particle),\n";
    csvFile << "PID,Count\n";
    for (const auto& pid_count : pre_trigger_pid_counts) {
        int pid = pid_count.first;
        int valid_count = pre_trigger_chi2pid_valid_counts[pid];
        csvFile << pid << "," << valid_count << "\n";
    }

    csvFile << "\nPre-Trigger PID Counts with Invalid Chi2PID (REC::Particle),\n";
    csvFile << "PID,Count\n";
    for (const auto& pid_count : pre_trigger_pid_counts) {
        int pid = pid_count.first;
        int invalid_count = pre_trigger_chi2pid_invalid_counts[pid];
        csvFile << pid << "," << invalid_count << "\n";
    }

    csvFile << "\nPre-Trigger Chi2PID Count Validation,\n";
    csvFile << "Sum of particles with valid Chi2PID," << sum_valid_chi2pid << "\n";
    csvFile << "Sum of particles with invalid Chi2PID," << sum_invalid_chi2pid << "\n";
    csvFile << "Total (valid + invalid)," << (sum_valid_chi2pid + sum_invalid_chi2pid) << "\n";
    csvFile << "Matches sum of pre-trigger PID counts?," << chi2pid_validation_message << "\n";

    // Add pre-trigger positive charged particle chi2pid statistics to CSV
    csvFile << "\nPre-Trigger Positive Charged Particles by PID and Chi2PID Validity,\n";
    csvFile << "PID,Total Positive,Valid Chi2PID,Invalid Chi2PID\n";
    for (const auto& pid : pids_to_visualize) {
        int total = pre_trigger_positive_pid_counts[pid];
        int valid = pre_trigger_positive_chi2pid_valid_counts[pid];
        int invalid = pre_trigger_positive_chi2pid_invalid_counts[pid];
        csvFile << pid << "," << total << "," << valid << "," << invalid << "\n";
    }

    // Add post-trigger region-separated positive charged particle chi2pid statistics to CSV
    csvFile << "\nPost-Trigger Positive Charged Particles by PID and Chi2PID Validity (Forward Region),\n";
    csvFile << "PID,Total Positive,Valid Chi2PID,Invalid Chi2PID\n";
    for (const auto& pid : pids_to_visualize) {
        int total = post_trigger_positive_pid_counts_forward[pid];
        int valid = post_trigger_positive_chi2pid_valid_counts_forward[pid];
        int invalid = post_trigger_positive_chi2pid_invalid_counts_forward[pid];
        csvFile << pid << "," << total << "," << valid << "," << invalid << "\n";
    }

    csvFile << "\nPost-Trigger Positive Charged Particles by PID and Chi2PID Validity (Central Region),\n";
    csvFile << "PID,Total Positive,Valid Chi2PID,Invalid Chi2PID\n";
    for (const auto& pid : pids_to_visualize) {
        int total = post_trigger_positive_pid_counts_central[pid];
        int valid = post_trigger_positive_chi2pid_valid_counts_central[pid];
        int invalid = post_trigger_positive_chi2pid_invalid_counts_central[pid];
        csvFile << pid << "," << total << "," << valid << "," << invalid << "\n";
    }

    csvFile.close();
    cout << "Statistics saved to output/particle_stats.csv" << endl;

    // Print summary
    cout << "Total events processed: " << total_events_processed << endl;
    cout << "Pre-trigger total particles: " << pre_trigger_total_particles << endl;
    cout << "Pre-trigger total neutral particles (charge == 0): " << pre_trigger_total_neutral_particles << endl;
    cout << "Pre-trigger total positive charged particles (charge > 0): " << pre_trigger_total_positive_charge << endl;
    cout << "Pre-trigger total negative charged particles (charge < 0): " << pre_trigger_total_negative_charge << endl;
    cout << "Post-trigger total positive charged particles: " << post_trigger_total_positive_charge << endl;
    cout << "Post-trigger total negative charged particles: " << post_trigger_total_negative_charge << endl;
    cout << "Post-trigger positive particles with invalid beta (beta <= 0 or beta > 1.2): " << post_trigger_positive_beta_invalid << endl;
    cout << "Post-trigger positive particles with valid beta and p: " << post_trigger_positive_beta_p_valid << endl;
    cout << "Post-trigger positive particles with 1 <= beta <= 1.2 (excluded from mass histogram): " << post_trigger_positive_beta_1_to_1p2 << endl;
    cout << "Pre-trigger total electrons (pid == 11): " << pre_trigger_electron_count << endl;
    cout << "Trigger electron events: " << trigger_elec_events << endl;

    cout << "\nPre-Trigger PID Counts (REC::Particle):" << endl;
    for (const auto& pid_count : pre_trigger_pid_counts) {
        cout << "PID " << pid_count.first << ": " << pid_count.second << " particles" << endl;
    }

    cout << "\nPre-Trigger PID Count Validation:" << endl;
    cout << "Sum of pre-trigger PID counts: " << sum_pid_counts << endl;
    cout << "Matches pre-trigger total particles? " << pid_validation_message << endl;

    cout << "\nTrigger Electron Cut Statistics:" << endl;
    cout << "Electrons with pid == 11: " << trigger_elec_electrons_pid_11 << endl;
    cout << "Electrons with pid == 11 and status < 0: " << trigger_elec_status_negative << endl;
    cout << "Electrons with pid == 11, status < 0, and vz_tele in [-20, 5]: " << trigger_elec_vz_range << endl;
    cout << "Electrons with pid == 11, status < 0, vz_tele in [-20, 5], and abs(status) / 1000 == 2: " << trigger_elec_final << endl;

    cout << "\nPost-Trigger Positive Pion Statistics (Events with Trigger Electron):" << endl;
    cout << "Pre-trigger total positive pions (pid == 211) across all events: " << pre_trigger_pid_counts[211] << endl;
    cout << "Post-trigger positive pions (pid == 211) in events with trigger electron (all regions): " << post_trigger_positive_pions_all_regions << endl;
    cout << "Post-trigger positive pions in Forward region (abs(status) / 2000 == 1): " << post_trigger_positive_pions_forward << endl;
    cout << "Post-trigger positive pions in Central region (abs(status) / 4000 == 1): " << post_trigger_positive_pions_central << endl;
    cout << "Sum of Forward and Central regions: " << sum_positive_pions_regions << endl;
    cout << "Sum matches post-trigger total positive pions in events with trigger electron? " << positive_pions_validation_message << endl;

    cout << "\nPost-Trigger Chi2PID Diagnostics for pid == 211:" << endl;
    cout << "Post-trigger number of particles with pid == 211 and chi2pid > 0: " << post_trigger_pid_211_chi2pid_positive << endl;
    cout << "Post-trigger number of particles with pid == 211 and chi2pid < 0: " << post_trigger_pid_211_chi2pid_negative << endl;
    cout << "Post-trigger number of particles with pid == 211 and chi2pid == 0: " << post_trigger_pid_211_chi2pid_zero << endl;

    cout << "\nPre-Trigger PID Counts with Valid Chi2PID (REC::Particle):" << endl;
    for (const auto& pid_count : pre_trigger_pid_counts) {
        int pid = pid_count.first;
        int valid_count = pre_trigger_chi2pid_valid_counts[pid];
        cout << "PID " << pid << ": " << valid_count << " particles" << endl;
    }

    cout << "\nPre-Trigger PID Counts with Invalid Chi2PID (REC::Particle):" << endl;
    for (const auto& pid_count : pre_trigger_pid_counts) {
        int pid = pid_count.first;
        int invalid_count = pre_trigger_chi2pid_invalid_counts[pid];
        cout << "PID " << pid << ": " << invalid_count << " particles" << endl;
    }

    cout << "\nPre-Trigger Chi2PID Count Validation:" << endl;
    cout << "Sum of pre-trigger particles with valid Chi2PID: " << sum_valid_chi2pid << endl;
    cout << "Sum of pre-trigger particles with invalid Chi2PID: " << sum_invalid_chi2pid << endl;
    cout << "Total (valid + invalid): " << (sum_valid_chi2pid + sum_invalid_chi2pid) << endl;
    cout << "Matches sum of pre-trigger PID counts? " << chi2pid_validation_message << endl;

    cout << "\nPre-Trigger Positive Charged Particles by PID and Chi2PID Validity:" << endl;
    cout << "PID\tTotal Positive\tValid Chi2PID\tInvalid Chi2PID\n";
    for (const auto& pid : pids_to_visualize) {
        int total = pre_trigger_positive_pid_counts[pid];
        int valid = pre_trigger_positive_chi2pid_valid_counts[pid];
        int invalid = pre_trigger_positive_chi2pid_invalid_counts[pid];
        cout << pid << "\t" << total << "\t\t" << valid << "\t\t" << invalid << endl;
    }

    cout << "\nPost-Trigger Positive Charged Particles by PID and Chi2PID Validity (Forward Region):" << endl;
    cout << "PID\tTotal Positive\tValid Chi2PID\tInvalid Chi2PID\n";
    for (const auto& pid : pids_to_visualize) {
        int total = post_trigger_positive_pid_counts_forward[pid];
        int valid = post_trigger_positive_chi2pid_valid_counts_forward[pid];
        int invalid = post_trigger_positive_chi2pid_invalid_counts_forward[pid];
        cout << pid << "\t" << total << "\t\t" << valid << "\t\t" << invalid << endl;
    }

    cout << "\nPost-Trigger Positive Charged Particles by PID and Chi2PID Validity (Central Region):" << endl;
    cout << "PID\tTotal Positive\tValid Chi2PID\tInvalid Chi2PID\n";
    for (const auto& pid : pids_to_visualize) {
        int total = post_trigger_positive_pid_counts_central[pid];
        int valid = post_trigger_positive_chi2pid_valid_counts_central[pid];
        int invalid = post_trigger_positive_chi2pid_invalid_counts_central[pid];
        cout << pid << "\t" << total << "\t\t" << valid << "\t\t" << invalid << endl;
    }

    cout << "Execution time: " << execution_time << " seconds" << endl;

    TCanvas* canvas = new TCanvas("canvas", "Electron Momentum", 800, 600);
    hElectronMomentum->Draw("HIST");
    canvas->Print("output/electron_momentum.pdf");
    cout << "Electron momentum histogram saved to output/electron_momentum.pdf" << endl;

    // Visualization 1: Stacked Bar Chart for All Particles (Neutral, Positive, Negative)
    TCanvas* canvasChargeBreakdown = new TCanvas("canvasChargeBreakdown", "Breakdown of Particles by Charge", 800, 600);
    THStack* hsCharge = new THStack("hsCharge", "Breakdown of Particles by Charge;Category;Counts");

    TH1F* hNeutral = new TH1F("hNeutral", "Neutral Particles", 1, 0, 1);
    hNeutral->SetBinContent(1, pre_trigger_total_neutral_particles);
    hNeutral->SetFillColor(kGray);
    hNeutral->SetLineColor(kBlack);

    TH1F* hPositive = new TH1F("hPositive", "Positive Particles", 1, 0, 1);
    hPositive->SetBinContent(1, post_trigger_total_positive_charge);
    hPositive->SetFillColor(kRed);
    hPositive->SetLineColor(kBlack);

    TH1F* hNegative = new TH1F("hNegative", "Negative Particles", 1, 0, 1);
    hNegative->SetBinContent(1, post_trigger_total_negative_charge);
    hNegative->SetFillColor(kBlue);
    hNegative->SetLineColor(kBlack);

    hsCharge->Add(hNeutral);
    hsCharge->Add(hPositive);
    hsCharge->Add(hNegative);

    hsCharge->Draw("HIST");
    hsCharge->GetXaxis()->SetBinLabel(1, "All Particles");

    TLegend* legendCharge = new TLegend(0.7, 0.7, 0.9, 0.9);
    legendCharge->AddEntry(hNeutral, Form("Neutral: %d", pre_trigger_total_neutral_particles), "f");
    legendCharge->AddEntry(hPositive, Form("Positive: %d", post_trigger_total_positive_charge), "f");
    legendCharge->AddEntry(hNegative, Form("Negative: %d", post_trigger_total_negative_charge), "f");
    legendCharge->Draw();

    canvasChargeBreakdown->Print("output/charge_breakdown_stacked.pdf");
    cout << "Stacked bar chart for charge breakdown saved to output/charge_breakdown_stacked.pdf" << endl;

    // Visualization 2: Grouped Bar Chart for Positive Particles Across Beta Ranges
    TCanvas* canvasBetaBreakdown = new TCanvas("canvasBetaBreakdown", "Breakdown of Positive Particles by Beta Range (Post-Trigger)", 800, 600);

    TH1F* hBetaInvalid = new TH1F("hBetaInvalid", "Invalid Beta", 3, 0, 3);
    TH1F* hBetaValidUsed = new TH1F("hBetaValidUsed", "Valid Beta (Used in Mass)", 3, 0, 3);
    TH1F* hBetaValidExcluded = new TH1F("hBetaValidExcluded", "Valid Beta (Excluded)", 3, 0, 3);

    hBetaInvalid->SetBinContent(1, post_trigger_positive_beta_invalid);
    hBetaValidUsed->SetBinContent(2, post_trigger_positive_beta_p_valid - post_trigger_positive_beta_1_to_1p2);
    hBetaValidExcluded->SetBinContent(3, post_trigger_positive_beta_1_to_1p2);

    hBetaInvalid->SetFillColor(kOrange);
    hBetaInvalid->SetLineColor(kBlack);
    hBetaValidUsed->SetFillColor(kGreen);
    hBetaValidUsed->SetLineColor(kBlack);
    hBetaValidExcluded->SetFillColor(kCyan);
    hBetaValidExcluded->SetLineColor(kBlack);

    hBetaInvalid->GetXaxis()->SetBinLabel(1, "Invalid Beta (beta <= 0 or > 1.2)");
    hBetaInvalid->GetXaxis()->SetBinLabel(2, "Valid Beta (beta > 0 and < 1)");
    hBetaInvalid->GetXaxis()->SetBinLabel(3, "Valid Beta (beta >= 1 and <= 1.2)");
    hBetaInvalid->GetYaxis()->SetTitle("Counts");

    hBetaInvalid->Draw("HIST");
    hBetaValidUsed->Draw("HIST SAME");
    hBetaValidExcluded->Draw("HIST SAME");

    TLegend* legendBeta = new TLegend(0.7, 0.7, 0.9, 0.9);
    legendBeta->AddEntry(hBetaInvalid, Form("Invalid Beta: %d", post_trigger_positive_beta_invalid), "f");
    legendBeta->AddEntry(hBetaValidUsed, Form("Valid Beta (Used): %d", post_trigger_positive_beta_p_valid - post_trigger_positive_beta_1_to_1p2), "f");
    legendBeta->AddEntry(hBetaValidExcluded, Form("Valid Beta (Excluded): %d", post_trigger_positive_beta_1_to_1p2), "f");
    legendBeta->Draw();

    canvasBetaBreakdown->Print("output/beta_breakdown_grouped.pdf");
    cout << "Grouped bar chart for beta breakdown saved to output/beta_breakdown_grouped.pdf" << endl;

    // Visualization 3: Grouped Bar Chart for Positive Charged Particles by PID and Chi2PID Validity (Forward Region)
    TCanvas* canvasChi2PIDForwardVis = new TCanvas("canvasChi2PIDForwardVis", "Breakdown of Positive Charged Particles by PID and Chi2PID Validity (Forward Region, Post-Trigger)", 1000, 600);

    int nPIDs = pids_to_visualize.size();
    TH1F* hChi2PIDValidForward = new TH1F("hChi2PIDValidForward", "Valid Chi2PID (Forward)", nPIDs * 2, 0, nPIDs * 2);
    TH1F* hChi2PIDInvalidForward = new TH1F("hChi2PIDInvalidForward", "Invalid Chi2PID (Forward)", nPIDs * 2, 0, nPIDs * 2);

    for (size_t i = 0; i < pids_to_visualize.size(); ++i) {
        int pid = pids_to_visualize[i];
        int bin = i * 2 + 1;

        hChi2PIDValidForward->SetBinContent(bin, post_trigger_positive_chi2pid_valid_counts_forward[pid]);
        hChi2PIDInvalidForward->SetBinContent(bin + 1, post_trigger_positive_chi2pid_invalid_counts_forward[pid]);

        hChi2PIDValidForward->GetXaxis()->SetBinLabel(bin, Form("PID %d", pid));
        hChi2PIDValidForward->GetXaxis()->SetBinLabel(bin + 1, Form("PID %d", pid));
    }

    hChi2PIDValidForward->SetFillColor(kGreen);
    hChi2PIDValidForward->SetLineColor(kBlack);
    hChi2PIDInvalidForward->SetFillColor(kRed);
    hChi2PIDInvalidForward->SetLineColor(kBlack);

    hChi2PIDValidForward->GetYaxis()->SetTitle("Counts");
    hChi2PIDValidForward->GetXaxis()->SetTitle("Particle ID (PID)");

    double maxHeightForward = hChi2PIDValidForward->GetMaximum();
    if (hChi2PIDInvalidForward->GetMaximum() > maxHeightForward) {
        maxHeightForward = hChi2PIDInvalidForward->GetMaximum();
    }
    hChi2PIDValidForward->SetMaximum(maxHeightForward * 1.3);

    hChi2PIDValidForward->Draw("HIST");
    hChi2PIDInvalidForward->Draw("HIST SAME");

    TLatex* latexForward = new TLatex();
    latexForward->SetTextSize(0.03);
    latexForward->SetTextAlign(22);

    for (size_t i = 0; i < pids_to_visualize.size(); ++i) {
        int pid = pids_to_visualize[i];
        int bin = i * 2 + 1;

        double validCount = post_trigger_positive_chi2pid_valid_counts_forward[pid];
        double xValid = hChi2PIDValidForward->GetBinCenter(bin);
        double yValid = validCount + maxHeightForward * 0.05;
        latexForward->SetTextColor(kGreen);
        latexForward->DrawLatex(xValid, yValid, Form("%d", (int)validCount));

        double invalidCount = post_trigger_positive_chi2pid_invalid_counts_forward[pid];
        double xInvalid = hChi2PIDInvalidForward->GetBinCenter(bin + 1);
        double yInvalid = invalidCount + maxHeightForward * 0.05;
        latexForward->SetTextColor(kRed);
        latexForward->DrawLatex(xInvalid, yInvalid, Form("%d", (int)invalidCount));
    }

    TLegend* legendChi2PIDForward = new TLegend(0.7, 0.7, 0.9, 0.9);
    int total_valid_forward = 0, total_invalid_forward = 0;
    for (const auto& pid : pids_to_visualize) {
        total_valid_forward += post_trigger_positive_chi2pid_valid_counts_forward[pid];
        total_invalid_forward += post_trigger_positive_chi2pid_invalid_counts_forward[pid];
    }
    legendChi2PIDForward->AddEntry(hChi2PIDValidForward, Form("Valid Chi2PID: %d", total_valid_forward), "f");
    legendChi2PIDForward->AddEntry(hChi2PIDInvalidForward, Form("Invalid Chi2PID: %d", total_invalid_forward), "f");
    legendChi2PIDForward->Draw();

    canvasChi2PIDForwardVis->Print("output/chi2pid_breakdown_by_pid_forward.pdf");
    cout << "Grouped bar chart for chi2pid breakdown by PID (Forward Region, Post-Trigger) saved to output/chi2pid_breakdown_by_pid_forward.pdf" << endl;

    // Visualization 4: Grouped Bar Chart for Positive Charged Particles by PID and Chi2PID Validity (Central Region)
    TCanvas* canvasChi2PIDCentralVis = new TCanvas("canvasChi2PIDCentralVis", "Breakdown of Positive Charged Particles by PID and Chi2PID Validity (Central Region, Post-Trigger)", 1000, 600);

    TH1F* hChi2PIDValidCentral = new TH1F("hChi2PIDValidCentral", "Valid Chi2PID (Central)", nPIDs * 2, 0, nPIDs * 2);
    TH1F* hChi2PIDInvalidCentral = new TH1F("hChi2PIDInvalidCentral", "Invalid Chi2PID (Central)", nPIDs * 2, 0, nPIDs * 2);

    for (size_t i = 0; i < pids_to_visualize.size(); ++i) {
        int pid = pids_to_visualize[i];
        int bin = i * 2 + 1;

        hChi2PIDValidCentral->SetBinContent(bin, post_trigger_positive_chi2pid_valid_counts_central[pid]);
        hChi2PIDInvalidCentral->SetBinContent(bin + 1, post_trigger_positive_chi2pid_invalid_counts_central[pid]);

        hChi2PIDValidCentral->GetXaxis()->SetBinLabel(bin, Form("PID %d", pid));
        hChi2PIDValidCentral->GetXaxis()->SetBinLabel(bin + 1, Form("PID %d", pid));
    }

    hChi2PIDValidCentral->SetFillColor(kGreen);
    hChi2PIDValidCentral->SetLineColor(kBlack);
    hChi2PIDInvalidCentral->SetFillColor(kRed);
    hChi2PIDInvalidCentral->SetLineColor(kBlack);

    hChi2PIDValidCentral->GetYaxis()->SetTitle("Counts");
    hChi2PIDValidCentral->GetXaxis()->SetTitle("Particle ID (PID)");

    double maxHeightCentral = hChi2PIDValidCentral->GetMaximum();
    if (hChi2PIDInvalidCentral->GetMaximum() > maxHeightCentral) {
        maxHeightCentral = hChi2PIDInvalidCentral->GetMaximum();
    }
    hChi2PIDValidCentral->SetMaximum(maxHeightCentral * 1.3);

    hChi2PIDValidCentral->Draw("HIST");
    hChi2PIDInvalidCentral->Draw("HIST SAME");

    TLatex* latexCentral = new TLatex();
    latexCentral->SetTextSize(0.03);
    latexCentral->SetTextAlign(22);

    for (size_t i = 0; i < pids_to_visualize.size(); ++i) {
        int pid = pids_to_visualize[i];
        int bin = i * 2 + 1;

        double validCount = post_trigger_positive_chi2pid_valid_counts_central[pid];
        double xValid = hChi2PIDValidCentral->GetBinCenter(bin);
        double yValid = validCount + maxHeightCentral * 0.05;
        latexCentral->SetTextColor(kGreen);
        latexCentral->DrawLatex(xValid, yValid, Form("%d", (int)validCount));

        double invalidCount = post_trigger_positive_chi2pid_invalid_counts_central[pid];
        double xInvalid = hChi2PIDInvalidCentral->GetBinCenter(bin + 1);
        double yInvalid = invalidCount + maxHeightCentral * 0.05;
        latexCentral->SetTextColor(kRed);
        latexCentral->DrawLatex(xInvalid, yInvalid, Form("%d", (int)invalidCount));
    }

    TLegend* legendChi2PIDCentral = new TLegend(0.7, 0.7, 0.9, 0.9);
    int total_valid_central = 0, total_invalid_central = 0;
    for (const auto& pid : pids_to_visualize) {
        total_valid_central += post_trigger_positive_chi2pid_valid_counts_central[pid];
        total_invalid_central += post_trigger_positive_chi2pid_invalid_counts_central[pid];
    }
    legendChi2PIDCentral->AddEntry(hChi2PIDValidCentral, Form("Valid Chi2PID: %d", total_valid_central), "f");
    legendChi2PIDCentral->AddEntry(hChi2PIDInvalidCentral, Form("Invalid Chi2PID: %d", total_invalid_central), "f");
    legendChi2PIDCentral->Draw();

    canvasChi2PIDCentralVis->Print("output/chi2pid_breakdown_by_pid_central.pdf");
    cout << "Grouped bar chart for chi2pid breakdown by PID (Central Region, Post-Trigger) saved to output/chi2pid_breakdown_by_pid_central.pdf" << endl;

    // Clean up
    delete hElectronMomentum;
    delete hMassPositiveAll;
    delete hChi2PIDPositivePionsForward;
    delete hChi2PIDPositivePionsCentral;
    delete canvas;
    delete canvasMassPositiveAll;
    delete canvasMassWithBkg;
    delete canvasChi2PIDForward;
    delete canvasChi2PIDCentral;
    delete gaussFit;
    delete gaussPlusBkg;
    delete chi2pidFitForward;
    delete chi2pidFitCentral;

    // Clean up visualization objects for charge breakdown
    delete hNeutral;
    delete hPositive;
    delete hNegative;
    delete hsCharge;
    delete legendCharge;
    delete canvasChargeBreakdown;

    // Clean up visualization objects for beta breakdown
    delete hBetaInvalid;
    delete hBetaValidUsed;
    delete hBetaValidExcluded;
    delete legendBeta;
    delete canvasBetaBreakdown;

    // Clean up visualization objects for chi2pid breakdown (Forward)
    delete hChi2PIDValidForward;
    delete hChi2PIDInvalidForward;
    delete legendChi2PIDForward;
    delete latexForward;
    delete canvasChi2PIDForwardVis;

    // Clean up visualization objects for chi2pid breakdown (Central)
    delete hChi2PIDValidCentral;
    delete hChi2PIDInvalidCentral;
    delete legendChi2PIDCentral;
    delete latexCentral;
    delete canvasChi2PIDCentralVis;

    return 0;
}