#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <filesystem>
#include <map>
#include "reader.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
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

    // Create output directory (relative to current directory)
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

    const int maxEvents = 10000000;
    int total_events_processed = 0;
    int event_count = 0;
    int trigger_electrons_events = 0;

    // Counters for particles in events with trigger electrons
    int total_positive_charge = 0;
    int total_negative_charge = 0;
    int positive_beta_invalid = 0;
    int positive_beta_1_to_1p2 = 0;
    int positive_beta_p_valid = 0;
    int pions_in_mass_range_updated = 0;
    int pions_pid_211_in_mass_range_updated = 0;
    int pions_pid_211_outside_mass_range = 0;
    int pions_in_mass_range_updated_forward = 0;
    int pions_pid_211_in_mass_range_updated_forward = 0;
    int pions_in_mass_range_updated_central = 0;
    int pions_pid_211_in_mass_range_updated_central = 0;
    int pid_211_with_chi2pid = 0;

    // Vectors to store data for mass and chi2pid analysis
    vector<float> positive_hadron_masses;
    vector<int> positive_hadron_pids;
    vector<float> positive_hadron_chi2pids;
    vector<bool> positive_hadron_is_forward;
    vector<bool> positive_hadron_is_central;
    vector<float> positive_hadron_momenta;

    // Histograms
    TH1F* hMassPositiveAll = new TH1F("hMassPositiveAll", "Mass of All Positive Hadrons;Mass (GeV/c^2);Counts", 100, 0, 1.5);
    TH1F* hChi2PIDPositivePionsAll = new TH1F("hChi2PIDPositivePionsAll", "Chi2PID for Positive Pions (PID 211, All Regions);Chi2PID;Counts", 100, -5, 5);
    TH1F* hChi2PIDPositivePionsForward = new TH1F("hChi2PIDPositivePionsForward", "Chi2PID for Positive Pions (PID 211, Forward);Chi2PID;Counts", 100, -5, 5);
    TH1F* hChi2PIDPositivePionsCentral = new TH1F("hChi2PIDPositivePionsCentral", "Chi2PID for Positive Pions (PID 211, Central);Chi2PID;Counts", 100, -5, 5);
    TH2F* hChi2PIDvsPAll = new TH2F("hChi2PIDvsPAll", "Chi2PID vs. Momentum for Positive Pions (All Regions);Momentum p (GeV);Chi2PID", 50, 0, 10, 50, -5, 5);

    // Define momentum bins from 1 to 8 GeV with 0.2 GeV width
    const double p_min = 1.0;
    const double p_max = 8.0;
    const double p_bin_width = 0.5;
    const int n_p_bins = static_cast<int>((p_max - p_min) / p_bin_width);
    vector<TH1F*> hChi2PIDvsPBins(n_p_bins);
    vector<double> p_bin_centers(n_p_bins);
    vector<double> chi2pid_means(n_p_bins);
    vector<double> chi2pid_sigmas(n_p_bins);

    // Initialize histograms for each momentum bin
    for (int i = 0; i < n_p_bins; ++i) {
        double p_low = p_min + i * p_bin_width;
        double p_high = p_low + p_bin_width;
        p_bin_centers[i] = (p_low + p_high) / 2.0;
        string hist_name = "hChi2PID_p_" + to_string(p_low) + "_to_" + to_string(p_high);
        string hist_title = "Chi2PID for Positive Pions (p = " + to_string(p_low) + " to " + to_string(p_high) + " GeV);Chi2PID;Counts";
        hChi2PIDvsPBins[i] = new TH1F(hist_name.c_str(), hist_title.c_str(), 100, -5, 5);
    }

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
            hipo::bank SCIN(dict.getSchema("REC::Scintillator"));
            hipo::bank CALO(dict.getSchema("REC::Calorimeter"));
            hipo::bank CHER(dict.getSchema("REC::Cherenkov"));
            hipo::event event;

            while (reader.next() && event_count < maxEvents) {
                reader.read(event);
                event.getStructure(PART);
                event.getStructure(SCIN);
                event.getStructure(CALO);
                event.getStructure(CHER);
                if (PART.getRows() == 0) continue;

                total_events_processed++;
                event_count++;

                bool hasTriggerElectron = false;

                // First pass: Check for a trigger electron in the event
                for (int i = 0; i < PART.getRows(); ++i) {
                    int pid = PART.getInt("pid", i);
                    int status = PART.getShort("status", i);
                    if (status == 0) continue;

                    if (pid == 11) { // Electron
                        bool is_trigger = (status < 0);
                        if (is_trigger) {
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

                // Second pass: If the event has a trigger electron, process particles for analysis
                if (hasTriggerElectron) {
                    trigger_electrons_events++;

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
                            total_positive_charge++;
                            
                            if (beta <= 0 || beta > 1.2) {
                                positive_beta_invalid++;
                            }
                            if (beta > 0 && beta <= 1.2 && p > 0) {
                                positive_beta_p_valid++;
                            }
                            if (beta >= 1 && beta <= 1.2) {
                                positive_beta_1_to_1p2++;
                            }

                            if (pid == 211 && chi2pid != 9999.0) {
                                pid_211_with_chi2pid++;
                                hChi2PIDvsPAll->Fill(p, chi2pid);
                                hChi2PIDPositivePionsAll->Fill(chi2pid);
                                if (isForward) {
                                    hChi2PIDPositivePionsForward->Fill(chi2pid);
                                }
                                if (isCentral) {
                                    hChi2PIDPositivePionsCentral->Fill(chi2pid);
                                }

                                // Fill the momentum-binned chi2pid histograms
                                if (p >= p_min && p < p_max) {
                                    int bin_index = static_cast<int>((p - p_min) / p_bin_width);
                                    if (bin_index >= 0 && bin_index < n_p_bins) {
                                        hChi2PIDvsPBins[bin_index]->Fill(chi2pid);
                                    }
                                }
                            }

                            if (beta > 0 && beta < 1) {
                                if (pid == -11 || pid == 45) continue;

                                float mass = (p / beta) * std::sqrt(1 - beta * beta);
                                hMassPositiveAll->Fill(mass);
                                positive_hadron_masses.push_back(mass);
                                positive_hadron_pids.push_back(pid);
                                positive_hadron_chi2pids.push_back(chi2pid);
                                positive_hadron_is_forward.push_back(isForward);
                                positive_hadron_is_central.push_back(isCentral);
                                positive_hadron_momenta.push_back(p);
                                if (positive_hadron_pids.size() <= 10) {
                                    cout << "Stored particle: PID = " << pid << ", Mass = " << mass << ", Momentum = " << p << endl;
                                }
                            }
                        } else if (charge < 0) {
                            total_negative_charge++;
                        }
                    }
                }

                if (event_count % 1000 == 0) {
                    cout << "Processed " << event_count << " events" << endl;
                }
            }
        }
    }

    // Fit with Gaussian + quadratic background for mass histogram
    TCanvas* canvasMassWithBkg = new TCanvas("canvasMassWithBkg", "Mass of Positive Hadrons with Background Fit", 800, 600);
    hMassPositiveAll->Draw("HIST");

    TF1* gaussPlusBkg = new TF1("gaussPlusBkg", "gaus(0) + pol2(3)", 0.0, 0.3);
    gaussPlusBkg->SetParameters(3000, 0.14, 0.02, 1000, -2000, 2000);
    gaussPlusBkg->SetParNames("Amplitude", "Mean", "Sigma", "Bkg_Const", "Bkg_Slope", "Bkg_Quad");
    hMassPositiveAll->Fit(gaussPlusBkg, "R");

    double pion_mass_mean_bkg = gaussPlusBkg->GetParameter(1);
    double pion_mass_sigma_bkg = gaussPlusBkg->GetParameter(2);
    double amplitude_bkg = gaussPlusBkg->GetParameter(0);
    double bkg_const = gaussPlusBkg->GetParameter(3);
    double bkg_slope = gaussPlusBkg->GetParameter(4);
    double bkg_quad = gaussPlusBkg->GetParameter(5);
    double pion_mass_low_bkg = pion_mass_mean_bkg - 3 * pion_mass_sigma_bkg;
    double pion_mass_high_bkg = pion_mass_mean_bkg + 3 * pion_mass_sigma_bkg;
    cout << "Pion peak fit (with background) - Mean: " << pion_mass_mean_bkg << " GeV/c^2, Sigma: " << pion_mass_sigma_bkg << " GeV/c^2" << endl;
    cout << "Background parameters - Const: " << bkg_const << ", Slope: " << bkg_slope << ", Quad: " << bkg_quad << endl;
    cout << "Pion mass selection range (with background, 3 sigma): [" << pion_mass_low_bkg << ", " << pion_mass_high_bkg << "] GeV/c^2" << endl;

    gaussPlusBkg->SetLineColor(kRed);
    gaussPlusBkg->SetLineWidth(2);
    gaussPlusBkg->Draw("SAME");

    TF1* pionComponent = new TF1("pionComponent", "gaus", 0.0, 0.3);
    pionComponent->SetParameters(amplitude_bkg, pion_mass_mean_bkg, pion_mass_sigma_bkg);
    pionComponent->SetLineColor(kGreen);
    pionComponent->SetLineStyle(2);
    pionComponent->Draw("SAME");

    TF1* bkgComponent = new TF1("bkgComponent", "pol2", 0.0, 0.3);
    bkgComponent->SetParameters(bkg_const, bkg_slope, bkg_quad);
    bkgComponent->SetLineColor(kBlack);
    bkgComponent->SetLineStyle(2);
    bkgComponent->Draw("SAME");

    canvasMassWithBkg->Update();
    canvasMassWithBkg->Print("output/mass_positive_all_with_fit_and_bkg.pdf");

    // Count particles in the mass range
    for (size_t i = 0; i < positive_hadron_masses.size(); ++i) {
        float mass = positive_hadron_masses[i];
        int pid = positive_hadron_pids[i];
        bool isForward = positive_hadron_is_forward[i];
        bool isCentral = positive_hadron_is_central[i];

        if (mass >= pion_mass_low_bkg && mass <= pion_mass_high_bkg) {
            pions_in_mass_range_updated++;
            if (pid == 211) {
                pions_pid_211_in_mass_range_updated++;
            }
            if (isForward) {
                pions_in_mass_range_updated_forward++;
                if (pid == 211) {
                    pions_pid_211_in_mass_range_updated_forward++;
                }
            }
            if (isCentral) {
                pions_in_mass_range_updated_central++;
                if (pid == 211) {
                    pions_pid_211_in_mass_range_updated_central++;
                }
            }
        }

        if (!(mass >= pion_mass_low_bkg && mass <= pion_mass_high_bkg)) {
            if (pid == 211) {
                pions_pid_211_outside_mass_range++;
            }
        }
    }

    // Contamination estimates for Gaussian + quadratic background fit
    double bkg_integral_bkg = (bkg_const * (pion_mass_high_bkg - pion_mass_low_bkg)) + 
                              (0.5 * bkg_slope * (pion_mass_high_bkg * pion_mass_high_bkg - pion_mass_low_bkg * pion_mass_low_bkg)) + 
                              (bkg_quad / 3.0) * (pion_mass_high_bkg * pion_mass_high_bkg * pion_mass_high_bkg - 
                                                  pion_mass_low_bkg * pion_mass_low_bkg * pion_mass_low_bkg);
    double gauss_integral_bkg = amplitude_bkg * sqrt(2 * TMath::Pi()) * pion_mass_sigma_bkg * 
                                (TMath::Erf((pion_mass_high_bkg - pion_mass_mean_bkg) / (sqrt(2) * pion_mass_sigma_bkg)) - 
                                 TMath::Erf((pion_mass_low_bkg - pion_mass_mean_bkg) / (sqrt(2) * pion_mass_sigma_bkg))) / 2;
    double total_integral_bkg = bkg_integral_bkg + gauss_integral_bkg;
    double contamination_fraction_bkg = (total_integral_bkg > 0) ? (bkg_integral_bkg / total_integral_bkg) * 100.0 : 0.0;

    double bkg_integral_bkg_forward = bkg_integral_bkg * (pions_in_mass_range_updated_forward / (double)pions_in_mass_range_updated);
    double bkg_integral_bkg_central = bkg_integral_bkg * (pions_in_mass_range_updated_central / (double)pions_in_mass_range_updated);
    double total_integral_bkg_forward = total_integral_bkg * (pions_in_mass_range_updated_forward / (double)pions_in_mass_range_updated);
    double total_integral_bkg_central = total_integral_bkg * (pions_in_mass_range_updated_central / (double)pions_in_mass_range_updated);
    double contamination_fraction_bkg_forward = (total_integral_bkg_forward > 0) ? (bkg_integral_bkg_forward / total_integral_bkg_forward) * 100.0 : 0.0;
    double contamination_fraction_bkg_central = (total_integral_bkg_central > 0) ? (bkg_integral_bkg_central / total_integral_bkg_central) * 100.0 : 0.0;

    // Fit and plot chi2pid histogram for all regions
    TCanvas* canvasChi2PIDAll = new TCanvas("canvasChi2PIDAll", "Chi2PID for Positive Pions (PID 211, All Regions)", 800, 600);
    hChi2PIDPositivePionsAll->Draw("HIST");

    TF1* chi2pidFitAll = new TF1("chi2pidFitAll", "gaus", -2, 2);
    chi2pidFitAll->SetParameters(1000, 0, 1);
    hChi2PIDPositivePionsAll->Fit(chi2pidFitAll, "R");

    double chi2pid_mean_all = chi2pidFitAll->GetParameter(1);
    double chi2pid_sigma_all = chi2pidFitAll->GetParameter(2);
    double chi2pid_low_cut_all = chi2pid_mean_all - 3 * chi2pid_sigma_all;
    double chi2pid_high_cut_all = chi2pid_mean_all + 3 * chi2pid_sigma_all;

    chi2pidFitAll->SetLineColor(kRed);
    chi2pidFitAll->Draw("SAME");
    canvasChi2PIDAll->Print("output/chi2pid_positive_pions_all.pdf");

    // Fit chi2pid histograms for forward and central regions
    TCanvas* canvasChi2PIDForward = new TCanvas("canvasChi2PIDForward", "Chi2PID for Positive Pions (PID 211, Forward)", 800, 600);
    hChi2PIDPositivePionsForward->Draw("HIST");

    TF1* chi2pidFitForward = new TF1("chi2pidFitForward", "gaus", -2, 2);
    chi2pidFitForward->SetParameters(1000, 0, 1);
    hChi2PIDPositivePionsForward->Fit(chi2pidFitForward, "R");

    double chi2pid_mean_forward = chi2pidFitForward->GetParameter(1);
    double chi2pid_sigma_forward = chi2pidFitForward->GetParameter(2);
    double chi2pid_low_cut_forward = chi2pid_mean_forward - 3 * chi2pid_sigma_forward;
    double chi2pid_high_cut_forward = chi2pid_mean_forward + 3 * chi2pid_sigma_forward;

    chi2pidFitForward->SetLineColor(kRed);
    chi2pidFitForward->Draw("SAME");
    canvasChi2PIDForward->Print("output/chi2pid_positive_pions_forward.pdf");

    TCanvas* canvasChi2PIDCentral = new TCanvas("canvasChi2PIDCentral", "Chi2PID for Positive Pions (PID 211, Central)", 800, 600);
    hChi2PIDPositivePionsCentral->Draw("HIST");

    TF1* chi2pidFitCentral = new TF1("chi2pidFitCentral", "gaus", -2, 2);
    chi2pidFitCentral->SetParameters(1000, 0, 1);
    hChi2PIDPositivePionsCentral->Fit(chi2pidFitCentral, "R");

    double chi2pid_mean_central = chi2pidFitCentral->GetParameter(1);
    double chi2pid_sigma_central = chi2pidFitCentral->GetParameter(2);
    double chi2pid_low_cut_central = chi2pid_mean_central - 3 * chi2pid_sigma_central;
    double chi2pid_high_cut_central = chi2pid_mean_central + 3 * chi2pid_sigma_central;

    chi2pidFitCentral->SetLineColor(kRed);
    chi2pidFitCentral->Draw("SAME");
    canvasChi2PIDCentral->Print("output/chi2pid_positive_pions_central.pdf");

    // Plot chi2pid vs. p histogram (all regions)
    TCanvas* canvasChi2PIDvsPAll = new TCanvas("canvasChi2PIDvsPAll", "Chi2PID vs. Momentum (All Regions)", 800, 600);
    hChi2PIDvsPAll->Draw("COLZ");
    gPad->SetLogz();
    canvasChi2PIDvsPAll->Print("output/chi2pid_vs_p_all.pdf");

    // Create chi2pid histograms for each momentum bin, fit with Gaussian, and save to PDF
    TCanvas* canvasChi2PIDBins = new TCanvas("canvasChi2PIDBins", "Chi2PID in Momentum Bins", 800, 600);
    canvasChi2PIDBins->Print("output/chi2pid_vs_p_bins.pdf["); // Open PDF file

    for (int i = 0; i < n_p_bins; ++i) {
        if (hChi2PIDvsPBins[i]->GetEntries() == 0) {
            chi2pid_means[i] = 0.0;
            chi2pid_sigmas[i] = 0.0;
            continue;
        }

        canvasChi2PIDBins->cd();
        hChi2PIDvsPBins[i]->Draw("HIST");

        // Fit with Gaussian
        TF1* chi2pidFitBin = new TF1(("chi2pidFitBin_" + to_string(i)).c_str(), "gaus", -1, 1);
        chi2pidFitBin->SetParameters(hChi2PIDvsPBins[i]->GetMaximum(), 0, 1);
        hChi2PIDvsPBins[i]->Fit(chi2pidFitBin, "R");

        chi2pid_means[i] = chi2pidFitBin->GetParameter(1);
        chi2pid_sigmas[i] = chi2pidFitBin->GetParameter(2);

        chi2pidFitBin->SetLineColor(kRed);
        chi2pidFitBin->Draw("SAME");

        // Add text with fit parameters
        TLatex* latex = new TLatex();
        latex->SetNDC();
        latex->SetTextSize(0.04);
        latex->DrawLatex(0.6, 0.85, Form("Mean: %.3f", chi2pid_means[i]));
        latex->DrawLatex(0.6, 0.80, Form("Sigma: %.3f", chi2pid_sigmas[i]));

        canvasChi2PIDBins->Update();
        canvasChi2PIDBins->Print("output/chi2pid_vs_p_bins.pdf");

        delete chi2pidFitBin;
        delete latex;
    }

    canvasChi2PIDBins->Print("output/chi2pid_vs_p_bins.pdf]"); // Close PDF file

    // Save chi2pid fit parameters to a new CSV file
    ofstream csvChi2PIDBins("output/chi2pid_vs_p_bins.csv");
    if (!csvChi2PIDBins.is_open()) {
        cerr << "Error: Could not open output/chi2pid_vs_p_bins.csv for writing!" << endl;
        return 1;
    }

    csvChi2PIDBins << "Momentum_Low (GeV),Momentum_High (GeV),Momentum_Center (GeV),Chi2PID_Mean,Chi2PID_Sigma\n";
    for (int i = 0; i < n_p_bins; ++i) {
        double p_low = p_min + i * p_bin_width;
        double p_high = p_low + p_bin_width;
        csvChi2PIDBins << p_low << "," << p_high << "," << p_bin_centers[i] << ","
                       << chi2pid_means[i] << "," << chi2pid_sigmas[i] << "\n";
    }
    csvChi2PIDBins.close();
    cout << "Chi2PID fit parameters saved to output/chi2pid_vs_p_bins.csv" << endl;

    // Apply chi2pid cut
    int pions_pid_211_in_mass_range_after_chi2pid = 0;
    int pions_pid_211_in_mass_range_after_chi2pid_forward = 0;
    int pions_pid_211_in_mass_range_after_chi2pid_central = 0;

    for (size_t i = 0; i < positive_hadron_masses.size(); ++i) {
        float mass = positive_hadron_masses[i];
        int pid = positive_hadron_pids[i];
        float chi2pid = positive_hadron_chi2pids[i];
        bool isForward = positive_hadron_is_forward[i];
        bool isCentral = positive_hadron_is_central[i];

        if (mass >= pion_mass_low_bkg && mass <= pion_mass_high_bkg && pid == 211 && chi2pid != 9999.0) {
            if (isForward && chi2pid >= chi2pid_low_cut_forward && chi2pid <= chi2pid_high_cut_forward) {
                pions_pid_211_in_mass_range_after_chi2pid_forward++;
                pions_pid_211_in_mass_range_after_chi2pid++;
            }
            if (isCentral && chi2pid >= chi2pid_low_cut_central && chi2pid <= chi2pid_high_cut_central) {
                pions_pid_211_in_mass_range_after_chi2pid_central++;
                if (!isForward || !(chi2pid >= chi2pid_low_cut_forward && chi2pid <= chi2pid_high_cut_forward)) {
                    pions_pid_211_in_mass_range_after_chi2pid++;
                }
            }
        }
    }

    // Recalculate contamination after chi2pid cut
    double total_after_chi2pid = pions_in_mass_range_updated - (pions_pid_211_in_mass_range_updated - pions_pid_211_in_mass_range_after_chi2pid);
    double contamination_fraction_after_chi2pid = (total_after_chi2pid > 0) ? (bkg_integral_bkg / total_after_chi2pid) * 100.0 : 0.0;

    double total_after_chi2pid_forward = pions_in_mass_range_updated_forward - (pions_pid_211_in_mass_range_updated_forward - pions_pid_211_in_mass_range_after_chi2pid_forward);
    double total_after_chi2pid_central = pions_in_mass_range_updated_central - (pions_pid_211_in_mass_range_updated_central - pions_pid_211_in_mass_range_after_chi2pid_central);
    double contamination_fraction_after_chi2pid_forward = (total_after_chi2pid_forward > 0) ? (bkg_integral_bkg_forward / total_after_chi2pid_forward) * 100.0 : 0.0;
    double contamination_fraction_after_chi2pid_central = (total_after_chi2pid_central > 0) ? (bkg_integral_bkg_central / total_after_chi2pid_central) * 100.0 : 0.0;

    timer.Stop();
    double execution_time = timer.RealTime();

    int mass_histogram_entries = hMassPositiveAll->GetEntries();

    // Write to stats CSV
    ofstream csvFile("output/stats.csv");
    if (!csvFile.is_open()) {
        cerr << "Error: Could not open output/stats.csv for writing!" << endl;
        return 1;
    }

    csvFile << "Statistic,Value\n";
    csvFile << "Total events processed," << total_events_processed << "\n";
    csvFile << "Events with trigger electrons," << trigger_electrons_events << "\n";
    csvFile << "Total positive charged particles," << total_positive_charge << "\n";
    csvFile << "Total negative charged particles," << total_negative_charge << "\n";
    csvFile << "Positive particles with invalid beta (beta <= 0 or beta > 1.2)," << positive_beta_invalid << "\n";
    csvFile << "Positive particles with valid beta and p (plotted in beta vs. p)," << positive_beta_p_valid << "\n";
    csvFile << "Positive particles with 1 <= beta <= 1.2 (excluded from mass histogram)," << positive_beta_1_to_1p2 << "\n";
    csvFile << "Entries in mass histogram (hMassPositiveAll)," << mass_histogram_entries << "\n";
    csvFile << "Pion peak mean (with background) (GeV/c^2)," << pion_mass_mean_bkg << "\n";
    csvFile << "Pion peak sigma (with background) (GeV/c^2)," << pion_mass_sigma_bkg << "\n";
    csvFile << "Pion mass selection range low (with background, 3 sigma) (GeV/c^2)," << pion_mass_low_bkg << "\n";
    csvFile << "Pion mass selection range high (with background, 3 sigma) (GeV/c^2)," << pion_mass_high_bkg << "\n";
    csvFile << "Number of particles in pion mass range (with background)," << pions_in_mass_range_updated << "\n";
    csvFile << "Number of particles with pid == 211 in mass range," << pions_pid_211_in_mass_range_updated << "\n";
    csvFile << "Estimated contamination fraction in pion sample (with background) (%)," << contamination_fraction_bkg << "\n";
    csvFile << "Number of particles with pid == 211 outside the mass range," << pions_pid_211_outside_mass_range << "\n";
    csvFile << "Number of particles in pion mass range (Forward)," << pions_in_mass_range_updated_forward << "\n";
    csvFile << "Number of particles with pid == 211 in mass range (Forward)," << pions_pid_211_in_mass_range_updated_forward << "\n";
    csvFile << "Estimated contamination fraction in pion sample (Forward) (%)," << contamination_fraction_bkg_forward << "\n";
    csvFile << "Number of particles in pion mass range (Central)," << pions_in_mass_range_updated_central << "\n";
    csvFile << "Number of particles with pid == 211 in mass range (Central)," << pions_pid_211_in_mass_range_updated_central << "\n";
    csvFile << "Estimated contamination fraction in pion sample (Central) (%)," << contamination_fraction_bkg_central << "\n";
    csvFile << "Number of particles with pid == 211 and valid chi2pid," << pid_211_with_chi2pid << "\n";
    csvFile << "Chi2PID mean for positive pions (All Regions)," << chi2pid_mean_all << "\n";
    csvFile << "Chi2PID sigma for positive pions (All Regions)," << chi2pid_sigma_all << "\n";
    csvFile << "Chi2PID cut range low (All Regions, 3 sigma)," << chi2pid_low_cut_all << "\n";
    csvFile << "Chi2PID cut range high (All Regions, 3 sigma)," << chi2pid_high_cut_all << "\n";
    csvFile << "Chi2PID mean for positive pions (Forward)," << chi2pid_mean_forward << "\n";
    csvFile << "Chi2PID sigma for positive pions (Forward)," << chi2pid_sigma_forward << "\n";
    csvFile << "Chi2PID cut range low (Forward, 3 sigma)," << chi2pid_low_cut_forward << "\n";
    csvFile << "Chi2PID cut range high (Forward, 3 sigma)," << chi2pid_high_cut_forward << "\n";
    csvFile << "Chi2PID mean for positive pions (Central)," << chi2pid_mean_central << "\n";
    csvFile << "Chi2PID sigma for positive pions (Central)," << chi2pid_sigma_central << "\n";
    csvFile << "Chi2PID cut range low (Central, 3 sigma)," << chi2pid_low_cut_central << "\n";
    csvFile << "Chi2PID cut range high (Central, 3 sigma)," << chi2pid_high_cut_central << "\n";
    csvFile << "Number of particles with pid == 211 in mass range after chi2pid cut (overall)," << pions_pid_211_in_mass_range_after_chi2pid << "\n";
    csvFile << "Estimated contamination fraction after chi2pid cut (overall) (%)," << contamination_fraction_after_chi2pid << "\n";
    csvFile << "Number of particles with pid == 211 in mass range after chi2pid cut (Forward)," << pions_pid_211_in_mass_range_after_chi2pid_forward << "\n";
    csvFile << "Estimated contamination fraction after chi2pid cut (Forward) (%)," << contamination_fraction_after_chi2pid_forward << "\n";
    csvFile << "Number of particles with pid == 211 in mass range after chi2pid cut (Central)," << pions_pid_211_in_mass_range_after_chi2pid_central << "\n";
    csvFile << "Estimated contamination fraction after chi2pid cut (Central) (%)," << contamination_fraction_after_chi2pid_central << "\n";
    csvFile << "Execution time (seconds)," << execution_time << "\n";

    csvFile.close();
    cout << "Statistics saved to output/stats.csv" << endl;

    // Save histograms to a ROOT file
    TFile* outputFile = new TFile("output/histograms.root", "RECREATE");
    hMassPositiveAll->Write();
    hChi2PIDPositivePionsForward->Write();
    hChi2PIDPositivePionsCentral->Write();
    hChi2PIDvsPAll->Write();
    hChi2PIDPositivePionsAll->Write();
    for (int i = 0; i < n_p_bins; ++i) {
        hChi2PIDvsPBins[i]->Write();
    }
    outputFile->Close();
    delete outputFile;

    // Clean up
    delete hMassPositiveAll;
    delete hChi2PIDPositivePionsForward;
    delete hChi2PIDPositivePionsCentral;
    delete hChi2PIDvsPAll;
    delete hChi2PIDPositivePionsAll;
    for (int i = 0; i < n_p_bins; ++i) {
        delete hChi2PIDvsPBins[i];
    }
    delete canvasMassWithBkg;
    delete canvasChi2PIDForward;
    delete canvasChi2PIDCentral;
    delete canvasChi2PIDvsPAll;
    delete canvasChi2PIDAll;
    delete canvasChi2PIDBins;

    return 0;
}