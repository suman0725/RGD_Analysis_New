#include <cstdlib>
#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <fstream>
#include <string>
#include <filesystem>
#include "reader.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include <cstddef>
#include <TStopwatch.h>
#include <TLegend.h>
#include <limits>
#include "TStyle.h"
#include <set>
#include "ParticleData.h"
#include "PidFunctions.h"
#include "Validation.h"

using namespace std;
namespace fs = std::filesystem;

int main() {
    TStopwatch timer;
    timer.Start();

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

    loadCCDBParams();

    const int maxEvents = 100000;
    int event_count = 0;
    int total_events_processed = 0;
    int directory_count = 0;
    int hipo_file_count = 0;
    set<int> unique_events;
    ParticleStats stats;
    vector<Mismatch> mismatches;
    vector<float> chi2pidList;
    vector<float> chi2pidCustomList;

    // Initialize histograms (unchanged)
    map<int, TH1F*> hMomentumRec;
    map<int, TH1F*> hMomentumCustom;
    vector<int> particle_types = {11, -11, 211, -211, 321, -321, 2212, -2212, 45};
    for (int pid : particle_types) {
        string name_rec = "hMomentumRec_" + to_string(pid);
        string name_custom = "hMomentumCustom_" + to_string(pid);
        string title_rec = "Momentum (REC::PID " + to_string(pid) + ");Momentum (GeV);Counts";
        string title_custom = "Momentum (Assigned PID " + to_string(pid) + ");Momentum (GeV);Counts";
        hMomentumRec[pid] = new TH1F(name_rec.c_str(), title_rec.c_str(), 100, 0, 10);
        hMomentumCustom[pid] = new TH1F(name_custom.c_str(), title_custom.c_str(), 100, 0, 10);
        hMomentumRec[pid]->SetLineColor(kBlue);
        hMomentumRec[pid]->SetFillColor(kBlue);
        hMomentumRec[pid]->SetFillStyle(3004);
        hMomentumCustom[pid]->SetLineColor(kRed);
        hMomentumCustom[pid]->SetFillColor(kRed);
        hMomentumCustom[pid]->SetFillStyle(3005);
    }

    // Initialize chi2pid histograms (unchanged)
    map<int, TH1F*> hChi2pidRec;
    map<int, TH1F*> hChi2pidCustom;
    for (int pid : particle_types) {
        string name_rec = "hChi2pidRec_" + to_string(pid);
        string name_custom = "hChi2pidCustom_" + to_string(pid);
        string title_rec = "Chi2pid (REC::PID " + to_string(pid) + ");Chi2pid;Counts";
        string title_custom = "Chi2pid_custom (Assigned PID " + to_string(pid) + ");Chi2pid;Counts";
        hChi2pidRec[pid] = new TH1F(name_rec.c_str(), title_rec.c_str(), 100, -15, 15);
        hChi2pidCustom[pid] = new TH1F(name_custom.c_str(), title_custom.c_str(), 100, -15, 15);
        hChi2pidRec[pid]->SetLineColor(kBlue);
        hChi2pidRec[pid]->SetFillColor(kBlue);
        hChi2pidRec[pid]->SetFillStyle(3004);
        hChi2pidCustom[pid]->SetLineColor(kRed);
        hChi2pidCustom[pid]->SetFillColor(kRed);
        hChi2pidCustom[pid]->SetFillStyle(3005);
    }

    // Initialize Delta T histograms (unchanged)
    map<int, TH1F*> hDtRec;
    map<int, TH1F*> hDtCustom;
    for (int pid : particle_types) {
        string name_rec = "hDtRec_" + to_string(pid);
        string name_custom = "hDtCustom_" + to_string(pid);
        string title_rec = "Delta T (REC::PID " + to_string(pid) + ");Delta T (ns);Counts";
        string title_custom = "Delta T (Assigned PID " + to_string(pid) + ");Delta T (ns);Counts";
        hDtRec[pid] = new TH1F(name_rec.c_str(), title_rec.c_str(), 100, -5, 5);
        hDtCustom[pid] = new TH1F(name_custom.c_str(), title_custom.c_str(), 100, -5, 5);
        hDtRec[pid]->SetLineColor(kBlue);
        hDtRec[pid]->SetFillColor(kBlue);
        hDtRec[pid]->SetFillStyle(3004);
        hDtCustom[pid]->SetLineColor(kRed);
        hDtCustom[pid]->SetFillColor(kRed);
        hDtCustom[pid]->SetFillStyle(3005);
    }

    // Initialize 2D histograms (unchanged)
    map<int, TH2F*> hDtVsPRec;
    map<int, TH2F*> hDtVsPCustom;
    for (int pid : particle_types) {
        string name_rec = "hDtVsPRec_" + to_string(pid);
        string name_custom = "hDtVsPCustom_" + to_string(pid);
        string title_rec = "Delta T vs Momentum (REC::PID " + to_string(pid) + ");Momentum (GeV);Delta T (ns)";
        string title_custom = "Delta T vs Momentum (Assigned PID " + to_string(pid) + ");Momentum (GeV);Delta T (ns)";
        hDtVsPRec[pid] = new TH2F(name_rec.c_str(), title_rec.c_str(), 100, 0, 10, 100, -5, 5);
        hDtVsPCustom[pid] = new TH2F(name_custom.c_str(), title_custom.c_str(), 100, 0, 10, 100, -5, 5);
        hDtVsPRec[pid]->SetMarkerColor(kBlue);
        hDtVsPCustom[pid]->SetMarkerColor(kRed);
    }

    map<int, TH2F*> hChi2pidVsPRec;
    map<int, TH2F*> hChi2pidVsPCustom;
    for (int pid : particle_types) {
        string name_rec = "hChi2pidVsPRec_" + to_string(pid);
        string name_custom = "hChi2pidVsPCustom_" + to_string(pid);
        string title_rec = "Chi2pid vs Momentum (REC::PID " + to_string(pid) + ");Momentum (GeV);Chi2pid";
        string title_custom = "Chi2pid vs Momentum (Assigned PID " + to_string(pid) + ");Momentum (GeV);Chi2pid";
        hChi2pidVsPRec[pid] = new TH2F(name_rec.c_str(), title_rec.c_str(), 100, 0, 10, 100, -15, 15);
        hChi2pidVsPCustom[pid] = new TH2F(name_custom.c_str(), title_custom.c_str(), 100, 0, 10, 100, -15, 15);
        hChi2pidVsPRec[pid]->SetMarkerColor(kBlue);
        hChi2pidVsPCustom[pid]->SetMarkerColor(kRed);
    }

    // Main processing loop (unchanged)
    for (const auto& dir : directories) {
        directory_count++;
        vector<string> hipoFiles;
        for (const auto& entry : fs::directory_iterator(dir)) {
            if (entry.is_regular_file() && entry.path().extension() == ".hipo") {
                hipoFiles.push_back(entry.path().string());
            }
        }
        if (hipoFiles.empty()) {
            cout << "No .hipo files found in directory: " << dir << endl;
            continue;
        }

        for (const auto& file : hipoFiles) {
            hipo_file_count++;
            cout << "Processing file: " << file << endl;
            hipo::reader reader;
            reader.open(file.c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);

            while (reader.next() && event_count < maxEvents) {
                event_count++;
                unique_events.insert(event_count);
                total_events_processed++;

                hipo::event event;
                reader.read(event);
                hipo::bank PART(factory.getSchema("REC::Particle"));
                hipo::bank EVENT(factory.getSchema("REC::Event"));
                hipo::bank SCIN(factory.getSchema("REC::Scintillator"));
                hipo::bank CHER(factory.getSchema("REC::Cherenkov"));
                hipo::bank CALO(factory.getSchema("REC::Calorimeter"));
                event.getStructure(PART);
                event.getStructure(EVENT);
                event.getStructure(SCIN);
                event.getStructure(CHER);
                event.getStructure(CALO);

                set<tuple<int, int, int, float, float, float>> unique_particles_per_event;
                for (int i = 0; i < PART.getRows(); i++) {
                    int pid = PART.getInt("pid", i);
                    int charge = PART.getByte("charge", i);
                    int status = PART.getShort("status", i);
                    float px = PART.getFloat("px", i);
                    float py = PART.getFloat("py", i);
                    float pz = PART.getFloat("pz", i);
                    unique_particles_per_event.insert(make_tuple(pid, charge, status, px, py, pz));
                }
                if (unique_particles_per_event.size() != PART.getRows()) {
                    cout << "Event " << event_count << ": Found duplicates! Unique particles = "
                         << unique_particles_per_event.size() << ", Total rows = " << PART.getRows() << endl;
                }

                IndexMap cherMap = loadMapByIndex(CHER, "pindex");
                IndexMap caloMap = loadMapByIndex(CALO, "pindex");
                IndexMap scinMap = loadMapByIndex(SCIN, "pindex");

                float startTime = EVENT.getFloat("startTime", 0);

                map<tuple<int, int, int, float, float, float>, int> particle_indices;
                vector<int> unique_indices;
                for (int i = 0; i < PART.getRows(); i++) {
                    int pid = PART.getInt("pid", i);
                    int charge = PART.getByte("charge", i);
                    int status = PART.getShort("status", i);
                    float px = PART.getFloat("px", i);
                    float py = PART.getFloat("py", i);
                    float pz = PART.getFloat("pz", i);
                    auto key = make_tuple(pid, charge, status, px, py, pz);
                    if (particle_indices.find(key) == particle_indices.end()) {
                        particle_indices[key] = i;
                        unique_indices.push_back(i);
                    }
                }

                vector<ParticleData> particles(unique_indices.size());
                int valid_particles = 0;
                for (int idx : unique_indices) {
                    ParticleData pd = getParticleData(idx, PART, CHER, CALO, SCIN, cherMap, caloMap, scinMap);
                    stats.total_particles_before_filtering++;
                    if (pd.pid == 0) stats.total_unidentified_before_filtering++;
                    if (pd.status == 0 || pd.charge == 0) {
                        stats.total_filtered_out++;
                        if (pd.pid == 0) stats.total_filtered_out_pid_zero++;
                        continue;
                    }
                    particles[valid_particles++] = pd;
                }
                particles.resize(valid_particles);

                assignPids(particles, startTime, event_count);
                validateAndCompare(particles, event_count, stats, mismatches, chi2pidList, chi2pidCustomList,
                                   hMomentumRec, hMomentumCustom, hChi2pidRec, hChi2pidCustom,
                                   hDtRec, hDtCustom, hDtVsPRec, hDtVsPCustom,
                                   hChi2pidVsPRec, hChi2pidVsPCustom);
            }
            if (event_count >= maxEvents) break;
        }
        if (event_count >= maxEvents) break;
    }

    timer.Stop();

    cout << "****************************************************************************************" << endl;
    cout << "Total number of events: " << event_count << endl;
    cout << "Total unique events: " << unique_events.size() << endl;
    cout << "Total particles before filtering: " << stats.total_particles_before_filtering << endl;
    cout << "Total unidentified before filtering (pid == 0): " << stats.total_unidentified_before_filtering << endl;
    cout << "Total particles filtered out (status == 0 || charge == 0): " << stats.total_filtered_out << endl;
    cout << "Total particles filtered out with pid == 0: " << stats.total_filtered_out_pid_zero << endl;

    // Print PID Counts Table for each region
    cout << "\nPID Counts by Region:" << endl;
    const std::map<Region, std::string> region_names = {
        {Region::Forward, "Forward"},
        {Region::Central, "Central"},
        {Region::Band, "Band"},
        {Region::Unknown, "Unknown"}
    };

    for (const auto& region_pair : region_names) {
        Region region = region_pair.first;
        cout << "\nRegion: " << region_pair.second << endl;
        cout << "pid\tPid from REC\tCustom Pid" << endl;
        for (int pid : particle_types) {
            cout << pid << "\t" << stats.pid_counts_by_region[region][pid] << "\t" << stats.assigned_pid_counts_by_region[region][pid] << endl;
        }
        cout << "0\t" << stats.pid_counts_by_region[region][0] << "\t" << stats.assigned_pid_counts_by_region[region][0] << endl;

        int total_particles_with_pid_rec = 0;
        int total_particles_with_pid_custom = 0;
        for (int pid : particle_types) {
            total_particles_with_pid_rec += stats.pid_counts_by_region[region][pid];
            total_particles_with_pid_custom += stats.assigned_pid_counts_by_region[region][pid];
        }
        cout << "Total particle with pid\t" << total_particles_with_pid_rec << "\t" << total_particles_with_pid_custom << endl;
    }

    // Print overall PID Counts Table (unchanged)
    cout << "\nOverall PID Counts Table:" << endl;
    cout << "pid\tPid from REC\tCustom Pid" << endl;
    for (int pid : particle_types) {
        cout << pid << "\t" << stats.pid_counts[pid] << "\t" << stats.assigned_pid_counts[pid] << endl;
    }
    cout << "0\t" << stats.pid_counts[0] << "\t" << stats.assigned_pid_counts[0] << endl;

    int total_particles_with_pid_rec = 0;
    int total_particles_with_pid_custom = 0;
    for (const auto& pair : stats.pid_counts) {
        if (pair.first != 0) total_particles_with_pid_rec += pair.second;
    }
    for (const auto& pair : stats.assigned_pid_counts) {
        if (pair.first != 0) total_particles_with_pid_custom += pair.second;
    }
    cout << "Total particle with pid\t" << total_particles_with_pid_rec << "\t" << total_particles_with_pid_custom << endl;

    // Print chi2pid and chi2pid_custom Statistics (unchanged)
    cout << "\nchi2pid Statistics:" << endl;
    cout << "valid: chi2pid != 9999.0\t" << stats.valid_chi2pid << endl;
    cout << "invalid: chi2pid = 9999.0\t" << stats.total_invalid_chi2pid << endl;
    cout << "valid: chi2pid != 9999.0 && pid == 0\t" << 0 << endl;
    cout << "invalid: chi2pid == 9999.0 && pid == 0\t" << stats.pid_counts[0] << endl;
    cout << "valid: chi2pid != 9999.0 && pid != 0\t" << (stats.valid_chi2pid - 0) << endl;
    cout << "invalid: chi2pid == 9999.0 && pid != 0\t" << stats.total_invalid_chi2pid_nonzero_pid << endl;

    cout << "\nchi2pid_custom Statistics:" << endl;
    cout << "valid: chi2pid_custom != 99999.0\t" << stats.valid_chi2pid_custom << endl;
    cout << "invalid: chi2pid_custom == 99999.0\t" << stats.total_invalid_chi2pid_custom << endl;
    cout << "valid: chi2pid_custom != 99999.0 && Assigned PID == 0\t" << stats.valid_chi2pid_custom_pid_zero << endl;
    cout << "invalid: chi2pid_custom == 99999.0 && Assigned PID == 0\t" << stats.invalid_chi2pid_custom_pid_zero << endl;
    cout << "valid: chi2pid_custom != 99999.0 && Assigned PID != 0\t" << stats.valid_chi2pid_custom_nonzero_pid << endl;
    cout << "invalid: chi2pid_custom == 99999.0 && Assigned PID != 0\t" << stats.total_invalid_chi2pid_custom_nonzero_pid << endl;

    cout << "\nTotal events processed: " << total_events_processed << endl;
    cout << "Real time: " << timer.RealTime() << " s, CPU time: " << timer.CpuTime() << " s" << endl;

    // Mismatch details (unchanged)
    if (!mismatches.empty()) {
        cout << "****************************************************************************************" << endl;
        cout << "Mismatch Details:" << endl;
        for (const auto& mismatch : mismatches) {
            cout << "Event " << mismatch.eventNum << ", Particle Index " << mismatch.particleIdx
                 << ": REC::PID=" << mismatch.recPid << ", Assigned PID=" << mismatch.assignedPid
                 << ", Charge=" << mismatch.charge << ", Momentum=" << mismatch.momentum
                 << ", NPHE_HTCC=" << mismatch.npheHtcc << ", NPHE_LTCC=" << mismatch.npheLtcc
                 << ", Status=" << mismatch.status
                 << ", Chi2pid (REC)=" << mismatch.chi2pid
                 << ", Chi2pid_custom (Computed for Assigned)=" << mismatch.chi2pidCustom
                 << ", Diff=" << (mismatch.chi2pidCustom - mismatch.chi2pid)
                 << ", Has_FTOF=" << (mismatch.hasFtof ? "Yes" : "No")
                 << ", Has_CTOF=" << (mismatch.hasCtof ? "Yes" : "No") << endl;
        }
        cout << "Total mismatches: " << mismatches.size() << endl;
    } else {
        cout << "No mismatches found." << endl;
    }

    // Plotting code (unchanged)
    TCanvas *canvasChi2pid = new TCanvas("canvasChi2pid", "Chi2pid Comparison for Charged Particles", 1200, 600);
    canvasChi2pid->Divide(3, 1);
    TH1F *hChi2pidAll = new TH1F("hChi2pidAll", "Chi2pid (REC::Particle);Chi2pid;Counts", 100, -15, 15);
    TH1F *hChi2pidCustomAll = new TH1F("hChi2pidCustomAll", "Chi2pid_custom (Assigned);Chi2pid;Counts", 100, -15, 15);
    TH1F *hChi2pidDiff = new TH1F("hChi2pidDiff", "Chi2pid_custom - Chi2pid (Charged Particles);Difference;Counts", 100, -15, 15);
    hChi2pidAll->SetLineColor(kBlue);
    hChi2pidAll->SetFillColor(kBlue);
    hChi2pidAll->SetFillStyle(3004);
    hChi2pidCustomAll->SetLineColor(kRed);
    hChi2pidCustomAll->SetFillColor(kRed);
    hChi2pidCustomAll->SetFillStyle(3005);
    hChi2pidDiff->SetLineColor(kGreen);
    hChi2pidDiff->SetFillColor(kGreen);
    hChi2pidDiff->SetFillStyle(3006);

    for (size_t i = 0; i < chi2pidList.size() && i < chi2pidCustomList.size(); ++i) {
        hChi2pidAll->Fill(chi2pidList[i]);
        hChi2pidCustomAll->Fill(chi2pidCustomList[i]);
        hChi2pidDiff->Fill(chi2pidCustomList[i] - chi2pidList[i]);
    }

    canvasChi2pid->cd(1);
    hChi2pidAll->Draw("HIST");
    canvasChi2pid->cd(2);
    hChi2pidCustomAll->Draw("HIST");
    canvasChi2pid->cd(3);
    hChi2pidAll->Draw("HIST");
    hChi2pidCustomAll->Draw("HIST SAME");
    hChi2pidDiff->Draw("HIST SAME");
    TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg1->AddEntry(hChi2pidAll, "REC::Particle Chi2pid", "f");
    leg1->AddEntry(hChi2pidCustomAll, "Chi2pid_custom (Assigned PID)", "f");
    leg1->AddEntry(hChi2pidDiff, "Chi2pid_custom - Chi2pid", "f");
    leg1->Draw();

    cout << "Statistical Comparison for Charged Particles:" << endl;
    cout << "Mean Chi2pid (REC): " << hChi2pidAll->GetMean() << ", RMS: " << hChi2pidAll->GetRMS() << endl;
    cout << "Mean Chi2pid_custom (Assigned): " << hChi2pidCustomAll->GetMean() << ", RMS: " << hChi2pidCustomAll->GetRMS() << endl;
    cout << "Mean Difference (Chi2pid_custom - Chi2pid): " << hChi2pidDiff->GetMean() << ", RMS: " << hChi2pidDiff->GetRMS() << endl;

    canvasChi2pid->Print("output.pdf");
    cout << "Chi2pid histograms saved to output.pdf" << endl;

    // Momentum Plots (unchanged)
    TCanvas *canvasMomentum = new TCanvas("canvasMomentum", "Momentum Comparison for Charged Particles", 1200, 1200);
    canvasMomentum->Divide(3, 3);
    int pad = 1;
    for (int pid : particle_types) {
        canvasMomentum->cd(pad++);
        hMomentumRec[pid]->Draw("HIST");
        hMomentumCustom[pid]->Draw("HIST SAME");
        TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        string rec_label = "REC::PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hMomentumRec[pid]->GetEntries())) + ")";
        string custom_label = "Assigned PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hMomentumCustom[pid]->GetEntries())) + ")";
        leg->AddEntry(hMomentumRec[pid], rec_label.c_str(), "f");
        leg->AddEntry(hMomentumCustom[pid], custom_label.c_str(), "f");
        leg->Draw();
    }
    canvasMomentum->Print("momentum_plots.pdf");
    cout << "Momentum plots saved to momentum_plots.pdf" << endl;

    // Chi2pid Plots for Each Particle Type (unchanged)
    TCanvas *canvasChi2pidPerParticle = new TCanvas("canvasChi2pidPerParticle", "Chi2pid Comparison per Particle Type", 1200, 1200);
    canvasChi2pidPerParticle->Divide(3, 3);
    pad = 1;
    for (int pid : particle_types) {
        canvasChi2pidPerParticle->cd(pad++);
        if (hChi2pidRec[pid]->GetEntries() > 0 || hChi2pidCustom[pid]->GetEntries() > 0) {
            hChi2pidRec[pid]->Draw("HIST");
            hChi2pidCustom[pid]->Draw("HIST SAME");
            TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
            string rec_label = "REC::PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hChi2pidRec[pid]->GetEntries())) + ")";
            string custom_label = "Assigned PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hChi2pidCustom[pid]->GetEntries())) + ")";
            leg->AddEntry(hChi2pidRec[pid], rec_label.c_str(), "f");
            leg->AddEntry(hChi2pidCustom[pid], custom_label.c_str(), "f");
            leg->Draw();

            cout << "Chi2pid Comparison for PID " << pid << ":" << endl;
            cout << "  Mean Chi2pid (REC): " << hChi2pidRec[pid]->GetMean() << ", RMS: " << hChi2pidRec[pid]->GetRMS() << endl;
            cout << "  Mean Chi2pid_custom (Assigned): " << hChi2pidCustom[pid]->GetMean() << ", RMS: " << hChi2pidCustom[pid]->GetRMS() << endl;
        }
    }
    canvasChi2pidPerParticle->Print("chi2pid_per_particle_plots.pdf");
    cout << "Chi2pid per particle plots saved to chi2pid_per_particle_plots.pdf" << endl;

    // Delta T Plots for Each Particle Type (unchanged)
    TCanvas *canvasDeltaTPerParticle = new TCanvas("canvasDeltaTPerParticle", "Delta T Comparison per Particle Type", 1200, 1200);
    canvasDeltaTPerParticle->Divide(3, 3);
    pad = 1;
    for (int pid : particle_types) {
        canvasDeltaTPerParticle->cd(pad++);
        if (hDtRec[pid]->GetEntries() > 0 || hDtCustom[pid]->GetEntries() > 0) {
            hDtRec[pid]->Draw("HIST");
            hDtCustom[pid]->Draw("HIST SAME");
            TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
            string rec_label = "REC::PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hDtRec[pid]->GetEntries())) + ")";
            string custom_label = "Assigned PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hDtCustom[pid]->GetEntries())) + ")";
            leg->AddEntry(hDtRec[pid], rec_label.c_str(), "f");
            leg->AddEntry(hDtCustom[pid], custom_label.c_str(), "f");
            leg->Draw();

            cout << "Delta T Comparison for PID " << pid << ":" << endl;
            cout << "  Mean Delta T (REC): " << hDtRec[pid]->GetMean() << " ns, RMS: " << hDtRec[pid]->GetRMS() << " ns" << endl;
            cout << "  Mean Delta T (Assigned): " << hDtCustom[pid]->GetMean() << " ns, RMS: " << hDtCustom[pid]->GetRMS() << " ns" << endl;
        }
    }
    canvasDeltaTPerParticle->Print("delta_t_per_particle_plots.pdf");
    cout << "Delta T per particle plots saved to delta_t_per_particle_plots.pdf" << endl;

    // Delta T vs Momentum Plots (unchanged)
    TCanvas *canvasDtVsPRec = new TCanvas("canvasDtVsPRec", "Delta T vs Momentum (REC::Particle)", 1200, 1200);
    canvasDtVsPRec->Divide(3, 3);
    pad = 1;
    for (int pid : particle_types) {
        canvasDtVsPRec->cd(pad++);
        if (hDtVsPRec[pid]->GetEntries() > 0) {
            hDtVsPRec[pid]->Draw("COLZ");
            TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
            string rec_label = "REC::PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hDtVsPRec[pid]->GetEntries())) + ")";
            leg->AddEntry(hDtVsPRec[pid], rec_label.c_str(), "f");
            leg->Draw();
        }
    }
    canvasDtVsPRec->Print("delta_t_vs_p_rec_plots.pdf");
    cout << "Delta T vs Momentum (REC::Particle) plots saved to delta_t_vs_p_rec_plots.pdf" << endl;

    TCanvas *canvasDtVsPCustom = new TCanvas("canvasDtVsPCustom", "Delta T vs Momentum (Assigned PID)", 1200, 1200);
    canvasDtVsPCustom->Divide(3, 3);
    pad = 1;
    for (int pid : particle_types) {
        canvasDtVsPCustom->cd(pad++);
        if (hDtVsPCustom[pid]->GetEntries() > 0) {
            hDtVsPCustom[pid]->Draw("COLZ");
            TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
            string custom_label = "Assigned PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hDtVsPCustom[pid]->GetEntries())) + ")";
            leg->AddEntry(hDtVsPCustom[pid], custom_label.c_str(), "f");
            leg->Draw();
        }
    }
    canvasDtVsPCustom->Print("delta_t_vs_p_custom_plots.pdf");
    cout << "Delta T vs Momentum (Assigned PID) plots saved to delta_t_vs_p_custom_plots.pdf" << endl;

    // Chi2pid vs Momentum Plots (unchanged)
    TCanvas *canvasChi2pidVsPRec = new TCanvas("canvasChi2pidVsPRec", "Chi2pid vs Momentum (REC::Particle)", 1200, 1200);
    canvasChi2pidVsPRec->Divide(3, 3);
    pad = 1;
    for (int pid : particle_types) {
        canvasChi2pidVsPRec->cd(pad++);
        if (hChi2pidVsPRec[pid]->GetEntries() > 0) {
            hChi2pidVsPRec[pid]->Draw("COLZ");
            TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
            string rec_label = "REC::PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hChi2pidVsPRec[pid]->GetEntries())) + ")";
            leg->AddEntry(hChi2pidVsPRec[pid], rec_label.c_str(), "f");
            leg->Draw();
        }
    }
    canvasChi2pidVsPRec->Print("chi2pid_vs_p_rec_plots.pdf");
    cout << "Chi2pid vs Momentum (REC::Particle) plots saved to chi2pid_vs_p_rec_plots.pdf" << endl;

    TCanvas *canvasChi2pidVsPCustom = new TCanvas("canvasChi2pidVsPCustom", "Chi2pid vs Momentum (Assigned PID)", 1200, 1200);
    canvasChi2pidVsPCustom->Divide(3, 3);
    pad = 1;
    for (int pid : particle_types) {
        canvasChi2pidVsPCustom->cd(pad++);
        if (hChi2pidVsPCustom[pid]->GetEntries() > 0) {
            hChi2pidVsPCustom[pid]->Draw("COLZ");
            TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
            string custom_label = "Assigned PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hChi2pidVsPCustom[pid]->GetEntries())) + ")";
            leg->AddEntry(hChi2pidVsPCustom[pid], custom_label.c_str(), "f");
            leg->Draw();
        }
    }
    canvasChi2pidVsPCustom->Print("chi2pid_vs_p_custom_plots.pdf");
    cout << "Chi2pid vs Momentum (Assigned PID) plots saved to chi2pid_vs_p_custom_plots.pdf" << endl;

    // Clean up (unchanged)
    delete canvasChi2pid;
    delete hChi2pidAll;
    delete hChi2pidCustomAll;
    delete hChi2pidDiff;
    delete canvasMomentum;
    delete canvasChi2pidPerParticle;
    delete canvasDeltaTPerParticle;
    delete canvasDtVsPRec;
    delete canvasDtVsPCustom;
    delete canvasChi2pidVsPRec;
    delete canvasChi2pidVsPCustom;
    for (int pid : particle_types) {
        delete hMomentumRec[pid];
        delete hMomentumCustom[pid];
        delete hChi2pidRec[pid];
        delete hChi2pidCustom[pid];
        delete hDtRec[pid];
        delete hDtCustom[pid];
        delete hDtVsPRec[pid];
        delete hDtVsPCustom[pid];
        delete hChi2pidVsPRec[pid];
        delete hChi2pidVsPCustom[pid];
    }

    cout << "****************************************************************************************" << endl;

    return 0;
}