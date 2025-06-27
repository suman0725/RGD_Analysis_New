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
#include <cstddef>
#include <TStopwatch.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <limits>
#include <algorithm>

using namespace std;
namespace fs = std::filesystem;

// Constants (matching Java code)
const float HTCC_PION_THRESHOLD = 4.9f;
const float LTCC_PION_THRESHOLD = 3.0f;
const int HTCC_DETECTOR = 15;
const int LTCC_DETECTOR = 16;
const vector<int> FTOF_LAYERS = {2, 1, 3}; // 1B (layer 2), 1A (layer 1), 2 (layer 3)

typedef std::map<int, std::vector<int>> IndexMap;

IndexMap loadMapByIndex(hipo::bank& fromBank, const char* idxVarName) {
    IndexMap map;
    if (fromBank.getRows() > 0) {
        for (int iFrom = 0; iFrom < fromBank.getRows(); ++iFrom) {
            int iTo = fromBank.getInt(idxVarName, iFrom);
            map[iTo].push_back(iFrom);
        }
    }
    return map;
}

float calculateEnergy(float p, float mass) { return sqrt(p * p + mass * mass); }
float calculateBeta(float p, float mass) { return p / calculateEnergy(p, mass); }
float calculateDeltaT(float path, float time, float beta, float startTime) {
    return time - (path / (beta * 29.9792458f)) - startTime;
}

int main() {
    TStopwatch timer;
    timer.Start();

    // Read directories
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

    TFile *outputFile = new TFile("pion_count_comparison.root", "RECREATE");
    TH1F* h_builtin = new TH1F("h_builtin", "Built-in PID Count;p [GeV];Counts", 
                             1, 2.15, 2.20);
    TH1F* h_custom = new TH1F("h_custom", "Custom PID Count;p [GeV];Counts", 
                            1, 2.15, 2.20);

    // Particle masses
    const float mass_pion = 0.13957f;
    const float mass_kaon = 0.49367f;
    const float mass_proton = 0.93827f;

    for (const auto& dir : directories) {
        vector<string> hipoFiles;
        for (const auto& entry : fs::directory_iterator(dir)) {
            if (entry.is_regular_file() && entry.path().extension() == ".hipo") {
                hipoFiles.push_back(entry.path().string());
            }
        }

        for (const auto& file : hipoFiles) {
            cout << "Processing: " << file << endl;
            hipo::reader reader;
            reader.open(file.c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);

            const int maxEvents = 10000;
            int eventCounter = 0;
            
            while (reader.next() && eventCounter++ < maxEvents) {
                hipo::event event;
                reader.read(event);
                
                hipo::bank PART(factory.getSchema("REC::Particle"));
                hipo::bank EVENT(factory.getSchema("REC::Event"));
                hipo::bank SCIN(factory.getSchema("REC::Scintillator")); 
                hipo::bank CHER(factory.getSchema("REC::Cherenkov"));

                event.getStructure(PART);
                event.getStructure(EVENT);
                event.getStructure(SCIN);
                event.getStructure(CHER);

                float startTime = EVENT.getFloat("startTime", 0);
                if (startTime < 0) continue;

                IndexMap scinMap = loadMapByIndex(SCIN, "pindex");
                IndexMap cherMap = loadMapByIndex(CHER, "pindex");

                // Find trigger electron
                int trigger_index = -1;
                for (int i = 0; i < PART.getRows(); i++) {
                    if (PART.getInt("pid", i) == 11 && 
                        PART.getInt("status", i) < 0 &&
                        PART.getFloat("chi2pid", i) > -5 && 
                        PART.getFloat("chi2pid", i) < 5) {
                        trigger_index = i;
                        break;
                    }
                }

                // Process particles
                for (int i = 0; i < PART.getRows(); i++) {
                    if (i == trigger_index) continue;

                    const int charge = PART.getByte("charge", i);
                    if (charge <= 0) continue;

                    // Calculate momentum
                    const float px = PART.getFloat("px", i);
                    const float py = PART.getFloat("py", i);
                    const float pz = PART.getFloat("pz", i);
                    const float p = sqrt(px*px + py*py + pz*pz);
                    if (p < 2.15 || p > 2.20) continue;

                    // Timing analysis
                    float best_dt_pion = numeric_limits<float>::max();
                    float best_dt_kaon = numeric_limits<float>::max();
                    float best_dt_proton = numeric_limits<float>::max();

                    if (scinMap.find(i) != scinMap.end()) {
                        for (int iScin : scinMap[i]) {
                            // Get detector and layer with proper types
                            const int detector = SCIN.getByte("detector", iScin);
                            const float layer_float = SCIN.getFloat("layer", iScin);
                            const int layer = static_cast<int>(layer_float);

                            if (detector == 12) { // FTOF
                                // Apply layer prioritization
                                if (find(FTOF_LAYERS.begin(), FTOF_LAYERS.end(), layer) != FTOF_LAYERS.end()) {
                                    const float path = SCIN.getFloat("path", iScin);
                                    const float time = SCIN.getFloat("time", iScin);
                                    
                                    const float beta_pion = calculateBeta(p, mass_pion);
                                    const float beta_kaon = calculateBeta(p, mass_kaon);
                                    const float beta_proton = calculateBeta(p, mass_proton);
                                    
                                    const float dt_pion = calculateDeltaT(path, time, beta_pion, startTime);
                                    const float dt_kaon = calculateDeltaT(path, time, beta_kaon, startTime);
                                    const float dt_proton = calculateDeltaT(path, time, beta_proton, startTime);
                                    
                                    // Track best delta-t values
                                    if (abs(dt_pion) < abs(best_dt_pion)) best_dt_pion = dt_pion;
                                    if (abs(dt_kaon) < abs(best_dt_kaon)) best_dt_kaon = dt_kaon;
                                    if (abs(dt_proton) < abs(best_dt_proton)) best_dt_proton = dt_proton;
                                }
                            }
                        }
                    }

                    // PID selection
                    int selected_pid = 211;
                    float min_dt = abs(best_dt_pion);
                    
                    // Check kaon hypothesis
                    if (abs(best_dt_kaon) < min_dt) {
                        selected_pid = 321;
                        min_dt = abs(best_dt_kaon);
                    }
                    
                    // Check proton hypothesis
                    if (abs(best_dt_proton) < min_dt) {
                        selected_pid = 2212;
                    }

                    // Cherenkov checks
                    bool passes_cherenkov = true;
                    
                    // Pion HTCC check
                    if (selected_pid == 211 && p > HTCC_PION_THRESHOLD) {
                        passes_cherenkov = false;
                        for (int iCher : cherMap[i]) {
                            if (CHER.getByte("detector", iCher) == HTCC_DETECTOR &&
                                CHER.getFloat("nphe", iCher) > 2) {
                                passes_cherenkov = true;
                                break;
                            }
                        }
                    }
                    
                    // Kaon/Proton LTCC veto
                    if ((selected_pid == 321 || selected_pid == 2212) && p > LTCC_PION_THRESHOLD) {
                        bool hasLTCC = false;
                        for (int iCher : cherMap[i]) {
                            if (CHER.getByte("detector", iCher) == LTCC_DETECTOR &&
                                CHER.getFloat("nphe", iCher) > 2) {
                                hasLTCC = true;
                                break;
                            }
                        }
                        if (hasLTCC) {
                            selected_pid = 211;
                            passes_cherenkov = false;
                        }
                    }

                    // Fill histograms
                    if (PART.getInt("pid", i) == 211) {
                        h_builtin->Fill(p);
                    }
                    if (selected_pid == 211 && passes_cherenkov) {
                        h_custom->Fill(p);
                    }
                }
            }
        }
    }

    // Save results
    TCanvas *c1 = new TCanvas("c1", "PID Comparison", 800, 600);
    h_builtin->SetLineColor(kBlue);
    h_custom->SetLineColor(kRed);
    
    h_builtin->Draw();
    h_custom->Draw("SAME");

    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(h_builtin, Form("Built-in (%d)", h_builtin->GetEntries()), "l");
    leg->AddEntry(h_custom, Form("Custom (%d)", h_custom->GetEntries()), "l");
    leg->Draw();

    c1->SaveAs("pid_comparison.png");
    outputFile->Write();
    outputFile->Close();

    timer.Stop();
    cout << "Real time: " << timer.RealTime() << " s, CPU time: " << timer.CpuTime() << " s" << endl;
    return 0;
}