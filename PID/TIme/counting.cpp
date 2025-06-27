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
#include "TH2F.h"
#include <cstddef>
#include <TStopwatch.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <limits>
#include "TStyle.h"
#include <unordered_map>

// Particle mass constants (GeV/c^2)
const float MASS_D  = 1.8756;  // Deuteron
const float MASS_P  = 0.9383;  // Proton
const float MASS_PI = 0.1396;  // Pion
const float MASS_K  = 0.4937;  // Kaon

int pion_count = 0;
int pion_count_eb = 0;
int total_electrons = 0;
int trigger_electron_count_before = 0;
int trigger_electron_count_after = 0;
int total_events_processed = 0;

// Map type to hold indices from detector banks keyed by particle index.
typedef std::map<int, std::vector<int>> IndexMap;

// Function to build a map from a bank based on a given index variable (e.g., "pindex")
IndexMap loadMapByIndex(hipo::bank& fromBank, const char* idxVarName) {
    IndexMap map;
    for (int iFrom = 0; iFrom < fromBank.getRows(); ++iFrom) {
        int iTo = fromBank.getInt(idxVarName, iFrom);
        map[iTo].push_back(iFrom);
    }
    return map;
}

// Compute β(p) for a given mass hypothesis
float betaFromMomentum(float p, float mass) {
    return p / sqrt(p * p + mass * mass);
}

// Assign PID based on minimizing |Δt|
int assignPID(float p, float beta_meas, float path, float time, float start_time) {
    struct Hypothesis {
        int pid;
        float mass;
    };

    std::vector<Hypothesis> hypotheses = {
        {45, MASS_D},  {2212, MASS_P}, {211, MASS_PI}, {321, MASS_K}
    };

    int best_pid = 0;
    float min_dt = std::numeric_limits<float>::max();

    for (const auto& hyp : hypotheses) {
        float beta_calc = betaFromMomentum(p, hyp.mass);
        float expected_time = path / (beta_calc * 29.9792);
        float dt = fabs(start_time - (time - expected_time));

        if (dt < min_dt) {
            min_dt = dt;
            best_pid = hyp.pid;
        }
    }
    return best_pid;
}

namespace fs = std::filesystem;
using namespace std;

int main() {
    fstream inputFile("directories.txt");
    if (!inputFile.is_open()) {
        cerr << "Error: Could not open directories.txt" << endl;
        return 1;
    }
    vector<string> directories;
    string dir;
    while (getline(inputFile, dir)) {
        if (!dir.empty()) {
            directories.push_back(dir);
        }
    }

    int maxEvents = 10000;
    int counter = 0;

    for (const auto& dir : directories) {
        vector<string> hipoFiles;
        /* for (const auto& entry : fs::directory_iterator(dir)) {
            if (entry.is_regular_file() && entry.path().extension() == ".hipo") {
                hipoFiles.push_back(entry.path().string());
            }
        } */
        for (const auto& entry : fs::directory_iterator(dir)) {
        if (entry.is_regular_file() && entry.path().filename() == "skim_run_018537.hipo") {
            hipoFiles.push_back(entry.path().string());
        }

        for (const auto& file : hipoFiles) {
            hipo::reader reader;
            reader.open(file.c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);

            while (reader.next() == true) {
                if (counter >= maxEvents) break;
                counter++;

                hipo::event event;
                reader.read(event);

                hipo::bank PART(factory.getSchema("REC::Particle"));
                hipo::bank EVENT(factory.getSchema("REC::Event"));
                hipo::bank SCIN(factory.getSchema("REC::Scintillator"));

                event.getStructure(PART);
                event.getStructure(EVENT);
                event.getStructure(SCIN);

                float start_time = EVENT.getFloat("startTime", 0);
                IndexMap scinMap = loadMapByIndex(SCIN, "pindex");

                // Identify trigger electron
                bool has_trigger_electron = false;
                int trigger_electron_index = -1;
                
                for (int i = 0; i < PART.getRows(); i++) {
                    int pid = PART.getInt("pid", i);
                     short status = PART.getShort("status", i);
                     if (abs(status)/1000 == 8){cout << "status: " << status << endl; }
                    if (pid == 11) total_electrons++;
                    if (pid == 11 && i == 0) {
                        trigger_electron_count_before++;}

                       
                      /*   float chi2pid = PART.getFloat("chi2pid", i);
                        float vz_tele = PART.getFloat("vz", i);
                        float px = PART.getFloat("px", i);
                        float py = PART.getFloat("py", i);
                        float pz = PART.getFloat("pz", i);
                        float p = sqrt(px * px + py * py + pz * pz);

                        if (status < 0 && chi2pid > -5 && chi2pid < 5 && vz_tele >= -20 && vz_tele <= 5 && abs(status)/1000 == 2) {
                            trigger_electron_count_after++;
                            has_trigger_electron = true;
                            trigger_electron_index = i;
                        } */

                       
                    }
                }

                /* if (!has_trigger_electron) continue; // Skip event if no trigger electron
                total_events_processed++;

                for (int i = 0; i < PART.getRows(); i++) {
                    float px = PART.getFloat("px", i);
                    float py = PART.getFloat("py", i);
                    float pz = PART.getFloat("pz", i);
                    float beta_meas = PART.getFloat("beta", i);
                    float p = sqrt(px * px + py * py + pz * pz);
                    int charge = PART.getInt("charge", i);
                    int pid = PART.getInt("pid", i); 
                    if (pid == 211){pion_count_eb++ ; }

                    int best_layer = 999;  // Initialize with a high number (no layer chosen yet)
                    float path = 0.0, time = 0.0;
                    bool valid_dt_found = false;


                    if (scinMap.find(i) != scinMap.end()) {
                        for (int iScinRow : scinMap[i]) {  
                            int detector = SCIN.getByte("detector", iScinRow);
                            int layer = SCIN.getByte("layer", iScinRow);  // Get FTOF layer

                            if (detector == 12) { // FTOF detector
                                if (layer == 2) {  // FTOF 1B (highest priority)
                                    best_layer = 2;
                                    path = SCIN.getFloat("path", iScinRow);
                                    time = SCIN.getFloat("time", iScinRow);
                                    valid_dt_found = true;
                                    break;  // Stop checking further since 1B is the best option
                                } else if (layer == 1 && best_layer != 2) { // FTOF 1A (second priority)
                                    best_layer = 1;
                                    path = SCIN.getFloat("path", iScinRow);
                                    time = SCIN.getFloat("time", iScinRow);
                                    valid_dt_found = true;
                                } else if (layer == 3 && best_layer != 2 && best_layer != 1) { // FTOF 2 (lowest priority)
                                    best_layer = 3;
                                    path = SCIN.getFloat("path", iScinRow);
                                    time = SCIN.getFloat("time", iScinRow);
                                    valid_dt_found = true;
                                }
                            }
                        }
                    }
                    if (valid_dt_found && charge > 0) {
                        int assigned_pid = assignPID(p, beta_meas, path, time, start_time);  
                        if (assigned_pid == 211){pion_count++;}                             
                    }
                } */
            }
            if (counter >= maxEvents) break;
        }
        if (counter >= maxEvents) break;
    }
       // cout << "Total number of pions using Event Builder(PID = 211): " << pion_count_eb << endl;
        //cout << "Total number of pions custom (PID = 211): " << pion_count << endl;
    
        return 0;

}
