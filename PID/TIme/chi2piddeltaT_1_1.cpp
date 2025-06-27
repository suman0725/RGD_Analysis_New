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

using namespace std;
namespace fs = std::filesystem;

// Constants (matching Java code)
const float HTCC_PION_THRESHOLD = 4.9f;
const float LTCC_PION_THRESHOLD = 3.0f;
const int HTCC_DETECTOR = 15;
const int LTCC_DETECTOR = 16;

// Map type to hold indices from detector banks keyed by particle index.
typedef std::map<int, std::vector<int>> IndexMap;

// Function to build a map from a bank based on a given index variable (e.g., "pindex")
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

// Utility functions for energy, beta, and deltaT
float calculateEnergy(float p, float mass) { 
    return sqrt(p * p + mass * mass); 
}

float calculateBeta(float p, float mass) { 
    return p / calculateEnergy(p, mass); 
}

// Note: 29.9792458 cm/ns is the speed of light in these units.
float calculateDeltaT(float path, float time, float beta, float startTime) {
    return time - (path / (beta * 29.9792458f)) - startTime;
}

int main() {
    TStopwatch timer;
    timer.Start();

    // Read directories from a file (one per line)
    ifstream inputFile("directories.txt");
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
    inputFile.close();
    if (directories.empty()) {
        cerr << "Error: No directories found" << endl;
        return 1;
    }

    // Create an output ROOT file and two histograms for ΔT
    TFile *outputFile = new TFile("particles_data_analysis.root", "RECREATE");

    TH1F* h1_deltaT_pion_negative_wpid = new TH1F("h1_deltaT_pion_negative_wpid", "-ve pions); #DeltaT (ns); Counts", 100, -5, 5);
    TH1F* h1_deltaT_pion_negative_wopid = new TH1F("h1_deltaT_pion_negative_wopid", "-ve pions; #DeltaT (ns); Counts", 100, -5, 5);
    TH1F* h1_deltaT_kaon_negative_wpid = new TH1F("h1_deltaT_kaon_negative_wpid", "-ve kaons; #DeltaT (ns); Counts", 100, -5, 5);
    TH1F* h1_deltaT_kaon_negative_wopid = new TH1F("h1_deltaT_kaon_negative_wopid", "-ve kaons; #DeltaT (ns); Counts", 100, -5, 5);

    TH1F* h1_deltaT_deuteron_wpid = new TH1F("h1_deltaT_deuteron_wpid", "Deuteron; #DeltaT (ns); Counts", 100, -5, 5);
    TH1F* h1_deltaT_deuteron_wopid = new TH1F("h1_deltaT_deuteron_wopid", "Deuteron; #DeltaT (ns); Counts", 100, -5, 5);

    TH1F* h1_deltaT_proton_wpid = new TH1F("h1_deltaT_antiproton_wpid", "- 2212; #DeltaT (ns); Counts", 100, -5, 5);
    TH1F* h1_deltaT_proton_wopid = new TH1F("h1_deltaT_antiproton_wopid", "-2212; #DeltaT (ns); Counts", 100, -40,140 );


    // required plots for trigger electrons
    TH1F* h1_nphe_tele = new TH1F ("h1_nphe_tele", "nphe_trigger_electron; nphe; counts", 200, 0, 60);
    TH1F* h1_nphe_tele_1 = new TH1F ("h1_nphe_tele", "nphe_trigger_electron 0-10; nphe; counts", 200, 0, 10); 
    TH1F* h1_charge_tele = new TH1F ("h1_charge_tele", "charge_trigger_electron; charge; counts", 200, -2, 0); 
    TH1F* h1_status_tele = new TH1F("h1_status_tele", "status_trigger_electron; status; counts", 200, -4000, -2000 );
    TH1F* h1_vz_tele = new TH1F ("h1_vz_tele", "vz_trigger_electron; vz (cm); counts", 200, - 25, 10); 
    TH1F* h1_chi2pid_tele = new TH1F ("h1_chi2pid_tele", "chi2pid_trigger_electron; chi2pid; counts", 200, -10, 10);
    TH1F* h1_energy_pcal_tele = new TH1F ("h1_energy_pcal_tele", "energy_pcal_trigger_electron; enrgy_pcal; counts", 200, 0, 2);
    TH1F* h1_energy_pcal_tel_1 = new TH1F ("h1_energy_pcal_tele_1", "energy_pcal_trigger_electron 0-0.3 GeV; enrgy_pcal; counts", 200, 0, 0.3);
    TH1F* h1_sampling_fraction_tele = new TH1F ("h1_sampling_fraction_tele", "sampling_fraction_trigger_electron; sampling_fraction; counts", 200, 0, 0.6);


    // Define particle masses (in GeV)
    const float mass_pion = 0.13957f;
    const float mass_kaon = 0.49367f;
    const float mass_proton = 0.93827f;
    const float mass_electron = 0.000511f;
    const float mass_positron = 0.000511f;
    const float mass_deuteron = 1.87561f; 


    
    int total_events_processed = 0;
    int total_electron_number = 0; 
    int pid, charge; 
    float chi2pid; 
    const int maxEvents = 100000;
    int event_count = 0;
    static int count_pid = 0, count_manual = 0;
    int positron_count = 0; 
    int  total_electrons = 0; 
    int trigger_positron_index = -1; 
   

    //*********For first row of REC::Particle Bank
    int trigger_electron_count_before = 0;
    int trigger_electron_count_after = 0;

    
    int trigger_electron_index = -1;
    float max_energy = -1.0;
    int directory_count = 0; 
    int hipo_file_count = 0; 


    for (const auto& dir : directories) {
        directory_count++; 
        vector<string> hipoFiles;
       /*  for (const auto& entry : fs::directory_iterator(dir)) {
            if (entry.is_regular_file() && entry.path().extension() == ".hipo") {
                hipoFiles.push_back(entry.path().string());
            }
        } */
        for (const auto& entry : fs::directory_iterator(dir)) {
        if (entry.is_regular_file() && entry.path().filename() == "skim_run_018536.hipo") {
            hipoFiles.push_back(entry.path().string());
        } 
        }
        
        if (hipoFiles.empty()) {
            cout << "  No .hipo files found in directory: " << dir << endl;
            continue;
        }

        for (const auto& file : hipoFiles) {
            hipo_file_count++;
            cout << "  Opening file: " << file << endl;
            hipo::reader reader;
            reader.open(file.c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);
           
            while (reader.next() == true ) {
                 
                if (event_count >= maxEvents) break;
                event_count++;
               

                hipo::event event;
                reader.read(event);
                
                // Load banks for this event
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

                // Build maps for scintillator and cherenkov hits keyed by particle index.
                IndexMap scinMap = loadMapByIndex(SCIN, "pindex");
                IndexMap cherMap = loadMapByIndex(CHER, "pindex");
                IndexMap caloMap = loadMapByIndex(CALO, "pindex");

                float startTime = EVENT.getFloat("startTime", 0);
                //if (startTime < 0) continue;

                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                // Identify trigger electron (to be skipped in PID processing)
                bool has_trigger_electron = false;
                bool has_trigger_positron = false;
                
                //cout << "Number of particles in this event: " << PART.getRows() << endl;
                for (int i = 0; i < PART.getRows(); i++) {
                    // Get particle information
                    pid = PART.getInt("pid", i);
                    if (pid ==11){total_electrons++; }
                    if (pid == 11 && i == 0) {                       
                        trigger_electron_count_before++;

                        short status = PART.getShort("status", i);
                        chi2pid = PART.getFloat("chi2pid", i); 
                        float vz_tele = PART.getFloat("vz", i);
                        charge = PART.getByte("charge",i);
                        float px = PART.getFloat("px", i);
                        float py = PART.getFloat("py", i);
                        float pz = PART.getFloat("pz", i);
                        float p = sqrt(px*px + py*py + pz*pz);

                        int det_cher,nphe_cher, det_cal, layer_cal;
                        float energy_pcal, energy_ecin, energy_ecout, energy_total, sampling_fraction ; 
                        if (cherMap.find(i) != cherMap.end()){
                            // Process associated CHER rows for this particle
                            for (int iCherRow : cherMap[i]){
                                det_cher = CHER.getByte("detector", iCherRow);
                                nphe_cher = CHER.getFloat("nphe", iCherRow); 
                            }
                        }
                        if (caloMap.find(i) != caloMap.end()){
                            // Process associated CHER rows for this particle
                            for (int iCalRow : caloMap[i]){
                                det_cal = CALO.getByte("detector", iCalRow);
                                layer_cal = CALO.getByte("layer", iCalRow);
                                if (/* det_cal == 7 && */ layer_cal == 1)
                                    energy_pcal = CALO.getFloat("energy", iCalRow); 
                                else if ( layer_cal == 4) 
                                    energy_ecin = CALO.getFloat("energy", iCalRow); 
                                else if (layer_cal == 7)
                                    energy_ecout = CALO.getFloat("energy", iCalRow); 
                            }
                        }
                        energy_total = energy_pcal + energy_ecin + energy_ecout; 
                        sampling_fraction = energy_total / p ; 
                        
                        // Apply trigger electron cuts
                        if (  status < 0 &&   chi2pid > -5 && chi2pid < 5 && vz_tele >= -20 && vz_tele<=5 && abs(status)/1000 ==  2   ) {
                            trigger_electron_count_after++; 
                            has_trigger_electron = true;
                            trigger_electron_index = i ; 
                            // Count it if all cuts pass           
                            h1_status_tele->Fill(status);
                            h1_vz_tele->Fill(vz_tele); 
                            h1_chi2pid_tele->Fill(chi2pid); 
                            h1_nphe_tele->Fill(nphe_cher);
                            h1_nphe_tele_1->Fill(nphe_cher); 
                            h1_energy_pcal_tele->Fill(energy_pcal);
                            h1_energy_pcal_tel_1->Fill(energy_pcal);
                            h1_sampling_fraction_tele->Fill(sampling_fraction);  
                            h1_charge_tele->Fill(charge);                           
                        }
                       
                    }

                   /* if (pid == -11){

                        positron_count++; 
                        trigger_positron_index = i ;

                        has_trigger_positron = true; 
                    } */
                    
                }

                
                ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


                // Skip this event if no trigger electron found
                if (!has_trigger_electron) continue; // Skip event if no trigger electrons       
                total_events_processed++;

               for (int i = 0; i < PART.getRows(); i++) {
                    if (i == trigger_electron_index) continue;
                    //if (i ==  trigger_positron_index) continue; 
                    const int charge = PART.getByte("charge", i);
                    const int status = PART.getInt("status", i); 
                    if (abs(status)/1000 != 2) continue;
                    //if (abs(status) < 2000 && abs(status)>= 4000) continue;
                    if (charge >= 0) continue;
                    // Get momentum components and compute magnitude
                    const float px = PART.getFloat("px", i);
                    const float py = PART.getFloat("py", i);
                    const float pz = PART.getFloat("pz", i);
                    const float p = sqrt(px*px + py*py + pz*pz);
                    float vt = PART.getFloat("vt",i);
                    // --- Bank PID filling (independent of ΔT) ---
                    int bank_pid = PART.getInt("pid", i);
                    //if (bank_pid == 211){count_pid++;}
                    
                   /* if (scinMap.find(i) != scinMap.end()) {
                    // Variables to store the best FTOF hit (prioritized by layer)
                    int best_layer = -1;
                    float best_path = 0.0;
                    float best_time = 0.0;

                    // First pass: Find the highest priority FTOF hit (1B > 1A > 2)
                    for (int iScinRow : scinMap[i]) {
                        if (SCIN.getByte("detector", iScinRow) == 12) { // FTOF
                            int layer = SCIN.getByte("layer", iScinRow);

                            // Priority check (1B > 1A > 2)
                            if (layer == 2) { // 1B (highest priority)
                                best_layer = 2;
                                best_path = SCIN.getFloat("path", iScinRow);
                                best_time = SCIN.getFloat("time", iScinRow);
                                break; // Stop after finding 1B
                            } 
                            else if (layer == 1 && best_layer < 2) { // 1A
                                best_layer = 1;
                                best_path = SCIN.getFloat("path", iScinRow);
                                best_time = SCIN.getFloat("time", iScinRow);
                            } 
                            else if (layer == 3 && best_layer < 1) { // 2 (lowest)
                                best_layer = 3;
                                best_path = SCIN.getFloat("path", iScinRow);
                                best_time = SCIN.getFloat("time", iScinRow);
                            }
                        }
                    }

                    // Second pass: Use the best FTOF hit for ΔT calculation
                    if (best_layer != -1) {
                        float mass = mass_pion; // Default to pion
                        if (bank_pid == -321) mass = mass_kaon;
                        else if (bank_pid == -2212) mass = mass_proton; */
                        //else if (bank_pid == 45) mass = mass_deuteron;

                        float dt_bank = vt - startTime ; 
                        
                        if (bank_pid == - 211) {
                            h1_deltaT_pion_negative_wpid->Fill(dt_bank);
                            count_pid++;
                        }
                        if (bank_pid == -321) {
                            h1_deltaT_kaon_negative_wpid->Fill(dt_bank);
                        }
                        if (bank_pid == - 2212) {
                            h1_deltaT_proton_wpid->Fill(dt_bank);
                        }
                        /*if (bank_pid == 45) {
                            h1_deltaT_deuteron_wpid->Fill(dt_bank);
                        } */
                   // }
                //}

                  
                   // --- Manual PID ΔT calculation with FTOF panel prioritization (by LAYER) ---
                    float best_dt_pion = numeric_limits<float>::max();
                    float best_dt_kaon = numeric_limits<float>::max();
                    float best_dt_proton = numeric_limits<float>::max();
                    float best_dt_deuteron = numeric_limits<float>::max();
                    float raw_dt_pion = 0;  // Store ACTUAL ΔT with sign
                    bool valid_dt_found = false;
                    float dt_pion = 0 ; 
                    float dt_kaon = 0 ; 
                    float dt_proton = 0 ; 
                    float dt_deuteron = 0 ; 
                    //float dt_kaon = 0; 

                    int highestPriorityLayer = -1;
                    float selectedPath = 0.0, selectedTime = 0.0;
                   

                    // First pass: Look for FTOF hits
                    if (scinMap.find(i) != scinMap.end()) {
                        for (int iScinRow : scinMap[i]) {
                            if (SCIN.getByte("detector", iScinRow) == 12) { // FTOF hits
                                int layer = SCIN.getByte("layer", iScinRow);
                                float path = SCIN.getFloat("path", iScinRow);
                                float time = SCIN.getFloat("time", iScinRow);

                                // Priority check (1B > 1A > 2)
                                if (layer == 2) { // 1B (highest priority)
                                    highestPriorityLayer = 2;
                                    selectedPath = path;
                                    selectedTime = time;
                                    valid_dt_found = true;
                                    break; // Stop after finding 1B
                                } else if (layer == 1 && highestPriorityLayer < 2) { // 1A
                                    highestPriorityLayer = 1;
                                    selectedPath = path;
                                    selectedTime = time;
                                    valid_dt_found = true;
                                } else if (layer == 3 && highestPriorityLayer < 2) { // 2 (lowest priority)
                                    highestPriorityLayer = 3;
                                    selectedPath = path;
                                    selectedTime = time;
                                    valid_dt_found = true;
                                }
                            }
                        }
                    }

                    // If no valid FTOF hit found, use ECAL hits
                    if (!valid_dt_found && caloMap.find(i) != caloMap.end()) {
                        int bestEcalLayer = -1;
                        float bestPath = 0.0, bestTime = 0.0;

                        for (int iCaloRow : caloMap[i]) {
                            if (CALO.getByte("detector", iCaloRow) == 7) { // ECAL hits
                                int layer = CALO.getByte("layer", iCaloRow);

                                // Priority order: PCAL (1) > ECIN (4) > ECOUT (7)
                                if ((layer == 1 && bestEcalLayer < 1) ||  // PCAL (highest priority)
                                    (layer == 4 && bestEcalLayer < 4) ||  // ECIN (middle priority)
                                    (layer == 7 && bestEcalLayer < 7)) {  // ECOUT (lowest priority)

                                    bestEcalLayer = layer;
                                    bestPath = CALO.getFloat("path", iCaloRow);
                                    bestTime = CALO.getFloat("time", iCaloRow);
                                }
                            }
                        }

                        // If a valid ECAL hit is found, use this data for ΔT calculation
                        if (bestEcalLayer != -1) {
                            selectedPath = bestPath;
                            selectedTime = bestTime;
                            highestPriorityLayer = bestEcalLayer;
                            valid_dt_found = true;
                        }
                    }

                    

                   // Now calculate ΔT using the selected path and time only if the beta is positive
                    float beta_pion = calculateBeta(p, mass_pion);
                    float beta_kaon = calculateBeta(p, mass_kaon);
                    float beta_proton = calculateBeta(p, mass_proton);
                    float beta_deuteron = calculateBeta(p, mass_deuteron);

                   

                    // Only proceed if beta values are positive
                    if (beta_pion > 0) {
                        dt_pion = calculateDeltaT(selectedPath, selectedTime, beta_pion, startTime);
                    }

                    if (beta_kaon > 0) {
                        dt_kaon = calculateDeltaT(selectedPath, selectedTime, beta_kaon, startTime);
                    }

                    if (beta_proton > 0) {
                        dt_proton = calculateDeltaT(selectedPath, selectedTime, beta_proton, startTime);
                    }

                    if (beta_deuteron > 0) {
                        dt_deuteron = calculateDeltaT(selectedPath, selectedTime, beta_deuteron, startTime);
                    }

                    // Store absolute values for PID selection
                    best_dt_pion = abs(dt_pion);
                    best_dt_kaon = abs(dt_kaon);
                    best_dt_proton = abs(dt_proton);
                    best_dt_deuteron = abs(dt_deuteron);

                    int manual_selected_pid = -211;
                    bool passes_cherenkov = true;

                    // Choose PID based on minimal |ΔT|
                    float min_dt = dt_pion;
                    

                    if (best_dt_kaon < min_dt) {
                        manual_selected_pid = -321;
                        min_dt = best_dt_kaon;
                    }
                    if (best_dt_proton < min_dt) {
                        manual_selected_pid = -2212;
                        min_dt = best_dt_proton;
                    }
                    

                   if (manual_selected_pid != -211 && p > HTCC_PION_THRESHOLD) {
                        bool htcc_signal = false;
                        if (cherMap.find(i) != cherMap.end()) {
                            for (int iCherRow : cherMap[i]) {
                                if (CHER.getByte("detector", iCherRow) == HTCC_DETECTOR &&
                                    CHER.getFloat("nphe", iCherRow) > 2) {
                                    htcc_signal = true;
                                    break;
                                }
                            }
                        }
                        if (htcc_signal) {
                            manual_selected_pid = -211;
        
                        }
                    }
                    else if ((manual_selected_pid == -321 || manual_selected_pid == -2212) && p > LTCC_PION_THRESHOLD) {
                        bool hasLTCCSignal = false;
                        if (cherMap.find(i) != cherMap.end()) {
                            for (int iCherRow : cherMap[i]) {
                                if (CHER.getByte("detector", iCherRow) == LTCC_DETECTOR &&
                                    CHER.getFloat("nphe", iCherRow) > 2) {
                                    hasLTCCSignal = true;
                                    break;
                                }
                            }
                        }
                        if (hasLTCCSignal) {
                            manual_selected_pid = -211;
                        }
                        
                    }
                    
                    // If no valid hits found, skip the particle
                    if (!valid_dt_found) {
                        cout << "No valid hits found for particle ID: " << i << ", skipping..." << endl;
                        continue;
                    }

                    // Fill histogram with ACTUAL ΔT (with sign)
                    if (manual_selected_pid == -211 && passes_cherenkov ) {
                        h1_deltaT_pion_negative_wopid->Fill(dt_pion);  // Use raw ΔT with sign
                        count_manual++;
                    }
                    if (manual_selected_pid == -321 && passes_cherenkov ) {
                        h1_deltaT_kaon_negative_wopid->Fill(dt_kaon);  // Use raw ΔT with sign
                        //count_manual++;
                    }
                    if (manual_selected_pid == -2212 && passes_cherenkov ) {
                        h1_deltaT_proton_wopid->Fill(dt_proton);  // Use raw ΔT with sign
                        //count_manual++;
                    }

                }
            
                
            }//end event loop 
            if (event_count >= maxEvents) break;
        } // end files loop 
        if (event_count >= maxEvents) break;
    } // end directories loop

    // Draw the histograms for comparison.
    TCanvas *c1 = new TCanvas("c1", "DeltaT Comparison", 800, 600);
    c1->Divide(2,2); 
    c1->cd(1); 

    h1_deltaT_pion_negative_wopid->SetLineColor(kRed);
    h1_deltaT_pion_negative_wpid->SetLineColor(kBlue);
    //h1_deltaT_pion_negative_wopid->SetStats(0); 
    //h1_deltaT_pion_negative_wpid->SetStats(0); 
    h1_deltaT_pion_negative_wopid->Draw();
    h1_deltaT_pion_negative_wpid->Draw("SAME");

   TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
   leg->SetTextSize(0.02); 
    leg->AddEntry(h1_deltaT_pion_negative_wopid, 
                Form("Custom PID (%.0f)", h1_deltaT_pion_negative_wopid->GetEntries()), 
                "l");
    leg->AddEntry(h1_deltaT_pion_negative_wpid, 
                Form("EB PID (%.0f)", h1_deltaT_pion_negative_wpid->GetEntries()), 
                "l");
    leg->Draw();

    c1->cd(2); 
    h1_deltaT_kaon_negative_wopid->SetLineColor(kRed);
    h1_deltaT_kaon_negative_wpid->SetLineColor(kBlue);
    //h1_deltaT_pion_negative_wopid->SetStats(0); 
    //h1_deltaT_pion_negative_wpid->SetStats(0); 
    h1_deltaT_kaon_negative_wopid->Draw();
    h1_deltaT_kaon_negative_wpid->Draw("SAME");

   TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
   leg1->SetTextSize(0.02); 
    leg1->AddEntry(h1_deltaT_kaon_negative_wopid, 
                Form("Custom PID (%.0f)", h1_deltaT_kaon_negative_wopid->GetEntries()), 
                "l");
    leg1->AddEntry(h1_deltaT_kaon_negative_wpid, 
                Form("EB PID (%.0f)", h1_deltaT_kaon_negative_wpid->GetEntries()), 
                "l");
    leg1->Draw();

    c1->cd(3);
    h1_deltaT_proton_wopid->SetLineColor(kRed);
    h1_deltaT_proton_wpid->SetLineColor(kBlue);
    h1_deltaT_proton_wopid->Draw();
    h1_deltaT_proton_wpid->Draw("SAME");

    TLegend *leg2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg2->SetTextSize(0.02); 
    leg2->AddEntry(h1_deltaT_proton_wopid, 
                Form("Custom PID (%.0f)", h1_deltaT_proton_wopid->GetEntries()), 
                "l");
    leg2->AddEntry(h1_deltaT_proton_wpid, 
                Form("EB PID (%.0f)", h1_deltaT_proton_wpid->GetEntries()), 
                "l");
    leg2->Draw();

    /* // Deuteron comparison
    c1->cd(4);
    h1_deltaT_deuteron_wopid->SetLineColor(kRed);
    h1_deltaT_deuteron_wpid->SetLineColor(kBlue);
    h1_deltaT_deuteron_wopid->Draw();
    h1_deltaT_deuteron_wpid->Draw("SAME");

    TLegend *leg3 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg3->SetTextSize(0.03);  
    leg3->AddEntry(h1_deltaT_deuteron_wopid, 
                Form("(%.0f)", h1_deltaT_deuteron_wopid->GetEntries()), 
                "l");
    leg3->AddEntry(h1_deltaT_deuteron_wpid, 
                Form("EB PID (%.0f)", h1_deltaT_deuteron_wpid->GetEntries()), 
                "l");
    leg3->Draw();*/
    c1->SaveAs("deltaT_comparison_negative_018536.pdf");
    gStyle->SetStatFont(42);        // 42 is the Helvetica font which is generally clear
    gStyle->SetStatFontSize(0.03);  

    TCanvas *c2 = new TCanvas("c2", "Criteria and cuts to select trigger electron", 800, 600); 
    c2->Divide (3,3); 
    c2->cd(1); 
    h1_charge_tele->Draw(); 
    c2->cd(2); 
    h1_status_tele->Draw(); 
    c2->cd(3); 
    h1_vz_tele->Draw(); 
    c2->cd(4); 
    h1_chi2pid_tele->Draw(); 
    c2->cd(5); 
    h1_nphe_tele->Draw(); 
    c2->cd(6); 
    h1_nphe_tele_1->Draw(); 
    c2->cd(7); 
    h1_energy_pcal_tele->Draw(); 
    c2->cd(8); 
    h1_energy_pcal_tel_1->Draw(); 
    c2->cd(9);
    h1_sampling_fraction_tele->Draw(); 

    c2->SaveAs("All_cuts_required_select_trigger_electron_018536.pdf");

    TCanvas *c3 = new TCanvas("c3", "DeltaT Comparison deuteron", 800, 600);
    c3->Divide(2, 2);

    // First plot
    c3->cd(1);
    //h1_deltaT_deuteron_wpid->GetXaxis()->SetRangeUser(-10, 10); // Set range before drawing
    h1_deltaT_proton_wpid->Draw(); 

    // Second plot
    c3->cd(2);
    //h1_deltaT_deuteron_wopid->GetXaxis()->SetAsixRange(-50, 50); // Set range before drawing
    h1_deltaT_proton_wopid->Draw();
    c3->cd(3);
    h1_deltaT_proton_wopid->Draw(); 
    h1_deltaT_proton_wpid->Draw("SAME"); 
    


    // Update canvas to apply the changes
    c3->Update(); 

    // Save the canvas
    c3->SaveAs("deltaT_comparison_018536_antiproton.pdf");

    gStyle->SetStatFont(42);        // 42 is the Helvetica font which is generally clear
    gStyle->SetStatFontSize(0.03);    // Adjust font size as needed

    


    outputFile->Write();
    outputFile->Close();

    timer.Stop();
    cout << "****************************************************************************************" << endl; 
    cout << "Total number of events: " << event_count << endl; 
   
    cout << "Entries: PID Event Builder = " << count_pid << ", Custom PID  = " << count_manual << endl;
    cout <<"Total electrons: " <<  total_electrons << endl; 
    cout << "Total first row elctrons before cut: " << trigger_electron_count_before << endl; 
    cout << "Total first row elctrons after cut: " << trigger_electron_count_after << endl; 
    cout << "Total events processed after selecting the trigger electrons after all cuts: " << total_events_processed << endl;
    cout << "Real time: " << timer.RealTime() << " s, CPU time: " << timer.CpuTime() << " s" << endl;
    cout << "****************************************************************************************" << endl; 

    
  
    return 0;
}
