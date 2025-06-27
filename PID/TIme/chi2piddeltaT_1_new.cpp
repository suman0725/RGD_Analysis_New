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
#include <set>

using namespace std;
namespace fs = std::filesystem;


std::set<int> unique_pids;

// Constants for HTCC&LTCC
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
    TH1F* h1_deltaT_pion_positive_wpid = new TH1F("h1_deltaT_pion_positive_measured", "+ve pions); #DeltaT (ns); Counts", 100, -5, 5);
    TH1F* h1_deltaT_pion_positive_wopid = new TH1F("h1_deltaT_pion_positive_theoretical", "+ve pions; #DeltaT (ns); Counts", 100, -5, 5);
    TH1F* h1_deltaT_kaon_positive_wpid = new TH1F("h1_deltaT_kaon_positive_measured", "+ve kaons; #DeltaT (ns); Counts", 100, -5, 5);
    TH1F* h1_deltaT_kaon_positive_wopid = new TH1F("h1_deltaT_kaon_positive_theoretical", "+ve kaons; #DeltaT (ns); Counts", 100, -5, 5);

    TH1F* h1_deltaT_deuteron_wpid = new TH1F("h1_deltaT_deuteron_measured", "Deuteron; #DeltaT (ns); Counts", 100, -5, 5);
    TH1F* h1_deltaT_deuteron_wopid = new TH1F("h1_deltaT_deuteron_theoretical", "Deuteron; #DeltaT (ns); Counts", 100, -40, 140);

    TH1F* h1_deltaT_proton_wpid = new TH1F("h1_deltaT_proton_measured", "Proton; #DeltaT (ns); Counts", 100, -5, 5);
    TH1F* h1_deltaT_proton_wopid = new TH1F("h1_deltaT_proton_theoretical", "Proton; #DeltaT (ns); Counts", 100, -5, 5);


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
    int count_all_positve = 0; 
    int count_positive_pion = 0; 
    int count_positive_kaon = 0; 
    int count_proton = 0 ; 
    int count_deuteron = 0; 
    int count_positron = 0 ; 
    int count_pidzero = 0; 
    int positron_count = 0; 
   
    int all_sum = 0; 
    int  total_electrons = 0; 
    int charge_number = 0; 

    int bank_pid_count = 0;
    int manual_pid_count = 0;

   

    //*********For first row of REC::Particle Bank
    int trigger_electron_count_before = 0;
    int trigger_electron_count_after = 0;

    
    int trigger_electron_index = -1;
    int trigger_positron_index = -1; 
    float max_energy = -1.0;
    int directory_count = 0; 
    int hipo_file_count = 0; 
    int total_positive_pions = 0;
    int total_negative_pions = 0;

    for (const auto& dir : directories) {
        directory_count++; 
        vector<string> hipoFiles;
        /* for (const auto& entry : fs::directory_iterator(dir)) {
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
                //PART.show(); 
               // SCIN.show(); 
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
                    if (pid == 211){total_positive_pions++;}
                    if (pid == -211){total_negative_pions++;}
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
                            // Process associated CALOW rows for this particle
                            for (int iCalRow : caloMap[i]){
                                det_cal = CALO.getByte("detector", iCalRow);
                                layer_cal = CALO.getByte("layer", iCalRow);
                                if ( layer_cal == 1)
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

                    if (pid == -11){

                        positron_count++; 
                        trigger_positron_index = i ;

                        has_trigger_positron = true; 
                    }
                    
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
                    if (charge <= 0) continue;
                    
                    count_all_positve++; 
                    // Get momentum components and compute magnitude
                    const float px = PART.getFloat("px", i);
                    const float py = PART.getFloat("py", i);
                    const float pz = PART.getFloat("pz", i);
                    const float p = sqrt(px*px + py*py + pz*pz);
                    float vt = PART.getFloat("vt",i);
                    // --- Bank PID filling (independent of ΔT) ---
                    int bank_pid = PART.getInt("pid", i);
                    unique_pids.insert(bank_pid);
                    
            

                    float dt_bank = vt - startTime ; 
                    
                    if (bank_pid == 211) {
                        h1_deltaT_pion_positive_wpid->Fill(dt_bank);
                        count_positive_pion++;
                    }
                    else if (bank_pid == 321) {
                        count_positive_kaon++; 
                        h1_deltaT_kaon_positive_wpid->Fill(dt_bank);
                    }
                    else if (bank_pid == 2212) {
                        count_proton++; 
                        h1_deltaT_proton_wpid->Fill(dt_bank);
                    }
                    else if (bank_pid == 45) {
                        count_deuteron++; 
                        h1_deltaT_deuteron_wpid->Fill(dt_bank);
                    }
                    else if (bank_pid == -11){
                        count_positron++; 

                    }
                    else if (bank_pid == 0){
                        //std::cout << "PID 0 found: Charge = " << charge << " PID = " << bank_pid << std::endl;
                        count_pidzero++;
                    }

                    
                    bool valid_dt_found = false;
                    

                    // FTOF panel prioritization: 1B (layer=2) > 1A (layer=1) > 2 (layer=3)
                    

                    // Debugging: Print the particle ID to track which particle we are processing
                    cout << "Processing particle ID: " << i << endl;

                    int highestPriorityLayer = -1;
                    float selectedPath = 0.0, selectedTime = 0.0;
                    std::string selectedDetector = "";  // To store the selected detector type
                   

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
                                    selectedDetector = "FTOF";
                                    break; // Stop after finding 1B
                                } else if (layer == 1 && highestPriorityLayer < 2) { // 1A
                                    highestPriorityLayer = 1;
                                    selectedPath = path;
                                    selectedTime = time;
                                    selectedDetector = "FTOF";
                                    valid_dt_found = true;
                                } else if (layer == 3 && highestPriorityLayer < 2) { // 2 (lowest priority)
                                    highestPriorityLayer = 3;
                                    selectedPath = path;
                                    selectedTime = time;
                                    selectedDetector = "FTOF";
                                    valid_dt_found = true;
                                }
                            }
                        }
                    }

                    // If no valid FTOF hit found, use ECAL hits
                   /*  if (!valid_dt_found && caloMap.find(i) != caloMap.end()) {
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
                            selectedDetector = "ECAL";
                            valid_dt_found = true;

                        
                        }
                    }  */

                    // Print the final selected layer (FTOF or ECAL)
                    if (valid_dt_found) {
                        float mass = mass_pion; // Default to pion
                        if (bank_pid == 321) mass = mass_kaon;
                        else if (bank_pid == 2212) mass = mass_proton;
                        else if (bank_pid == 45) mass = mass_deuteron; 

                        float dt_measured= calculateDeltaT(selectedPath, selectedTime, calculateBeta(p, mass), startTime);

                        // Fill histogram with ACTUAL ΔT (with sign)
                        if (bank_pid == 211  ) {
                            h1_deltaT_pion_positive_wopid->Fill(dt_measured);  
                            count_manual++;
                        }
                        if (bank_pid == 321 ) {
                            h1_deltaT_kaon_positive_wopid->Fill(dt_measured);  
                            //count_manual++;
                        }
                        if (bank_pid == 2212  ) {
                            h1_deltaT_proton_wopid->Fill(dt_measured); 
                            //count_manual++;
                        }
                        if (bank_pid == 45 ) {
                        h1_deltaT_deuteron_wopid->Fill(dt_measured);  //
                            //count_manual++;
                        }

                    }
                }
                //all_sum = count_positive_pion + count_positive_kaon + count_proton + count_deuteron + count_positron + count_pidzero ; 
                all_sum = count_positive_pion + count_positive_kaon + count_proton + count_deuteron + count_positron  ; 
                 
                            
                
            }//end event loop 
            if (event_count >= maxEvents) break;
        } // end files loop 
        if (event_count >= maxEvents) break;
    } // end directories loop

    // Draw the histograms for comparison.
    TCanvas *c1 = new TCanvas("c1", "DeltaT Comparison", 800, 600);
    c1->Divide(2,2); 
    c1->cd(1); 

    h1_deltaT_pion_positive_wopid->SetLineColor(kRed);
    h1_deltaT_pion_positive_wpid->SetLineColor(kBlue);
    //h1_deltaT_pion_positive_wopid->SetStats(0); 
    //h1_deltaT_pion_positive_wpid->SetStats(0); 
    h1_deltaT_pion_positive_wopid->Draw();
    h1_deltaT_pion_positive_wpid->Draw("SAME");

    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->SetTextSize(0.02); 
    leg->AddEntry(h1_deltaT_pion_positive_wopid, 
                Form("Theoretical (%.0f)", h1_deltaT_pion_positive_wopid->GetEntries()), 
                "l");
    leg->AddEntry(h1_deltaT_pion_positive_wpid, 
                Form("Measured (%.0f)", h1_deltaT_pion_positive_wpid->GetEntries()), 
                "l");
    

    leg->Draw();

    c1->cd(2); 
    h1_deltaT_kaon_positive_wopid->SetLineColor(kRed);
    h1_deltaT_kaon_positive_wpid->SetLineColor(kBlue);
    //h1_deltaT_pion_positive_wopid->SetStats(0); 
    //h1_deltaT_pion_positive_wpid->SetStats(0); 
    h1_deltaT_kaon_positive_wopid->Draw();
    h1_deltaT_kaon_positive_wpid->Draw("SAME");

   TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
   leg1->SetTextSize(0.02); 
    leg1->AddEntry(h1_deltaT_kaon_positive_wopid, 
                Form("Theoretical (%.0f)", h1_deltaT_kaon_positive_wopid->GetEntries()), 
                "l");
    leg1->AddEntry(h1_deltaT_kaon_positive_wpid, 
                Form("Measured (%.0f)", h1_deltaT_kaon_positive_wpid->GetEntries()), 
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
                Form("Theoretical (%.0f)", h1_deltaT_proton_wopid->GetEntries()), 
                "l");
    leg2->AddEntry(h1_deltaT_proton_wpid, 
                Form("Measured (%.0f)", h1_deltaT_proton_wpid->GetEntries()), 
                "l");
    leg2->Draw();

    // Deuteron comparison
    c1->cd(4);
    h1_deltaT_deuteron_wopid->SetLineColor(kRed);
    h1_deltaT_deuteron_wpid->SetLineColor(kBlue);
    h1_deltaT_deuteron_wopid->Draw();
    h1_deltaT_deuteron_wpid->Draw("SAME");

    TLegend *leg3 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg3->SetTextSize(0.02);  
    leg3->AddEntry(h1_deltaT_deuteron_wopid, 
                Form("Theoretical (%.0f)", h1_deltaT_deuteron_wopid->GetEntries()), 
                "l");
    leg3->AddEntry(h1_deltaT_deuteron_wpid, 
                Form("Measured (%.0f)", h1_deltaT_deuteron_wpid->GetEntries()), 
                "l");
    leg3->Draw();

    c1->SaveAs("deltaT_comparison_018536_measured.pdf");
    //c1->SaveAs("deltaT_comparison_018537.pdf");
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
    //c2->SaveAs("All_cuts_required_select_trigger_electron_018537.pdf");

    TCanvas *c3 = new TCanvas("c3", "DeltaT Comparison deuteron", 800, 600);
    c3->Divide(2, 2);

    // First plot
    c3->cd(1);
    //h1_deltaT_deuteron_wpid->GetXaxis()->SetRangeUser(-10, 10); // Set range before drawing
    h1_deltaT_deuteron_wpid->Draw(); 

    // Second plot
    c3->cd(2);
    //h1_deltaT_deuteron_wopid->GetXaxis()->SetAsixRange(-50, 50); // Set range before drawing
    h1_deltaT_deuteron_wopid->Draw();
    c3->cd(3);
    h1_deltaT_deuteron_wopid->Draw(); 
    h1_deltaT_deuteron_wpid->Draw("SAME"); 
    


    // Update canvas to apply the changes
    c3->Update(); 

    // Save the canvas
    c3->SaveAs("deltaT_comparison_018536_deuteron_measured.pdf");



    gStyle->SetStatFont(42);        // 42 is the Helvetica font which is generally clear
    gStyle->SetStatFontSize(0.03);    // Adjust font size as needed



    outputFile->Write();
    outputFile->Close();

    timer.Stop();


    for (int pid : unique_pids) {
                    std::cout << "Unique PID found: " << pid << std::endl;
                }
    cout << "****************************************************************************************" << endl; 
    cout << "Total number of events: " << event_count << endl; 
    cout << "Total number of +ve pions: " << total_positive_pions << endl; 
    cout << "Total number of -ve pions: " << total_negative_pions << endl; 

    cout << "****************************************************************************************" << endl; 
    cout << "Total number +ve charge particles after charge cut: " << count_all_positve << endl; 
    cout << "Total number of  +ve pions with pid 211: " << count_positive_pion << endl; 
    cout << "Total number of  +ve kaons with pid 321: " << count_positive_kaon << endl; 
    cout << "Total number of protons with pid 2212: " << count_proton << endl; 
    cout << "Total number of  deuterons with pid 45: " << count_deuteron << endl; 
    cout << "Total number of positrons with pid -11: " << count_positron << endl; 
    cout << "Total number of particles with pid 0: " << count_pidzero << endl;
    cout << "all sum +ve particles: " << all_sum << endl; 
    cout << "****************************************************************************************" << endl; 
    
    
    
    
    
    

   
    cout << "Entries: PID Event Builder = " << count_pid << ", Custom PID  = " << count_manual << endl;
    cout <<"Total electrons: " <<  total_electrons << endl; 
    cout << "Total positron: " << positron_count << endl; 
    cout << "Total first row elctrons before cut: " << trigger_electron_count_before << endl; 
    cout << "Total first row elctrons after cut: " << trigger_electron_count_after << endl; 
    cout << "Total events processed after selecting the trigger electrons after all cuts: " << total_events_processed << endl;
    cout << "Real time: " << timer.RealTime() << " s, CPU time: " << timer.CpuTime() << " s" << endl;
    cout << "****************************************************************************************" << endl; 

    
  
    return 0;
}
