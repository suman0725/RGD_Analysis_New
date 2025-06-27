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

using namespace std;
namespace fs = std::filesystem;

// Typedef for the index map
typedef std::map<int, std::vector<int>> IndexMap;

// Function to create a map from index values
IndexMap loadMapByIndex(hipo::bank& fromBank, const char* idxVarName) {
    IndexMap map;
    if (fromBank.getRows() > 0) {
        for (int iFrom = 0; iFrom < fromBank.getRows(); ++iFrom) {
            int iTo = fromBank.getInt(idxVarName, iFrom);
            if (map.find(iTo) == map.end()) {
                map[iTo] = vector<int>();
            }
            map[iTo].push_back(iFrom);
        }
    }
    return map;
}

// Function to calculate beta for a given particle
/* float calculateBeta(float p, float mass) {
    float E = sqrt(p * p + mass * mass);
    return p / E;
} */
float calculateEnergy(float p, float mass) {
    return sqrt(p * p + mass * mass);  // Energy: E = sqrt(p^2 + m^2)
}

float calculateBeta(float p, float mass) {
    return p / calculateEnergy(p, mass);  // Beta: β = p / E
}

// Function to calculate delta_t for a given particle
float calculateDeltaT(float path, float time, float beta, float startTime) {
    return time - (path / (beta * 29.9792458)) - startTime;
}

int main() {

    //Starting Timeer
    TStopwatch timer; 
    // Path to the text file containing directory paths
    string dirListFile = "directories.txt";
    ifstream inputFile(dirListFile);

    if (!inputFile.is_open()) {
        cerr << "Error: Could not open " << dirListFile << endl;
        return 1;
    }

    // Read directories from the file
    vector<string> directories;
    string dir;
    while (getline(inputFile, dir)) {
        if (!dir.empty()) {
            directories.push_back(dir);
        }
    }
    inputFile.close();

    if (directories.empty()) {
        cerr << "Error: No directories found in " << dirListFile << endl;
        return 1;
    }

    // Particle masses in GeV/c^2
   /*  const float mass_pion = 0.13957;
    const float mass_kaon = 0.49367;
    const float mass_proton = 0.93827;
    const float mass_electron = 0.000511; */

    // Create ROOT trees for storing delta t and momentum for different particles
  /*   TTree* tree_pion = new TTree("pion", "Pion+ Tree");
    TTree* tree_kaon = new TTree("kaon", "Kaon+ Tree");
    TTree* tree_proton = new TTree("proton", "Proton Tree") */;

    TTree* tree_pparticles = new TTree("pparticles", "Positive Particles Tree");


    // Variables to hold delta t and momentum for each particle type
   /*  float delta_t_pion, momentum_pion, beta_pion;
    float delta_t_kaon, momentum_kaon, beta_kaon;
    float delta_t_proton, momentum_proton, beta_proton;
 */
    float chi2pid, momentum, beta; 
    int pid, charge; 

    // Vectors to store Cherenkov hits dynamically
    vector<int> nphe; 
    vector<int> detector; 
    
    // Create branches for Positive Particles Tree only for chi2pid for now
    tree_pparticles->Branch("pid", &pid, "pid/I");
    tree_pparticles->Branch("charge", &charge, "charge/I"); 
    tree_pparticles->Branch("chi2pid", &chi2pid, "chi2pid/F");
    tree_pparticles->Branch("momentum", &momentum,"momentum/F");
    tree_pparticles->Branch("beta", &beta, "beta/F");
    tree_pparticles->Branch("nphe", &nphe);
    tree_pparticles->Branch("detector", &detector);

   /*  // Create branches for the trees
    tree_pion->Branch("delta_t", &delta_t_pion, "delta_t/F");
    tree_pion->Branch("momentum", &momentum_pion, "momentum/F"); 
    //tree_pion->Branch("beta", &beta_pion, "beta/F");

    tree_kaon->Branch("delta_t", &delta_t_kaon, "delta_t/F");
    tree_kaon->Branch("momentum", &momentum_kaon, "momentum/F");
    //tree_kaon->Branch("beta", &beta_kaon, "beta/F");

    tree_proton->Branch("delta_t", &delta_t_proton, "delta_t/F");
    tree_proton->Branch("momentum", &momentum_proton, "momentum/F");
   // tree_proton->Branch("beta", &beta_proton, "beta/F"); */

    /* TH2F* h2_momentum_deltaT_pion = new TH2F("h2_momentum_deltaT_pion", "Pion Momentum vs DeltaT;Momentum (GeV/c);Delta T (ns)", 100, 0, 10, 100, -10, 10);
    TH2F* h2_momentum_deltaT_kaon = new TH2F("h2_momentum_deltaT_kaon", "Kaon Momentum vs DeltaT;Momentum (GeV/c);Delta T (ns)", 100, 0, 10, 100, -10, 10);
    TH2F* h2_momentum_deltaT_proton = new TH2F("h2_momentum_deltaT_proton", "Proton Momentum vs DeltaT;Momentum (GeV/c);Delta T (ns)", 100, 0, 10, 100, -10, 10); */

    // chi2pid and momentum from REC::Particle bank
    TH2F* h2_momentum_chi2pid_positive = new TH2F("h2_momentum_chi2pid_positive", "Momentum vs chi2pid for positive particles;Momentum (GeV);chi2pid", 100, 0, 10, 100, -10, 10);
    TH2F* h2_momentum_chi2pid_negative = new TH2F("h2_momentum_chi2pid_negative", "Momentum vs chi2pid for negative particles;Momentum (GeV);chi2pid", 100, 0, 10, 100, -10, 10);

    // nphe from REC::Cherenkov Bank 
    // plot of momentum Vs number of photo electrons for charged hadrons for LTCC and HTCC
    TH2F* h2_momentum_nphe_positive = new TH2F("h2_momentum_nphe_positive", "Momentum vs nphe for positive particles LTCC+HTCC ;Momentum (GeV);nphe (counts)", 100, 0, 10, 100, 0, 100);

    TH2F* h2_momentum_nphe_positive_LTCC = new TH2F("h2_momentum_nphe_positive_LTCC", "Momentum vs nphe for positive particles LTCC ;Momentum (GeV);nphe (counts)", 100, 0, 10, 100, 0, 100);
    TH2F* h2_momentum_nphe_positive_pions_LTCC = new TH2F("h2_momentum_nphe_positive_pions_LTCC", "Momentum vs nphe for positive pions LTCC ;Momentum (GeV);nphe (counts)", 100, 0, 10, 100, 0, 100);
    TH2F* h2_momentum_nphe_positive_kaons_LTCC = new TH2F("h2_momentum_nphe_positive_kaons_LTCC", "Momentum vs nphe for positive kaons LTCC ;Momentum (GeV);nphe (counts)", 100, 0, 10, 100, 0, 100);

    TH2F* h2_momentum_nphe_positive_HTCC = new TH2F("h2_momentum_nphe_positive_HTCC", "Momentum vs nphe for positive particles HTCC ;Momentum (GeV);nphe (counts)", 100, 0, 10, 100, 0, 100);
    TH2F* h2_momentum_nphe_positive_pions_HTCC = new TH2F("h2_momentum_nphe_positive_pions_HTCC", "Momentum vs nphe for positive pions HTCC ;Momentum (GeV);nphe (counts)", 100, 0, 10, 100, 0, 100);
    TH2F* h2_momentum_nphe_positive_kaons_HTCC = new TH2F("h2_momentum_nphe_positive_kaons_HTCC", "Momentum vs nphe for positive kaons HTCC ;Momentum (GeV);nphe (counts)", 100, 0, 10, 100, 0, 100);

    TH2F* h2_momentum_nphe_negative = new TH2F("h2_momentum_nphe_negative", "Momentum vs nphe for negative particles LTCC+HTCC ;Momentum (GeV);nphe(counts)", 100, 0, 10, 100, 0, 100);

    TH2F* h2_momentum_nphe_negative_LTCC = new TH2F("h2_momentum_nphe_negative_LTCC", "Momentum vs nphe for negative particles LTCC ;Momentum (GeV);nphe(counts)", 100, 0, 10, 100, 0, 100);
    TH2F* h2_momentum_nphe_negative_pions_LTCC = new TH2F("h2_momentum_nphe_negative_pions_LTCC", "Momentum vs nphe for negative pions LTCC ;Momentum (GeV);nphe(counts)", 100, 0, 10, 100, 0, 100);
    TH2F* h2_momentum_nphe_negative_kaons_LTCC = new TH2F("h2_momentum_nphe_negative_kaons_LTCC", "Momentum vs nphe for negative kaons LTCC ;Momentum (GeV);nphe(counts)", 100, 0, 10, 100, 0, 100);

    TH2F* h2_momentum_nphe_negative_HTCC = new TH2F("h2_momentum_nphe_negative_HTCC", "Momentum vs nphe for negative particles HTCC ;Momentum (GeV);nphe(counts)", 100, 0, 10, 100, 0, 100);
    TH2F* h2_momentum_nphe_negative_pions_HTCC = new TH2F("h2_momentum_nphe_negative_pions_HTCC", "Momentum vs nphe for negative pions HTCC ;Momentum (GeV);nphe(counts)", 100, 0, 10, 100, 0, 100);
    TH2F* h2_momentum_nphe_negative_kaons_HTCC = new TH2F("h2_momentum_nphe_negative_kaons_HTCC", "Momentum vs nphe for negative kaons HTCC ;Momentum (GeV);nphe(counts)", 100, 0, 10, 100, 0, 100);

    // Plot of momentum for charged hadrons
    TH1F* h1_momentum_positive = new TH1F("h1_momentum_positive", "Momentum for positive particles ; Momentum (GeV); counts", 100, 0, 10); 
    TH1F* h1_momentum_negative = new TH1F("h1_momentum_negative", "Momentum for negative particles; Momentum (GeV); counts", 100, 0, 10); 

    
    TH1F* h1_nphe_positive = new TH1F ("h1_nphe_positive", "nphe for positive charge particles;nphe;counts", 100, 0, 100);
    TH1F* h1_nphe_negative = new TH1F ("h1_nphe_negative", "nphe for negative charge particles;nphe;counts", 100, 0, 100);

    TH1F* h1_chi2pid_proton = new TH1F ("h1_chi2pid_proton", "chi2pid for protons; chi2pid;counts", 100, -15, 15); 
    TH1F* h1_chi2pid_pion_positive = new TH1F ("h1_chi2pid_pion_positive", "chi2pid for pion positive; chi2pid;counts", 100, -15, 15); 
    TH1F* h1_chi2pid_pion_negative = new TH1F ("h1_chi2pid_pion_negative", "chi2pid for pion negative; chi2pid;counts", 100, -15, 15);
    TH1F* h1_chi2pid_kaon_positive = new TH1F ("h1_chi2pid_kaon_positive", "chi2pid for kaon positive; chi2pid;counts", 100, -15, 15); 
    TH1F* h1_chi2pid_kaon_negative = new TH1F ("h1_chi2pid_kaon_negative", "chi2pid for kaon negative; chi2pid;counts", 100, -15, 15); 


    TH1F* h1_detector = new TH1F("h1_detector", "detector; detector; counts", 100, 0, 25);
    TH1F* h1_pid = new TH1F("h1_pid", "pid; pid; counts", 100, -5000, 5000);

    TH2F* h2_momentum_beta_positive = new TH2F("h2_momentum_beta_positive", "Momentum vs beta for positive particles; Momentum (GeV); beta", 100, 0, 10, 100, -5, 5);
    TH2F* h2_momentum_beta_positive_pion = new TH2F("h2_momentum_beta_positive_pion", "Momentum vs beta for positive pions; Momentum (GeV); beta", 100, 0, 10, 100, -5, 5);
    TH2F* h2_momentum_beta_positive_kaon = new TH2F("h2_momentum_beta_positive_kaon", "Momentum vs beta for positive kaons; Momentum (GeV); beta", 100, 0, 10, 100, -5, 5);
    TH2F* h2_momentum_beta_proton = new TH2F("h2_momentum_beta_proton", "Momentum vs beta for protons; Momentum (GeV); beta", 100, 0, 10, 100, -5, 5);
    TH2F* h2_momentum_beta_negative = new TH2F("h2_momentum_beta_negative", "Momentum vs beta for negative particles; Momentum (GeV); beta", 100, 0, 10, 100, -5, 5);
    TH2F* h2_momentum_beta_negative_pion = new TH2F("h2_momentum_beta_negative_pion", "Momentum vs beta for negative pions; Momentum (GeV); beta", 100, 0, 10, 100, -5, 5);
    TH2F* h2_momentum_beta_negative_kaon = new TH2F("h2_momentum_beta_negative_kaon", "Momentum vs beta for negative kaons; Momentum (GeV); beta", 100, 0, 10, 100, -5, 5);







    //////////////////////////////////////////////////////////////////////
    //TH1F
    //////////////////////////////////////////////////////////////////////
    

    // Iterate over directories
    for (const auto& dir : directories) {
        cout << "Processing directory: " << dir << endl;

        // Find all .hipo files in the directory
        vector<string> hipoFiles;
        for (const auto& entry : fs::directory_iterator(dir)) {
            if (entry.is_regular_file() && entry.path().extension() == ".hipo") {
                hipoFiles.push_back(entry.path().string());
            }
        }

        if (hipoFiles.empty()) {
            cout << "  No .hipo files found in directory: " << dir << endl;
            continue;
        }

        // Iterate over files in the directory
        for (const auto& file : hipoFiles) {
            cout << "  Opening file: " << file << endl;

            // Open HIPO file
            hipo::reader reader;
            reader.open(file.c_str());

            hipo::dictionary factory;
            reader.readDictionary(factory);

            hipo::event event;
            int counter = 0;
            int electron_count = 0; 
            int selected_events = 0; 

            // Load banks
            hipo::bank PART(factory.getSchema("REC::Particle"));
            hipo::bank EVENT(factory.getSchema("REC::Event"));
            hipo::bank SCIN(factory.getSchema("REC::Scintillator"));
            hipo::bank CHER(factory.getSchema("REC::Cherenkov"));
            //const int maxEvents = 10000000; 
            //const int maxEvents = 1000000;
	        const int maxEvents = 10000; 
            // Process events in the file
            while (reader.next() && counter < maxEvents) {
                reader.read(event);
                event.getStructure(PART);
                event.getStructure(EVENT);
                event.getStructure(SCIN);
                event.getStructure(CHER); 

                // Map SCIN and Particle banks using `pindex`
                IndexMap scinMap = loadMapByIndex(SCIN, "pindex");
                IndexMap cherMap = loadMapByIndex(CHER,"pindex");

               // int nrows_event = EVENT.getRows();
                float startTime = EVENT.getFloat("startTime", 0);
                if (startTime < 0) continue;

                int nrows_part = PART.getRows();
                int trigger_electron_index = -1;

                // Step 1: Select trigger electron
                bool has_trigger_electron = false;
                for (int i = 0; i < PART.getRows(); i++) {
                    // Get particle information
                     pid = PART.getInt("pid", i);       // Particle ID
                     float status = PART.getInt("status", i); // Status
                     chi2pid = PART.getFloat("chi2pid", i); // Chi2 PID
                    
                    // Apply trigger electron cuts
                    if (pid == 11 && status < 0 && chi2pid > -5 && chi2pid < 5) {
                        has_trigger_electron = true;
                        trigger_electron_index = i; // Store the index of the trigger electron
                        break;     
                    }
                }

                // Count events with a valid trigger electron
                if (has_trigger_electron) {
                    selected_events++;

                    // Step 2: Loop over particles for pions, kaons, and protons
                    for (int i = 0; i < nrows_part; i++) {
                        if (i == trigger_electron_index) continue; // Skip the trigger electron

                         charge = PART.getByte("charge", i);
                        if (charge == 0) continue; // **Charge cut**: Skip neutral or negatively charged particles
                        int status = PART.getShort("status",i); 
                        if (abs(status/2000) != 1) continue;

                         pid = PART.getInt("pid", i);
                         chi2pid = PART.getFloat("chi2pid", i);
                         beta = PART.getFloat("beta", i); 
                        float px = PART.getFloat("px", i);
                        float py = PART.getFloat("py", i);
                        float pz = PART.getFloat("pz", i);
                        float p = sqrt(px * px + py * py + pz * pz);
                        h1_pid->Fill(pid);
                        

                        if (charge > 0){
                            
                            h1_momentum_positive->Fill(p);
                            h2_momentum_beta_positive->Fill(p,beta);
                            h2_momentum_chi2pid_positive->Fill(p,chi2pid);

                            if (pid == +2212 /*proton*/){
                                h1_chi2pid_proton->Fill(chi2pid);
                                h2_momentum_beta_proton->Fill(p,chi2pid);
                            }
                            if (pid == 211 /* positive pion*/){
                                h1_chi2pid_pion_positive->Fill(chi2pid);
                                h2_momentum_beta_positive_pion->Fill(p,beta);
                            }

                            if (pid == 321 /* positive kaon*/){
                                h1_chi2pid_kaon_positive->Fill(chi2pid);
                                h2_momentum_beta_positive_kaon->Fill(p,beta);
                            }
                        } else if(charge < 0){
                            h1_momentum_negative->Fill(p);
                            h2_momentum_beta_negative->Fill(p,beta);
                            h2_momentum_chi2pid_negative->Fill(p,chi2pid);
                            if (pid == -211 /* positive pion*/){
                                h1_chi2pid_pion_negative->Fill(chi2pid);
                                h2_momentum_beta_negative_pion->Fill(p,beta);
                            }

                            if (pid == -321){
                                h1_chi2pid_kaon_negative->Fill(chi2pid);
                                h2_momentum_beta_negative_kaon->Fill(p,beta);
                            }
                        }


                        // Clear vectors for Cherenkov hits at the beginning of each particle loop
                        nphe.clear();
                        detector.clear();

                        momentum = p;
                       
                        
                        


                        ///////////////////////////////////////////////////////////////////////
                        // For Cherenkov Signal REC::Cherenkov 
                        ///////////////////////////////////////////////////////////////////////

                        if (cherMap.find(i) != cherMap.end()){

                            // Process associated CHER rows for this particle
                            for (int iCherRow : cherMap[i]){
                                int det = CHER.getByte("detector", iCherRow);
                                int nphe_val = CHER.getFloat("nphe", iCherRow); 
                                detector.push_back(det);
                                nphe.push_back(static_cast<int>(nphe_val));
                                h1_detector->Fill(det);
                                if (charge > 0) {
                                    h1_nphe_positive->Fill(nphe_val);
                                    h2_momentum_nphe_positive->Fill(p,nphe_val);
                                    if (det == 15 /* HTCC*/){
                                        h2_momentum_nphe_positive_HTCC->Fill(p,nphe_val);
                                        if (pid == 211){
                                            h2_momentum_nphe_positive_pions_HTCC->Fill(p,nphe_val);

                                        }
                                        if (pid == 321){
                                            h2_momentum_nphe_positive_kaons_HTCC->Fill(p,nphe_val);

                                        }
                                    }  
                                    if (det == 16 /*LTCC*/){
                                        h2_momentum_nphe_positive_LTCC->Fill(p,nphe_val);
                                        if (pid ==221) {
                                            h2_momentum_nphe_positive_pions_LTCC->Fill(p,nphe_val);
                                        }
                                        if (pid == 321){
                                            h2_momentum_nphe_positive_kaons_LTCC->Fill(p,nphe_val);
                                        }
                                    }
                                }

                                if (charge < 0) {
                                    h1_nphe_negative->Fill(nphe_val);
                                    h2_momentum_nphe_negative->Fill(p,nphe_val);
                                    if (det == 15){
                                        h2_momentum_nphe_negative_HTCC->Fill(p,nphe_val);
                                        if (pid == -211){
                                            h2_momentum_nphe_negative_pions_HTCC->Fill(p,nphe_val);
                                        }
                                        if (pid = -321){
                                        
                                            h2_momentum_nphe_negative_kaons_HTCC->Fill(p,nphe_val);
                                        }
                                    }
                                    if (det == 16){
                                        h2_momentum_nphe_negative_LTCC->Fill(p,nphe_val);
                                        if (pid == -211 ){
                                            h2_momentum_nphe_negative_pions_LTCC->Fill(p,nphe_val);
                                        }
                                        if (pid == -321){
                                            h2_momentum_nphe_negative_kaons_LTCC->Fill(p,nphe_val);
                                        }
                                    }
                                }
                            }

                            /* if (pid != 211 && pid != 321 && pid != 2212) continue; // Process only pions, kaons, protons
                            // Process associated SCIN rows for this particle
                            if (scinMap.find(i) != scinMap.end()) {
                                for (int iScinRow : scinMap[i]) {
                                    int detector = SCIN.getByte("detector", iScinRow);
                                    if (detector == 12) { // Match to desired detector
                                        float path = SCIN.getFloat("path", iScinRow);
                                        float time = SCIN.getFloat("time", iScinRow);

                                        float px = PART.getFloat("px", i);
                                        float py = PART.getFloat("py", i);
                                        float pz = PART.getFloat("pz", i);
                                        float p = sqrt(px * px + py * py + pz * pz);

                                        // Step 3: Fill trees and histograms for respective particles
                                        if (pid == 211) {
                                            beta_pion = calculateBeta(p, mass_pion);
                                            delta_t_pion = calculateDeltaT(path, time, beta_pion, startTime);
                                            momentum_pion = p;
                                            tree_pion->Fill();
                                            h2_momentum_deltaT_pion->Fill(p, delta_t_pion);
                                        } else if (pid == 321) {
                                            beta_kaon = calculateBeta(p, mass_kaon);
                                            delta_t_kaon = calculateDeltaT(path, time, beta_kaon, startTime);
                                            momentum_kaon = p;
                                            tree_kaon->Fill();
                                            h2_momentum_deltaT_kaon->Fill(p, delta_t_kaon);
                                        } else if (pid == 2212) {
                                            beta_proton = calculateBeta(p, mass_proton);
                                            delta_t_proton = calculateDeltaT(path, time, beta_proton, startTime);
                                            momentum_proton = p;
                                            tree_proton->Fill();
                                            h2_momentum_deltaT_proton->Fill(p, delta_t_proton);
                                        }
                                    }
                                }
                            } */
                        }
                        tree_pparticles->Fill();
                    }  
                } 
                counter++;
            }

        }
    }

    // Save the trees to a ROOT file
    TFile *outputFile = new TFile("particles_data_chi2pid_beta_10M.root", "RECREATE");
    tree_pparticles->Write();
    /* tree_kaon->Write();
    tree_proton->Write(); */
    /* h2_momentum_deltaT_pion->Write();  // Write the pion histogram to the file
    h2_momentum_deltaT_kaon->Write();  // Write the kaon histogram to the file
    h2_momentum_deltaT_proton->Write(); 
    */
   h1_chi2pid_kaon_negative->Write();
   h1_chi2pid_kaon_positive->Write();
   h1_chi2pid_pion_negative->Write();
   h1_chi2pid_pion_positive->Write();
   h1_chi2pid_proton->Write();
   h1_detector->Write();
   h1_momentum_negative->Write();
   h1_momentum_positive->Write();
   h1_nphe_negative->Write();
   h1_nphe_positive->Write();

   h2_momentum_beta_negative->Write();
   h2_momentum_beta_negative_kaon->Write();
   h2_momentum_beta_negative_pion->Write();
   h2_momentum_beta_positive->Write();
   h2_momentum_beta_positive_kaon->Write();
   h2_momentum_beta_positive_pion->Write();
   h2_momentum_beta_proton->Write();

   h2_momentum_chi2pid_negative->Write();
   h2_momentum_chi2pid_positive->Write();


   h2_momentum_nphe_negative->Write();
   h2_momentum_nphe_negative_HTCC->Write();
   h2_momentum_nphe_negative_LTCC->Write();


   h2_momentum_nphe_positive->Write();
   h2_momentum_nphe_positive_HTCC->Write();
   h2_momentum_nphe_positive_LTCC->Write();

   h2_momentum_nphe_negative_kaons_HTCC->Write();
   h2_momentum_nphe_negative_kaons_LTCC->Write();
   h2_momentum_nphe_negative_pions_HTCC->Write();
   h2_momentum_nphe_negative_pions_LTCC->Write();
   h2_momentum_nphe_positive_kaons_HTCC->Write();
   h2_momentum_nphe_positive_kaons_LTCC->Write();
   h2_momentum_nphe_positive_pions_HTCC->Write();
   h2_momentum_nphe_positive_pions_LTCC->Write();

   h1_pid->Write(); 

  /* TCanvas *canvas = new TCanvas("canvas_name", "Canvas Title", 800, 600);
     h1_chi2pid_pion_negative->Draw();  // Draw the histogram

    // Save the canvas as a PNG file
    canvas->SaveAs("h1_chi2pid_pion_negative.png");
   
 */


    // Close the ROOT file
    outputFile->Close();

    timer.Stop(); // Stop the timer

    std::cout << "Real time: " << timer.RealTime() << " s" << std::endl;
    std::cout << "CPU time: " << timer.CpuTime() << " s" << std::endl;

    return 0;
}
