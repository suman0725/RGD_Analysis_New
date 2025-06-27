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
float calculateBeta(float p, float mass) {
    float E = sqrt(p * p + mass * mass);
    return p / E;
}

// Function to calculate delta_t for a given particle
float calculateDeltaT(float path, float time, float beta, float startTime) {
    return time - (path / (beta * 29.9792458)) - startTime;
}

int main() {
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
    const float mass_pion = 0.13957;
    const float mass_kaon = 0.49367;
    const float mass_proton = 0.93827;

    // Create ROOT trees for storing delta t and momentum for different particles
    TTree* tree_pion = new TTree("pion", "Pion+ Tree");
    TTree* tree_kaon = new TTree("kaon", "Kaon+ Tree");
    TTree* tree_proton = new TTree("proton", "Proton Tree");

    // Variables to hold delta t and momentum for each particle type
    float delta_t_pion, momentum_pion, beta_pion;
    float delta_t_kaon, momentum_kaon, beta_kaon;
    float delta_t_proton, momentum_proton, beta_proton;

    // Create branches for the trees
    tree_pion->Branch("delta_t", &delta_t_pion, "delta_t/F");
    tree_pion->Branch("momentum", &momentum_pion, "momentum/F");
    //tree_pion->Branch("beta", &beta_pion, "beta/F");

    tree_kaon->Branch("delta_t", &delta_t_kaon, "delta_t/F");
    tree_kaon->Branch("momentum", &momentum_kaon, "momentum/F");
    //tree_kaon->Branch("beta", &beta_kaon, "beta/F");

    tree_proton->Branch("delta_t", &delta_t_proton, "delta_t/F");
    tree_proton->Branch("momentum", &momentum_proton, "momentum/F");
    //tree_proton->Branch("beta", &beta_proton, "beta/F");

    TH2F* h2_momentum_deltaT_pion = new TH2F("h2_momentum_deltaT_pion", "Pion Momentum vs DeltaT;Momentum (GeV/c);Delta T (ns)", 100, 0, 10, 100, -10, 10);
    TH2F* h2_momentum_deltaT_kaon = new TH2F("h2_momentum_deltaT_kaon", "Kaon Momentum vs DeltaT;Momentum (GeV/c);Delta T (ns)", 100, 0, 10, 100, -10, 10);
    TH2F* h2_momentum_deltaT_proton = new TH2F("h2_momentum_deltaT_proton", "Proton Momentum vs DeltaT;Momentum (GeV/c);Delta T (ns)", 100, 0, 10, 100, -10, 10);

    

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
            int selected_events = 0;

            // Load banks
            hipo::bank PART(factory.getSchema("REC::Particle"));
            hipo::bank EVENT(factory.getSchema("REC::Event"));
            hipo::bank SCIN(factory.getSchema("REC::Scintillator"));
             const int maxEvents = 10000000; 
	        // const int maxEvents = 1; 


            // Process events in the file
            while (reader.next() && counter < maxEvents ) {

                reader.read(event);
                event.getStructure(PART);
                event.getStructure(EVENT);
                event.getStructure(SCIN);
                //SCIN.show(); 

                // Load the map between SCIN and Particle banks using the `pindex` variable
                IndexMap scinMap = loadMapByIndex(SCIN, "pindex");


                bool has_trigger_electron = false;
                for (int i = 0; i < PART.getRows(); i++) {
                    // Get particle information
                    int pid = PART.getInt("pid", i);       // Particle ID
                    int status = PART.getInt("status", i); // Status
                    float chi2pid = PART.getFloat("chi2pid", i); // Chi2 PID
                    
                    // Apply trigger electron cuts
                    if (pid == 11 && status < 0 && chi2pid > -5 && chi2pid < 5) {
                        has_trigger_electron = true;
                       
                    }
                }

                // Count events with a valid trigger electron
                if (has_trigger_electron) {
                    selected_events++;
                
    
                    //int nerows = EVENT.getRows();
                    float startTime = EVENT.getFloat("startTime", 0);
                    if (startTime < 0) continue; 
                    
                    int nrows = PART.getRows();
                    for (int i = 0; i < nrows; i++) {

                        int charge = PART.getInt("charge", i);
                        // Skip neutral or negatively charged particles
                        if (charge <= 0) continue;

                        // Get particle properties
                        float px = PART.getFloat("px", i);
                        float py = PART.getFloat("py", i);
                        float pz = PART.getFloat("pz", i);
                        float p = sqrt(px*px + py*py + pz*pz);

                        // Calculate beta for each particle based on momentum and mass
                        beta_pion = calculateBeta(p, mass_pion);
                        beta_kaon = calculateBeta(p, mass_kaon);
                        beta_proton = calculateBeta(p, mass_proton);

                        // Process associated SCIN rows for this particle
                        if (scinMap.find(i) != scinMap.end()) {
                            for (int iScinRow : scinMap[i]) {
                                int detector = SCIN.getByte("detector",iScinRow);
                                //cout << "detector = " << detector << endl;
                                if ( detector == 12) {
                                    // Get path and time for beta calculation
                                    float path = SCIN.getFloat("path", iScinRow);
                                    float time = SCIN.getFloat("time", iScinRow);
                                    
                                
                                    // Calculate delta_t for each particle type using the calculated beta
                                    delta_t_pion = calculateDeltaT(path, time, beta_pion, startTime);
                                    delta_t_kaon = calculateDeltaT(path, time, beta_kaon, startTime);
                                    delta_t_proton = calculateDeltaT(path, time, beta_proton, startTime);

                                    // Fill trees with delta_t, beta, and momentum values for each particle type
                                    momentum_pion = p;
                                    tree_pion->Fill();

                                    momentum_kaon = p;
                                    tree_kaon->Fill();

                                    momentum_proton = p;
                                    tree_proton->Fill();

                                    
                                    // Inside the event loop, fill the histograms for each particle type
                                    h2_momentum_deltaT_pion->Fill(p, delta_t_pion);  // For pion
                                    h2_momentum_deltaT_kaon->Fill(p, delta_t_kaon);  // For kaon
                                    h2_momentum_deltaT_proton->Fill(p, delta_t_proton);  // For proton
                                }
                            }
                        } 
                    }
                    counter++;
                }
            }
        }
    }

    // Save the trees to a ROOT file
    TFile *outputFile = new TFile("particles_data.root", "RECREATE");
    tree_pion->Write();
    tree_kaon->Write();
    tree_proton->Write();
    h2_momentum_deltaT_pion->Write();  // Write the pion histogram to the file
    h2_momentum_deltaT_kaon->Write();  // Write the kaon histogram to the file
    h2_momentum_deltaT_proton->Write(); 

    // Close the ROOT file
    outputFile->Close();

    return 0;
}
