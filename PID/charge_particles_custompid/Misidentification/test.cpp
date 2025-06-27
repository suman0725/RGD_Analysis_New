#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include "reader.h"
#include "ParticleData.h"
#include "TFile.h"
#include "TTree.h"
#include <TH2F.h>
#include <TCanvas.h>

using namespace std;
namespace fs = std::filesystem;

int main() {
    // Load CCDB parameters
    loadCCDBParams();

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
        cerr << "Error: No directories found in directories.txt" << endl;
        return 1;
    }

    // Loop over directories and HIPO files
    int event_count = 0;
    int events_with_trigger = 0;
    int maxEvents = 10000; // Limit to 1M events
    int electron_count = 0;
    int positron_count = 0;
    int trigger_electron_count = 0; 
    for (const auto& dir : directories) {
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
            cout << "Processing file: " << file << endl;
            // Open HIPO file
            hipo::reader reader;
            reader.open(file.c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);

            // Loop over events
            while (reader.next() && event_count < maxEvents) {
                event_count++;
                hipo::event event;
                reader.read(event);

                // Load necessary banks
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

                // Create index maps
                IndexMap cherMap = loadMapByIndex(CHER, "pindex");
                IndexMap caloMap = loadMapByIndex(CALO, "pindex");
                IndexMap scinMap = loadMapByIndex(SCIN, "pindex");

                // Extract particles
                vector<ParticleData> particles(PART.getRows());
                for (int i = 0; i < PART.getRows(); ++i) {
                    particles[i] = getParticleData(i, PART, CHER, CALO, SCIN, cherMap, caloMap, scinMap);
                }
////////////////////////////////////////////////////////////////////////////////////////////////////
                // Count electrons, positrons, and trigger electrons with additional cuts
                bool hadTrigger = false;
                for (const auto& p : particles) {
                    if (isSimpleElectron(p)) {
                        if (p.charge == -1) electron_count++;
                        else if (p.charge == 1) positron_count++;
                    }
                    if (isTriggerElectron(p)) {
                        float sf = p.energy_total / p.p;
                        float mean = getSamplingFractionMean(p.sector, p.energy_total);
                        float sigma = getSamplingFractionSigma(p.sector, p.energy_total);
                        float chi2pid = getSamplingFractionNSigma(sf, mean, sigma);
                        if (abs(chi2pid) < 5 && p.vz >= -10.55 && p.vz <= 5 && abs(p.status) / 1000 == 2) {//p.vz >= -20 && p.vz <= 5 for LD2 electron
                            trigger_electron_count++;
                            hadTrigger = true;
                        }
                    }
                } 
////////////////////////////////////////////////////////////////////////////////////////////////////
                /* for (const auto& p : particles) {
                        if (p.pid == 11) {electron_count++;
                            if (abs(p.chi2pid) < 5 && p.vz >= -20 && p.vz <= 5 && abs(p.status) / 1000 == 2 && p.status <0 ) {
                            trigger_electron_count++;
                            hadTrigger = true;
                        }
                        }
                        if (p.pid == -11)  positron_count++;
                        

                } */
////////////////////////////////////////////////////////////////////////////////////////////////////
        
////////////////////////////////////////////////////////////////////////////////////////////////////
                if (hadTrigger) {

                    events_with_trigger++;
                    for (const auto& p : particles) {
                        if (abs(p.status) / 2000 != 1) continue; 
                        if (isSimpleElectron(p)) continue; // Skip electrons/positrons
                        if (p.charge <= 0) continue; // Only positive hadrons
                    
                        //float beta = 
                        if ( p.pid== 321 ) { // Positive pions

                        PART.show(); 
                            float beta_theory = getTheoryBeta(p, PION_MASS); 
                           // cout << "beta_theory for pion" << beta_theory << endl;
                            cout << "beta from rec::particle bank for pion: " << p.beta << endl; 
                            float mass_pion = getCalculatedMass(p); 
                            cout << "reconsructed mass for pion: " << mass_pion << endl; 
                        }
                        
                    }
                }
            }
            if (event_count >= maxEvents) break;
        }
        if (event_count >= maxEvents) break;
    }

    cout << "Total events processed: " << event_count << endl;
    cout << "Total electron count: " << electron_count << endl;
    cout << "Total positron count: " << positron_count << endl;
    cout << "Total trigger electron count: " << trigger_electron_count << endl;
    cout << "Events with trigger electrons: " << events_with_trigger << endl; 

    

    return 0;
}