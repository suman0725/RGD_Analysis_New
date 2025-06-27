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

    // Create ROOT file and trees with updated names
    //TFile* outputFile = new TFile("piontree.root", "RECREATE");
    TFile* outputFile = new TFile("pkptreeCxC_5.root", "RECREATE");
    TTree* treeEBPidPions = new TTree("EB_pid_pions", "Trigger particles with pid == 211");
    TTree* treeEBPidKaons = new TTree("EB_pid_kaons", "Trigger particles with pid == 321");
    TTree* treeEBPidProtons = new TTree("EB_pid_protons", "Trigger particles with pid == 2212");
    TTree* treePionHypothesis = new TTree("pion_hypothesis", "Trigger particles with pion hypothesis");
    TTree* treeKaonHypothesis = new TTree("kaon_hypothesis", "Trigger particles with kaon hypothesis");
    TTree* treeProtonHypothesis = new TTree("proton_hypothesis", "Trigger particles with proton hypothesis");
    

    // Variables to store in the trees
    float chi2pid, momentum, mass, beta;
    int pid; 
    float chi2pid_kaon = 99999.0f;
    float chi2pid_proton = 99999.0f;

    treeEBPidPions->Branch("chi2pid", &chi2pid, "chi2pid/F");
    //treeEBPidPions->Branch("mass", &mass, "mass/F"); 
    treeEBPidPions->Branch("p", &momentum, "p/F");
    treeEBPidPions->Branch("beta", &beta, "beta/F");

    treePionHypothesis->Branch("chi2pid", &chi2pid, "chi2pid/F");
    treePionHypothesis->Branch("p", &momentum, "p/F");
    //treePionHypothesis->Branch("mass", &mass, "mass/F"); 
    //treePionHypothesis->Branch("pid", &pid, "pid/I");

    // EB Kaons (PID == 321) / Kaons hypotheis (Comparing all par)
    treeEBPidKaons->Branch("chi2pid", &chi2pid, "chi2pid/F");
    treeEBPidKaons->Branch("p", &momentum, "p/F");
    treeEBPidKaons->Branch("beta", &beta, "beta/F");
    treeKaonHypothesis->Branch("chi2pid_kaon", &chi2pid_kaon, "chi2pid_kaon/F");
    treeKaonHypothesis->Branch("p", &momentum, "p/F");
    

    treeEBPidProtons->Branch("chi2pid", &chi2pid, "chi2pid/F");
    treeEBPidProtons->Branch("p", &momentum, "p/F");
    treeEBPidProtons->Branch("beta", &beta, "beta/F");
    treeProtonHypothesis->Branch("chi2pid_proton", &chi2pid_proton, "chi2pid_proton/F");
    treeProtonHypothesis->Branch("p", &momentum, "p/F");

    
    
    
    

    // 2D Histograms for chi2pid vs momentum
    //TH2F* h2_truePions = new TH2F("h2_truePions", "EB Pions (pid == 211);p [GeV];chi2pid", 200, 0, 10, 200, -15, 15); // Match main.cpp binning
    //TH2F* h2_pionHypothesis = new TH2F("h2_pionHypothesis", "Pion Hypothesis; p [GeV]; chi2pid", 200, 0, 10, 200, -15, 15); // Match main.cpp binning

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
    //int maxEvents = 100000;
    int maxEvents = 100000000; // Limit to 1M events
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
                hipo::bank CON(factory.getSchema("RUN::config")); 
                event.getStructure(PART);
                event.getStructure(EVENT);
                event.getStructure(SCIN);
                event.getStructure(CHER);
                event.getStructure(CALO);
                event.getStructure(CON);

                // Create index maps
                IndexMap cherMap = loadMapByIndex(CHER, "pindex");
                IndexMap caloMap = loadMapByIndex(CALO, "pindex");
                IndexMap scinMap = loadMapByIndex(SCIN, "pindex");

                // Extract particles
                vector<ParticleData> particles(PART.getRows());
                for (int i = 0; i < PART.getRows(); ++i) {
                    particles[i] = getParticleData(i, PART, CHER, CALO, SCIN, cherMap, caloMap, scinMap);
                }

                float torus = CON.getFloat("torus", 0); 
                if (torus != 1.0) continue;
////////////////////////////////////////////////////////////////////////////////////////////////////
                //manual way to select electrons/positrons
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
                //EB pid for electrons/positrons
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
                        momentum = p.p;
                        pid = p.pid;
                        beta = p.beta; 
                        mass = getCalculatedMass(p); 
                        chi2pid = p.chi2pid; 
                        if (p.pid == 211) { // Positive pions    
                            treeEBPidPions->Fill();
                            //h2_truePions->Fill(momentum, chi2pid);
                        }
                        if (p.pid == 321) { // Positive pions    
                            treeEBPidKaons->Fill();
                            //h2_trueKaons->Fill(momentum, chi2pid);
                        }
                        if (p.pid == 2212) { // Positive pions    
                            treeEBPidProtons->Fill();
                            //h2_trueKaons->Fill(momentum, chi2pid);
                        }
                        // Compute chi2pid for all positive hadrons (including pid == 211)
                        chi2pid = computeChi2pid(p, PION_MASS);
                        chi2pid_kaon = computeChi2pid(p, KAON_MASS); 
                        chi2pid_proton = computeChi2pid(p, PROTON_MASS); 
                        if (chi2pid != 99999.0f) { // Only fill if chi2pid is valid
                            treePionHypothesis->Fill();
                            treeKaonHypothesis->Fill();
                            treeProtonHypothesis->Fill(); 
                            //h2_pionHypothesis->Fill(momentum, chi2pid);
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

    /* // Create canvases
    TCanvas* c1 = new TCanvas("c1", "EB PID 211 Pions", 800, 600);
    h2_truePions->SetStats(0); // Match main.cpp
    h2_truePions->Draw("COLZ");
    gPad->SetLogz(); // Match main.cpp
    c1->SaveAs("chi2pid_vs_p_truePions.png");

    TCanvas* c2 = new TCanvas("c2", "Pion Hypothesis", 800, 600);
    h2_pionHypothesis->SetStats(0); // Match main.cpp
    h2_pionHypothesis->Draw("COLZ");
    gPad->SetLogz(); // Match main.cpp
    c2->SaveAs("chi2pid_vs_p_pionHypothesis.png"); */

    outputFile->Write();
    outputFile->Close();
    delete outputFile;

    return 0;
}