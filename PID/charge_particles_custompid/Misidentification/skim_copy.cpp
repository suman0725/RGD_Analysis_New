#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <chrono>
#include <iterator>
#include "reader.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;
using namespace std::chrono;
namespace fs = std::filesystem;

map<int, float> PDG_MASS = {
    {11, 0.0005}, {-11, 0.0005}, {211, 0.13957018}, {-211, 0.13957018},
    {321, 0.49367716}, {-321, 0.49367716}, {2212, 0.938272046}, {-2212, 0.938272046},
    {45, 1.87705}
};

const vector<pair<int, vector<int>>> chargedBetaDetectors = {
    {FTOF_DETECTOR, {2, 1, 3}}, {ECAL_DETECTOR, {1, 4, 7}}
};

// Define ParticleData struct
struct ParticleData {
    int pid;
    float p;
    float beta;
    float chi2pid;
    float vz;
    int status;
};

// Define getParticleData function
ParticleData getParticleData(int index, hipo::bank& bank) {
    ParticleData pd;
    pd.pid = bank.getInt("pid", index);
    pd.p = sqrt(pow(bank.getFloat("px", index), 2) +
                pow(bank.getFloat("py", index), 2) +
                pow(bank.getFloat("pz", index), 2));
    pd.beta = bank.getFloat("beta", index);
    pd.chi2pid = bank.getFloat("chi2pid", index);
    pd.vz = bank.getFloat("vz", index);
    pd.status = bank.getInt("status", index);
    return pd;
}

ParticleData getParticleData(int partIdx, hipo::bank& PART, hipo::bank& CHER, hipo::bank& CALO, hipo::bank& SCIN,
    IndexMap& cherMap, IndexMap& caloMap, IndexMap& scinMap) {
    ParticleData pd;

    pd.charge = PART.getByte("charge", partIdx);
    pd.pid = PART.getInt("pid", partIdx);
    float px = PART.getFloat("px", partIdx);
    float py = PART.getFloat("py", partIdx);
    float pz = PART.getFloat("pz", partIdx);
    pd.p = std::sqrt(px * px + py * py + pz * pz);
    pd.vt = PART.getFloat("vt", partIdx);
    pd.beta = PART.getFloat("beta", partIdx);
    pd.status = PART.getShort("status", partIdx);
    pd.chi2pid = PART.getFloat("chi2pid", partIdx);
    pd.start_time = pd.vt;

    if (caloMap.find(partIdx) != caloMap.end()) {
        for (int iCalRow : caloMap[partIdx]) {
            if (CALO.getByte("detector", iCalRow) == ECAL_DETECTOR) {
                pd.sector = CALO.getByte("sector", iCalRow);
                int layer = CALO.getByte("layer", iCalRow);
                float energy = CALO.getFloat("energy", iCalRow);
                float time = CALO.getFloat("time", iCalRow);
                float path = CALO.getFloat("path", iCalRow);
                if (layer == 1) pd.energy_pcal = energy;
                pd.energy_total += energy;
                pd.hits.push_back({ECAL_DETECTOR, layer, pd.sector, 0, time, path});
            }
        }
    }

    if (cherMap.find(partIdx) != cherMap.end()) {
        for (int iCherRow : cherMap[partIdx]) {
            int detector = CHER.getByte("detector", iCherRow);
            float nphe = CHER.getFloat("nphe", iCherRow);
            if (detector == HTCC_DETECTOR) pd.nphe_htcc = nphe;
            else if (detector == LTCC_DETECTOR) pd.nphe_ltcc = nphe;
        }
    }

    if (scinMap.find(partIdx) != scinMap.end()) {
        for (int iScinRow : scinMap[partIdx]) {
            int detector = SCIN.getByte("detector", iScinRow);
            int layer = SCIN.getByte("layer", iScinRow);
            int sector = SCIN.getByte("sector", iScinRow);
            int component = SCIN.getShort("component", iScinRow);
            float time = SCIN.getFloat("time", iScinRow);
            float path = SCIN.getFloat("path", iScinRow);
            pd.hits.push_back({detector, layer, sector, component, time, path});
        }
    }

    return pd;
}

float getTheoryBeta(float p, float mass) {
    return p / sqrt(p * p + mass * mass);
}


float computeDeltaT(const ParticleData& p, int pid) {
    float delta_t = 99999.0f;
    bool found = false;

    for (const auto& det : chargedBetaDetectors) {
        for (int layer : det.second) {
            if (hasHit(p, det.first, layer)) {
                for (const auto& hit : p.hits) {
                    if (hit.detector == det.first && hit.layer == layer) {
                        float mass = PDG_MASS[pid];
                        float beta_theory = p.p / sqrt(p.p * p.p + mass * mass);
                        if (beta_theory <= 0) return 99999.0f;
                        float vt = hit.time - hit.path / (SPEED_OF_LIGHT * beta_theory);
                        delta_t = vt - p.start_time;
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
            if (found) break;
        }
        if (found) break;
    }
    return delta_t;
}

int main() {
    // Start timing
    auto startTime = high_resolution_clock::now();

    // File to save processed HIPO files
    ofstream processedFilesStream("processed_hipo_files.txt");
    if (!processedFilesStream.is_open()) {
        cerr << "Error: Could not open processed_hipo_files.txt" << endl;
        return 1;
    }

    // Create ROOT file and trees with updated names
    TFile* outputFile = new TFile("pkptreeCxC_7.root", "RECREATE");
    TTree* treeEBPidPions = new TTree("EB_pid_pions", "Trigger particles with pid == 211");
    TTree* treeEBPidKaons = new TTree("EB_pid_kaons", "Trigger particles with pid == 321");
    TTree* treeEBPidProtons = new TTree("EB_pid_protons", "Trigger particles with pid == 2212");

    // Variables to store in the trees
    float chi2pid, momentum, beta, vz;
    bool hadTrigger = false;

    treeEBPidPions->Branch("chi2pid", &chi2pid, "chi2pid/F");
    treeEBPidPions->Branch("p", &momentum, "p/F");
    treeEBPidPions->Branch("beta", &beta, "beta/F");
    treeEBPidPions->Branch("vz", &vz, "vz/F");


    treeEBPidKaons->Branch("chi2pid", &chi2pid, "chi2pid/F");
    treeEBPidKaons->Branch("p", &momentum, "p/F");
    treeEBPidKaons->Branch("beta", &beta, "beta/F");
    treeEBPidKaons->Branch("vz", &vz, "vz/F");

    treeEBPidProtons->Branch("chi2pid", &chi2pid, "chi2pid/F");
    treeEBPidProtons->Branch("p", &momentum, "p/F");
    treeEBPidProtons->Branch("beta", &beta, "beta/F");
    treeEBPidProtons->Branch("vz", &vz, "vz/F");


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

    if (directories.empty()) {
        cerr << "Error: No directories found in directories.txt" << endl;
       
        return 1;
    }

    // Loop over directories and HIPO files
    int event_count = 0;
    int events_with_trigger = 0;
    int maxEvents = 300000000;
    int electron_count = 0;
    int positron_count = 0;
    int trigger_electron_count = 0;
    int totalFiles = 0;
    int processedFiles = 0;

    // Count total HIPO files for progress estimation
    for (const auto& dir : directories) {
        for (const auto& entry : fs::directory_iterator(dir)) {
            if (entry.is_regular_file() && entry.path().extension() == ".hipo") {
                totalFiles++;
            }
        }
    }

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
            // Check torus value for outbending
            hipo::reader reader;
            reader.open(file.c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);
            hipo::event event;
            hipo::bank RUN(factory.getSchema("RUN::config"));

            bool isOutbending = false;
            if (reader.next()) {
                reader.read(event);
                event.getStructure(RUN);
                float torus = RUN.getFloat("torus", 0);
                if (torus == 1.0) {
                    isOutbending = true;
                    processedFiles++;
                    processedFilesStream << file << endl;
                    cout << "Processing outbending file: " << file << endl;
                } else {
                    cout << "Skipping inbending file: " << file << " (torus = " << torus << ")" << endl;
                    
                    continue;
                }
            } else {
                cout << "Skipping empty file: " << file << endl;
                
                continue;
            }

            // Reset reader to start of file
            
            reader.open(file.c_str());
            reader.readDictionary(factory);

            // Loop over events
            while (reader.next() && event_count < maxEvents) {
                event_count++;
                hadTrigger = false;
                reader.read(event);
                hipo::bank PART(factory.getSchema("REC::Particle"));
                event.getStructure(PART);

                // Extract particles
                vector<ParticleData> particles(PART.getRows());
                for (int i = 0; i < PART.getRows(); ++i) {
                    particles[i] = getParticleData(i, PART);
                }

                for (const auto& p : particles) {
                    if (p.pid == 11) {
                        electron_count++;
                        if (abs(p.chi2pid) < 5 && p.vz >= -10.56 && p.vz <= 5 && abs(p.status) / 1000 == 2 && p.status < 0) {
                            trigger_electron_count++;
                            hadTrigger = true;
                        }
                    } else if (p.pid == -11) {
                        positron_count++;
                    }
                }

                if (hadTrigger) {
                    events_with_trigger++;
                    for (const auto& p : particles) {
                        if (abs(p.status) / 1000 != 2) continue;
                        if (p.pid != 211 && p.pid != 321 && p.pid != 2212) continue;
                        momentum = p.p;
                        beta = p.beta;
                        chi2pid = p.chi2pid;
                        vz=p.vz; 
                        if (p.pid == 211 ) {
                            treeEBPidPions->Fill();
                        }
                        if (p.pid == 321) {
                            treeEBPidKaons->Fill();
                        }
                        if (p.pid == 2212) {
                            treeEBPidProtons->Fill();
                        }
                    }
                }

                // Estimate remaining time every 100,000 events
                if (event_count % 100000 == 0) {
                    auto currentTime = high_resolution_clock::now();
                    auto elapsed = duration_cast<seconds>(currentTime - startTime).count();
                    double eventsPerSecond = event_count / (elapsed + 1e-6); // Avoid division by zero
                    int remainingEvents = maxEvents - event_count;
                    double remainingSeconds = remainingEvents / eventsPerSecond;
                    int hours = remainingSeconds / 3600;
                    int minutes = (remainingSeconds - hours * 3600) / 60;
                    int seconds = remainingSeconds - hours * 3600 - minutes * 60;
                    cout << "Processed " << event_count << " events (" << (event_count * 100.0 / maxEvents) << "%) "
                         << "in " << elapsed << " seconds. Estimated time remaining: "
                         << hours << "h " << minutes << "m " << seconds << "s" << endl;
                    cout << "Processed " << processedFiles << " of " << totalFiles << " files" << endl;
                }
            }
            
            if (event_count >= maxEvents) break;
        }
        if (event_count >= maxEvents) break;
    }

    // Final time estimation
    auto endTime = high_resolution_clock::now();
    auto totalElapsed = duration_cast<seconds>(endTime - startTime).count();
    cout << "Total events processed: " << event_count << endl;
    cout << "Total electron count: " << electron_count << endl;
    cout << "Total positron count: " << positron_count << endl;
    cout << "Total trigger electron count: " << trigger_electron_count << endl;
    cout << "Events with trigger electrons: " << events_with_trigger << endl;
    cout << "Total files processed: " << processedFiles << " of " << totalFiles << endl;
    cout << "Total processing time: " << totalElapsed / 3600 << "h "
         << (totalElapsed % 3600) / 60 << "m " << (totalElapsed % 60) << "s" << endl;

    // Writput
    outputFile->Write();
    outputFile->Close();
    delete outputFile;
   

    return 0;
}