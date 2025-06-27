#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <chrono>
#include <map>
#include <cmath>
#include "reader.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "CCDB/Calibration.h"
#include "CCDB/CalibrationGenerator.h"

using namespace std;
using namespace std::chrono;
namespace fs = std::filesystem;

using IndexMap = std::map<int, std::vector<int>>;

const double C = 29.9792458; // Speed of light in cm/ns
const int FTOF_DETECTOR = 12;

const float MASS_PION = 0.139570;

map<tuple<int, int, int>, float> tresCache;

struct ParticleData {
    int pid, status, charge;
    float p, beta, chi2pid, vz, vt;
    vector<tuple<int, int, int>> hits;
    int hit_sector, hit_layer, hit_component;
    ParticleData() : hit_sector(-1), hit_layer(-1), hit_component(-1) {}
};

ParticleData getParticleData(int partIdx, hipo::bank& PART, hipo::bank& SCIN, IndexMap& scinMap) {
    ParticleData pd;
    pd.pid = PART.getInt("pid", partIdx);
    float px = PART.getFloat("px", partIdx);
    float py = PART.getFloat("py", partIdx);
    float pz = PART.getFloat("pz", partIdx);
    pd.p = sqrt(px * px + py * py + pz * pz); 
    pd.beta = PART.getFloat("beta", partIdx);
    pd.chi2pid = PART.getFloat("chi2pid", partIdx);
    pd.vz = PART.getFloat("vz", partIdx);
    pd.vt = PART.getFloat("vt", partIdx);
    pd.status = PART.getShort("status", partIdx);
    pd.charge = PART.getByte("charge", partIdx); // Added charge field

    if (scinMap.find(partIdx) != scinMap.end()) {
        for (int iScin : scinMap[partIdx]) {
            int detector = SCIN.getByte("detector", iScin);
            int layer = SCIN.getByte("layer", iScin);
            if (detector == FTOF_DETECTOR && (layer == 1 || layer == 2 || layer == 3)) {
                pd.hits.push_back({detector, layer, iScin});
            } 
        }
    }

    return pd;
}

void loadCCDBParams(int runNumber) {
    ccdb::Calibration *calib = ccdb::CalibrationGenerator::CreateCalibration(
        "mysql://clas12reader@clasdb.jlab.org/clas12", runNumber, "default");
    if (!calib) {
        cerr << "Failed to connect to CCDB!" << endl;
        return;
    }

    std::vector<std::vector<double>> tresValues;
    if (calib->GetCalib(tresValues, "/calibration/ftof/tres")) {
        tresCache.clear();
        for (const auto& row : tresValues) {
            int sector = static_cast<int>(row[0]);
            int layer = static_cast<int>(row[1]);
            int component = static_cast<int>(row[2]);
            float tres = row[3];
            tresCache[make_tuple(sector, layer, component)] = tres;
        }
        cout << "Cached " << tresCache.size() << " tres values for run " << runNumber << endl;
    } else {
        cerr << "Failed to load /calibration/ftof/tres for run " << runNumber << endl;
    }
    delete calib;
}

IndexMap loadMapByIndex(hipo::bank& fromBank, const char* idxVarName) {
    IndexMap map;
    for (int i = 0; i < fromBank.getRows(); ++i) {
        int pindex = fromBank.getInt(idxVarName, i);
        if (pindex >= 0) map[pindex].push_back(i);
    }
    return map;
}

float computeDeltaT(ParticleData& p, hipo::bank& SCIN) {
    float delta_t = 99999.0f;
    float mass = MASS_PION;
    float beta_theory = p.p / sqrt(p.p * p.p + mass * mass);
    if (beta_theory <= 0) return delta_t;

    vector<int> layer_priority = {2, 1, 3};
    for (int layer : layer_priority) {
        for (const auto& hit : p.hits) {
            int hit_detector, hit_layer, scin_idx;
            tie(hit_detector, hit_layer, scin_idx) = hit;
            if (hit_detector == FTOF_DETECTOR && hit_layer == layer) {
                float hit_time = SCIN.getFloat("time", scin_idx);
                float hit_path = SCIN.getFloat("path", scin_idx);
                if (hit_time == 0.0f || hit_path == 0.0f || isinf(hit_time)) continue;
                p.hit_sector = SCIN.getByte("sector", scin_idx);
                p.hit_layer = hit_layer;
                p.hit_component = SCIN.getInt("component", scin_idx);
                float vertex_time = hit_time - hit_path / (C * beta_theory);
                delta_t = vertex_time - p.vt;
                return delta_t;
            }
        }
    }
    return delta_t;
}

float computeChi2PidPion(const ParticleData& p, float delta_t, float& sigma_out) {
    if (delta_t >= 99999.0f || p.hit_sector < 1 || p.hit_layer < 1 || p.hit_component < 1) {
        sigma_out = -1.0f;
        return 99999.0f;
    }
    if (!tresCache.count(make_tuple(p.hit_sector, p.hit_layer, p.hit_component))) {
        sigma_out = -1.0f;
        return 99999.0f;
    }
    float sigma = tresCache[make_tuple(p.hit_sector, p.hit_layer, p.hit_component)];
    if (sigma <= 0) {
        sigma_out = -1.0f;
        return 99999.0f;
    }
    float q = delta_t / sigma;
    sigma_out = sigma;
    return q;
}

int main() {
    auto startTime = high_resolution_clock::now();
    
    ofstream processedFilesStream("processed_hipo_files.txt");
    if (!processedFilesStream.is_open()) {
        cerr << "Error: Could not open processed_hipo_files.txt" << endl;
        return 1;
    }

    ofstream csvFile("detector_info.csv");
    if (!csvFile.is_open()) {
        cerr << "Error: Could not open detector_info.csv" << endl;
        return 1;
    }
    csvFile << "event_number,pid,charge,orig_chi2pid,recomputed_chi2pid,sector,layer,component,sigma,dt\n"; // Added charge to CSV

    TFile* outputFile = new TFile("pkptreeCxC_9_test_with_pos_charge_pion_assumed.root", "RECREATE");
    TTree* treeEBPidPions = new TTree("EB_pid_pions", "Pions with pid == 211");
    TTree* treeEBPidKaons = new TTree("EB_pid_kaons", "Kaons with pid == 321");
    TTree* treeEBPidProtons = new TTree("EB_pid_protons", "Protons with pid == 2212");
    TTree* treeEBAllPionAssumed = new TTree("EB_all_pion_assumed", "All particles assumed as pions (PID 211, 321, 2212)");
    TTree* treeEBPosPionAssumed = new TTree("EB_pos_pion_assumed", "Positive charge particles assumed as pions");

    float orig_chi2pid, recomputed_chi2pid, momentum, beta, vz, dt;
    // Branches for separate trees
    treeEBPidPions->Branch("orig_chi2pid", &orig_chi2pid, "orig_chi2pid/F");
    treeEBPidPions->Branch("recomputed_chi2pid", &recomputed_chi2pid, "recomputed_chi2pid/F");
    treeEBPidPions->Branch("p", &momentum, "p/F");
    treeEBPidPions->Branch("beta", &beta, "beta/F");
    treeEBPidPions->Branch("vz", &vz, "vz/F");
    treeEBPidPions->Branch("dt", &dt, "dt/F");

    treeEBPidKaons->Branch("orig_chi2pid", &orig_chi2pid, "orig_chi2pid/F");
    treeEBPidKaons->Branch("recomputed_chi2pid", &recomputed_chi2pid, "recomputed_chi2pid/F");
    treeEBPidKaons->Branch("p", &momentum, "p/F");
    treeEBPidKaons->Branch("beta", &beta, "beta/F");
    treeEBPidKaons->Branch("vz", &vz, "vz/F");
    treeEBPidKaons->Branch("dt", &dt, "dt/F");
    
    treeEBPidProtons->Branch("orig_chi2pid", &orig_chi2pid, "orig_chi2pid/F");
    treeEBPidProtons->Branch("recomputed_chi2pid", &recomputed_chi2pid, "recomputed_chi2pid/F");
    treeEBPidProtons->Branch("p", &momentum, "p/F");
    treeEBPidProtons->Branch("beta", &beta, "beta/F");
    treeEBPidProtons->Branch("vz", &vz, "vz/F");
    treeEBPidProtons->Branch("dt", &dt, "dt/F");

    // Branches for the existing combined tree (PID 211, 321, 2212)
    treeEBAllPionAssumed->Branch("recomputed_chi2pid", &recomputed_chi2pid, "recomputed_chi2pid/F");
    treeEBAllPionAssumed->Branch("p", &momentum, "p/F");
    treeEBAllPionAssumed->Branch("beta", &beta, "beta/F");
    treeEBAllPionAssumed->Branch("vz", &vz, "vz/F");
    treeEBAllPionAssumed->Branch("dt", &dt, "dt/F");

    // Branches for the new positive charge particles tree
    treeEBPosPionAssumed->Branch("recomputed_chi2pid", &recomputed_chi2pid, "recomputed_chi2pid/F");
    treeEBPosPionAssumed->Branch("p", &momentum, "p/F");
    treeEBPosPionAssumed->Branch("beta", &beta, "beta/F");
    treeEBPosPionAssumed->Branch("vz", &vz, "vz/F");
    treeEBPosPionAssumed->Branch("dt", &dt, "dt/F");

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

    int event_count = 0;
    int events_with_trigger = 0;
    int maxEvents = 1000000;
    int electron_count = 0;
    int positron_count = 0;
    int trigger_electron_count = 0;
    int totalFiles = 0;
    int processedFiles = 0;
    int particle_count_pid = 0; // Counter for PID-based tree
    int particle_count_charge = 0; // Counter for charge-based tree
    const int maxParticles = 100;

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
            hipo::reader reader;
            reader.open(file.c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);
            hipo::event event;
            hipo::bank RUN(factory.getSchema("RUN::config"));
            hipo::bank PART(factory.getSchema("REC::Particle"));
            hipo::bank SCIN(factory.getSchema("REC::Scintillator"));

            bool isOutbending = false;
            int runNumber = -1;
            if (reader.gotoEvent(0)) {
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
                runNumber = RUN.getInt("run", 0);
                loadCCDBParams(runNumber);
            } else {
                cout << "Skipping empty file: " << file << endl;
                continue;
            }

            reader.rewind();

            while (reader.next() && event_count < maxEvents) {
                event_count++;
                reader.read(event);
                event.getStructure(PART);
                event.getStructure(SCIN);
                IndexMap scinMap = loadMapByIndex(SCIN, "pindex");

                vector<ParticleData> particles(PART.getRows());
                for (int i = 0; i < PART.getRows(); ++i) {
                    particles[i] = getParticleData(i, PART, SCIN, scinMap);
                }

                bool hasTriggerElectron = false;
                for (const auto& pd : particles) {
                    if (pd.pid == 11) {
                        electron_count++;
                        if (abs(pd.chi2pid) < 5 && pd.vz >= -10.56 && pd.vz <= 5 && abs(pd.status) / 1000 == 2 && pd.status < 0) {
                            trigger_electron_count++;
                            hasTriggerElectron = true;
                        }
                    } else if (pd.pid == -11) {
                        positron_count++;
                    }
                }

                if (hasTriggerElectron) {
                    events_with_trigger++;
                    for (auto& pd : particles) {
                        if (abs(pd.status) / 1000 != 2) continue;
                       
                        orig_chi2pid = pd.chi2pid;
                        momentum = pd.p;
                        beta = pd.beta;
                        vz = pd.vz;
                        dt = computeDeltaT(pd, SCIN);
                        float sigma;
                        recomputed_chi2pid = computeChi2PidPion(pd, dt, sigma);

                        // For PID-based tree (211, 321, 2212)
                        if (pd.pid == 211 || pd.pid == 321 || pd.pid == 2212) {

                            csvFile << event_count << "," << pd.pid << "," << pd.charge << "," << orig_chi2pid << "," << recomputed_chi2pid << ","
                                    << pd.hit_sector << "," << pd.hit_layer << "," << pd.hit_component << ","
                                    << sigma << "," << dt << "\n";

                            treeEBAllPionAssumed->Fill();
                            particle_count_pid++;

                            // Fill the separate trees as before
                            if (pd.pid == 211) {
                                treeEBPidPions->Fill();
                            } else if (pd.pid == 321) {
                                treeEBPidKaons->Fill();
                            } else if (pd.pid == 2212) {
                                treeEBPidProtons->Fill();
                            }
                        }

                        // For charge-based tree (charge > 0)
                        if (pd.charge > 0) {
                            treeEBPosPionAssumed->Fill();
                            particle_count_charge++;
                        }
                    }
                }

                if (event_count % 100000 == 0) {
                    auto currentTime = high_resolution_clock::now();
                    auto elapsed = duration_cast<seconds>(currentTime - startTime).count();
                    double eventsPerSecond = event_count / (elapsed + 1e-6);
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

            if (event_count >= maxEvents ) break;
        }   
        if (event_count >= maxEvents ) break;
    }

    auto endTime = high_resolution_clock::now();
    auto totalElapsed = duration_cast<seconds>(endTime - startTime).count();
    cout << "Total events processed: " << event_count << endl;
    cout << "Total particles saved in PID-based tree: " << particle_count_pid << endl;
    cout << "Total particles saved in charge-based tree: " << particle_count_charge << endl;
    cout << "Total electron count: " << electron_count << endl;
    cout << "Total positron count: " << positron_count << endl;
    cout << "Total trigger electron count: " << trigger_electron_count << endl;
    cout << "Events with trigger electrons: " << events_with_trigger << endl;
    cout << "Total files processed: " << processedFiles << " of " << totalFiles << endl;
    cout << "Total processing time: " << totalElapsed / 3600 << "h "
         << (totalElapsed % 3600) / 60 << "m " << (totalElapsed % 60) << "s" << endl;

    outputFile->Write();
    outputFile->Close();
    delete outputFile;
    processedFilesStream.close();
    csvFile.close();

    return 0;
}