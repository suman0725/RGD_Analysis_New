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
#include "CCDB/Calibration.h"
#include "CCDB/CalibrationGenerator.h"

using namespace std;
using namespace std::chrono;
namespace fs = std::filesystem;

// Define IndexMap as std::map<int, std::vector<int>>
using IndexMap = std::map<int, std::vector<int>>;

const double C = 29.9792458; // Speed of light in cm/ns
const int FTOF_DETECTOR = 12; // FTOF detector ID

// PDG masses in GeV/c^2
const float MASS_PION = 0.139570; // Pion mass (pid=211)
const float MASS_KAON = 0.493677; // Kaon mass (pid=321)

// Cache for FTOF time resolution (tres)
map<tuple<int, int, int>, float> tresCache; // (sector, layer, component) -> tres

// Load CCDB parameters
void loadCCDBParams() {
    ccdb::Calibration *calib = ccdb::CalibrationGenerator::CreateCalibration(
        "mysql://clas12reader@clasdb.jlab.org/clas12", 18536, "default");
    if (!calib) {
        cerr << "Failed to connect to CCDB! Using default sigma=0.1 ns" << endl;
        return;
    }

    std::vector<std::vector<double>> tresValues;
    if (calib->GetCalib(tresValues, "/calibration/ftof/tres")) {
        for (const auto& row : tresValues) {
            int sector = static_cast<int>(row[0]);
            int layer = static_cast<int>(row[1]);
            int component = static_cast<int>(row[2]);
            float tres = row[3];
            tresCache[make_tuple(sector, layer, component)] = tres;
        }
        cout << "Cached " << tresCache.size() << " tres values from CCDB" << endl;
    } else {
        cerr << "Failed to load /calibration/ftof/tres from CCDB! Using default sigma=0.1 ns" << endl;
    }

    delete calib;
}

// Define ParticleData struct
struct ParticleData {
    int pid, status;
    float p, beta, chi2pid, vz, vt;
    vector<tuple<int, int, int>> hits; // {detector, layer, index in scin bank}
    int hit_sector, hit_layer, hit_component; // FTOF hit used for dt
    ParticleData() : hit_sector(-1), hit_layer(-1), hit_component(-1) {} // Initialize to invalid
};

// Improved IndexMap loading
IndexMap loadMapByIndex(hipo::bank& fromBank, const char* idxVarName) {
    IndexMap map;
    for (int i = 0; i < fromBank.getRows(); ++i) {
        int pindex = fromBank.getInt(idxVarName, i);
        if (pindex >= 0) map[pindex].push_back(i);
    }
    return map;
}

// Get particle data with FTOF hit info
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

// Compute deltaT for a kaon under the pion hypothesis
float computeDeltaT(ParticleData& p, hipo::bank& SCIN) {
    float delta_t = 99999.0f;
    float mass = MASS_PION;
    float beta_theory = p.p / sqrt(p.p * p.p + mass * mass);
    if (beta_theory <= 0 || beta_theory >= 1.0f) return delta_t;

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

// Compute chi2pid under pion hypothesis
float computeChi2PidPion(const ParticleData& p, float delta_t) {
    if (delta_t >= 99999.0f || p.hit_sector < 1 || p.hit_layer < 1 || p.hit_component < 1) {
        return 99999.0f;
    }
    float sigma = 0.1f; // Default sigma
    if (tresCache.count(make_tuple(p.hit_sector, p.hit_layer, p.hit_component))) {
        sigma = tresCache[make_tuple(p.hit_sector, p.hit_layer, p.hit_component)];
    }
    if (sigma <= 0) return 99999.0f;
    float q = delta_t / sigma;
    return q ;
}

int main() {
    auto startTime = high_resolution_clock::now();

    // Load CCDB parameters
    loadCCDBParams();

    // Open output file
    TFile* outputFile = new TFile("pkptreeCxC_7.root", "RECREATE");
    TTree* treeEBPidPions = new TTree("EB_pid_pions", "Pions with pid == 211");
    TTree* treeEBPidKaons = new TTree("EB_pid_kaons", "Kaons with pid == 321");

    float chi2pid, momentum, beta, vz, dt, chi2pid_pion;
    treeEBPidPions->Branch("chi2pid", &chi2pid, "chi2pid/F");
    treeEBPidPions->Branch("p", &momentum, "p/F");
    treeEBPidPions->Branch("beta", &beta, "beta/F");
    treeEBPidPions->Branch("vz", &vz, "vz/F");
    treeEBPidKaons->Branch("chi2pid", &chi2pid, "chi2pid/F");
    treeEBPidKaons->Branch("p", &momentum, "p/F");
    treeEBPidKaons->Branch("beta", &beta, "beta/F");
    treeEBPidKaons->Branch("vz", &vz, "vz/F");
    treeEBPidKaons->Branch("dt", &dt, "dt/F");
    treeEBPidKaons->Branch("chi2pid_pion", &chi2pid_pion, "chi2pid_pion/F");

    // Read directories
    ifstream inputFile("directories.txt");
    if (!inputFile.is_open()) {
        cerr << "Error: Could not open directories.txt" << endl;
        return 1;
    }
    vector<string> directories;
    string dir;
    while (getline(inputFile, dir)) if (!dir.empty()) directories.push_back(dir);

    int processedFiles = 0;

    for (const auto& dir : directories) {
        for (const auto& entry : fs::directory_iterator(dir)) {
            if (entry.path().extension() != ".hipo") continue;

            string file = entry.path().string();
            hipo::reader reader;
            reader.open(file.c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);
            hipo::event event;
            hipo::bank RUN(factory.getSchema("RUN::config"));
            hipo::bank PART(factory.getSchema("REC::Particle"));
            hipo::bank SCIN(factory.getSchema("REC::Scintillator"));

            if (reader.gotoEvent(0)) {
                reader.read(event);
                event.getStructure(RUN);
                float torus = RUN.getFloat("torus", 0);
                if (torus != 1) {
                    cout << "Skipping inbending file: " << file << " (torus = " << torus << ")" << endl;
                    continue;
                }
                cout << "Processing outbending file: " << file << " (torus = " << torus << ")" << endl;
                processedFiles++;

                reader.rewind();
                int event_count = 0;
                while (reader.next() && event_count < 1000) {
                    reader.read(event);
                    event.getStructure(PART);
                    event.getStructure(SCIN);

                    IndexMap scinMap = loadMapByIndex(SCIN, "pindex");

                    bool hasTriggerElectron = false;
                    for (int i = 0; i < PART.getRows(); ++i) {
                        ParticleData pd = getParticleData(i, PART, SCIN, scinMap);
                        if (pd.pid == 11 && abs(pd.chi2pid) < 5 && pd.vz >= -10.56 && pd.vz <= 5 &&
                            abs(pd.status) / 1000 == 2 && pd.status < 0) {
                            hasTriggerElectron = true;
                            break;
                        }
                    }

                    if (hasTriggerElectron) {
                        for (int i = 0; i < PART.getRows(); ++i) {
                            ParticleData pd = getParticleData(i, PART, SCIN, scinMap);
                            if (pd.hits.empty()) continue;

                            chi2pid = pd.chi2pid;
                            momentum = pd.p;
                            beta = pd.beta;
                            vz = pd.vz;

                            if (pd.pid == 211) {
                                treeEBPidPions->Fill();
                            } else if (pd.pid == 321) {
                                dt = computeDeltaT(pd, SCIN);
                                chi2pid_pion = computeChi2PidPion(pd, dt);
                                treeEBPidKaons->Fill();
                            }
                        }
                    }
                    event_count++;
                }
            } else {
                cout << "Skipping empty file: " << file << endl;
            }
        }
    }

    if (processedFiles == 0) {
        cerr << "Error: No outbending files (torus = 1) found in directories.txt" << endl;
    }

    outputFile->Write();
    outputFile->Close();
    delete outputFile;

    auto endTime = high_resolution_clock::now();
    cout << "Processed " << processedFiles << " files in " << duration_cast<seconds>(endTime - startTime).count() << "s" << endl;
    return 0;
}