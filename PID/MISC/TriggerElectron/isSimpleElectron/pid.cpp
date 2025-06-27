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
#include "TH1F.h"
#include <TStopwatch.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <limits>
#include "TStyle.h"
#include <set>
#include "CCDB/Calibration.h"
#include "CCDB/CalibrationGenerator.h"

using namespace std;
namespace fs = std::filesystem;

// Constants
const float MIN_PCAL_ENERGY = 0.06f;
const float NSIGMA_CUT = 5.0f;
const int HTCC_DETECTOR = 15;
const int LTCC_DETECTOR = 16;
const int ECAL_DETECTOR = 7;
const int FTOF_DETECTOR = 12;
const int CTOF_DETECTOR = 4;
const int CND_DETECTOR = 3;
const int BAND_DETECTOR = 21;
const set<int> PID_POSITIVE = {-11, 211, 321, 2212, 45};
const set<int> PID_NEGATIVE = {11, -211, -321, -2212};
const set<int> PID_NEUTRAL = {22, 2112};
const float SPEED_OF_LIGHT = 29.9792458; // cm/ns
const float HTCC_PION_THRESHOLD = 4.9f;
const float LTCC_PION_THRESHOLD = 3.0f;

// Particle masses (GeV/c^2)
map<int, float> PDG_MASS = {
    {11, 0.0005}, {-11, 0.0005}, {211, 0.13957018}, {-211, 0.13957018},
    {321, 0.49367716}, {-321, 0.49367716}, {2212, 0.938272046}, {-2212, 0.938272046},
    {45, 1.87705}, {22, 0.0f}, {2112, 0.939565379}
};

// Detector priorities
const vector<pair<int, vector<int>>> chargedBetaDetectors = {
    {FTOF_DETECTOR, {2, 1, 3}}, {CTOF_DETECTOR, {1}}, {ECAL_DETECTOR, {1, 4, 7}}
};
const vector<pair<int, vector<int>>> neutralBetaDetectors = {
    {ECAL_DETECTOR, {1, 4, 7}}, {CND_DETECTOR, {1, 2, 3}}, {CTOF_DETECTOR, {1}}, {BAND_DETECTOR, {1}}
};

// Structures
typedef std::map<int, std::vector<int>> IndexMap;

struct DetectorHit {
    int detector;
    int layer;
    float time;
    float path;
};

struct ParticleData {
    int charge;
    float p;
    int sector;
    float nphe_htcc;
    float nphe_ltcc;
    float energy_pcal;
    float energy_total;
    int pid;
    float vt;
    float beta; // REC::Particle beta (not used for neutrals)
    vector<DetectorHit> hits;
    bool is_trigger;
    int assigned_pid;
    float start_time;
    int status;
};

struct SamplingFractionParams {
    float sf[4];
    float sfs[4];
};
std::vector<SamplingFractionParams> sfParams(7);
float htccNpheCut = 2.0f;
float ltccNpheCut = 2.0f;
float neutronMaxBeta = 0.9f;
float cndNeutronMaxBeta = 0.8f;

// CCDB Loading (unchanged)
void loadCCDBParams() {
    ccdb::Calibration *calib = ccdb::CalibrationGenerator::CreateCalibration("mysql://clas12reader@clasdb.jlab.org/clas12", 18536, "default");
    if (!calib) { cerr << "Failed to create CCDB Calibration object!" << endl; exit(1); }

    std::vector<std::vector<double>> sfValues;
    if (!calib->GetCalib(sfValues, "/calibration/eb/electron_sf")) {
        cerr << "Failed to get electron_sf! Using defaults." << endl;
        for (int sector = 1; sector <= 6; sector++) {
            sfParams[sector].sf[0] = 0.23372 + (sector - 1) * 0.01;
            sfParams[sector].sf[1] = 1.0;
            sfParams[sector].sf[2] = -0.01726 - (sector - 1) * 0.005;
            sfParams[sector].sf[3] = 0.00050 + (sector - 1) * 0.0001;
            sfParams[sector].sfs[0] = 0.01763 - (sector - 1) * 0.001;
            sfParams[sector].sfs[1] = 1.0;
            sfParams[sector].sfs[2] = 0.0;
            sfParams[sector].sfs[3] = 0.0;
        }
    } else {
        for (const auto& row : sfValues) {
            int sector = static_cast<int>(row[0]);
            if (sector < 1 || sector > 6 || row.size() < 11) continue;
            sfParams[sector].sf[0] = row[3]; sfParams[sector].sf[1] = row[4];
            sfParams[sector].sf[2] = row[5]; sfParams[sector].sf[3] = row[6];
            sfParams[sector].sfs[0] = row[7]; sfParams[sector].sfs[1] = row[8];
            sfParams[sector].sfs[2] = row[9]; sfParams[sector].sfs[3] = row[10];
        }
    }

    std::vector<std::vector<double>> htccValues;
    if (!calib->GetCalib(htccValues, "/calibration/eb/htcc_matching") || htccValues.empty() || htccValues[0].size() < 7) {
        cerr << "Failed to get HTCC nphe cut! Using 2.0" << endl;
        htccNpheCut = 2.0f;
    } else {
        htccNpheCut = htccValues[0][6];
    }
    ltccNpheCut = 2.0f;

    std::vector<std::vector<double>> neutronValues;
    if (!calib->GetCalib(neutronValues, "/calibration/eb/neutron_maxbeta") || neutronValues.empty()) {
        cerr << "Failed to get neutron_maxbeta! Using 0.9" << endl;
        neutronMaxBeta = 0.9f;
    } else {
        neutronMaxBeta = neutronValues[0][0];
    }
    if (!calib->GetCalib(neutronValues, "/calibration/eb/cnd_neutron_maxbeta") || neutronValues.empty()) {
        cerr << "Failed to get cnd_neutron_maxbeta! Using 0.8" << endl;
        cndNeutronMaxBeta = 0.8f;
    } else {
        cndNeutronMaxBeta = neutronValues[0][0];
    }
    delete calib;
}

// Sampling Fraction Functions (unchanged)
float getSamplingFractionMean(int sector, float measuredEnergy) {
    if (sector < 1 || sector > 6) return 0.0f;
    const auto& p = sfParams[sector].sf;
    return p[0] * (p[1] + p[2] / measuredEnergy + p[3] * std::pow(measuredEnergy, -2));
}

float getSamplingFractionSigma(int sector, float measuredEnergy) {
    if (sector < 1 || sector > 6) return 0.0f;
    const auto& p = sfParams[sector].sfs;
    return p[0] * (p[1] + p[2] / measuredEnergy + p[3] * std::pow(measuredEnergy, -2));
}

float getSamplingFractionNSigma(float samplingFraction, float mean, float sigma) {
    if (sigma == 0) return 0.0f;
    return (samplingFraction - mean) / sigma;
}

// Utility Functions
IndexMap loadMapByIndex(hipo::bank& fromBank, const char* idxVarName) {
    IndexMap map;
    for (int i = 0; i < fromBank.getRows(); ++i) {
        int iTo = fromBank.getInt(idxVarName, i);
        map[iTo].push_back(i);
    }
    return map;
}

ParticleData getParticleData(int partIdx, hipo::bank& PART, hipo::bank& CHER, hipo::bank& CALO, hipo::bank& SCIN,
                             IndexMap& cherMap, IndexMap& caloMap, IndexMap& scinMap) {
    ParticleData pd = {0, 0.0f, -1, 0.0f, 0.0f, 0.0f, 0.0f, 0, 0.0f, -99.0f, {}, false, 0, 0.0f, 0};
    pd.charge = PART.getByte("charge", partIdx);
    pd.pid = PART.getInt("pid", partIdx);
    float px = PART.getFloat("px", partIdx);
    float py = PART.getFloat("py", partIdx);
    float pz = PART.getFloat("pz", partIdx);
    pd.p = std::sqrt(px * px + py * py + pz * pz);
    pd.vt = PART.getFloat("vt", partIdx);
    pd.beta = PART.getFloat("beta", partIdx);
    pd.status = PART.getShort("status", partIdx);
    pd.is_trigger = (pd.status < 0);
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
                pd.hits.push_back({ECAL_DETECTOR, layer, time, path});
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
            float time = SCIN.getFloat("time", iScinRow);
            float path = SCIN.getFloat("path", iScinRow);
            pd.hits.push_back({detector, layer, time, path});
        }
    }
    return pd;
}

bool isSimpleElectron(const ParticleData& pd) {
    if (pd.charge != -1 && pd.charge != 1) return false;
    if (pd.sector < 1 || pd.sector > 6) return false;
    if (pd.nphe_htcc < htccNpheCut) return false;
    if (pd.energy_pcal < MIN_PCAL_ENERGY) return false;
    float sf = pd.energy_total / pd.p;
    float mean = getSamplingFractionMean(pd.sector, pd.energy_total);
    float sigma = getSamplingFractionSigma(pd.sector, pd.energy_total);
    float sfNSigma = getSamplingFractionNSigma(sf, mean, sigma);
    return std::abs(sfNSigma) <= NSIGMA_CUT;
}

float getTheoryBeta(float p, float mass) {
    return p / sqrt(p * p + mass * mass);
}

bool hasHit(const ParticleData& p, int detector, int layer = -1) {
    for (const auto& hit : p.hits) {
        if (hit.detector == detector && (layer == -1 || hit.layer == layer)) {
            return true;
        }
    }
    return false;
}

float getNeutralBeta(const ParticleData& p, int detector, const std::vector<int>& layers, float startTime, float vertexTime) {
    float beta = -1.0f;
    float correctedStartTime = (vertexTime > 0) ? vertexTime : startTime; // Use vertex if valid
    for (int layer : layers) {
        for (const auto& hit : p.hits) {
            if (hit.detector == detector && hit.layer == layer) {
                float deltaT = hit.time - correctedStartTime;
                beta = (deltaT > 0) ? hit.path / (deltaT * SPEED_OF_LIGHT) : -1.0f;
                cout << "Beta Calc: pindex=" << p.pid << ", detector=" << detector 
                     << ", layer=" << layer << ", time=" << hit.time << ", path=" << hit.path 
                     << ", deltaT=" << deltaT << ", beta=" << beta << endl;
                if (beta >= -0.1 && beta <= 1.1) return beta; // Allow β < 0 (neutrons)
                cout << "Unphysical beta=" << beta << ", using REC::Particle.beta=" << p.beta << endl;
                return p.beta; // Fallback
            }
        }
    }
    cout << "No valid hit for detector=" << detector << ", using REC::Particle.beta=" << p.beta << endl;
    return p.beta; // Fallback
}

int bestPidFromTiming(const ParticleData& p, float eventStartTime) {
    float startTime = eventStartTime;
    if (p.charge != 0 && p.vt != -9999.0f) startTime = p.vt;

    if (p.charge == 0) {
        float beta = -1.0f;
        int assigned_pid = 0;
        string detector_used = "None";
        bool is_forward = (abs(p.status) / 2000 == 1);
        bool is_central = (abs(p.status) / 4000 == 1);
        bool is_band = (abs(p.status) / 8000 == 1);

        for (const auto& det : neutralBetaDetectors) {
            if (hasHit(p, det.first)) {
                beta = getNeutralBeta(p, det.first, det.second, startTime, p.vt);
                detector_used = (det.first == ECAL_DETECTOR) ? "ECAL" :
                               (det.first == CND_DETECTOR) ? "CND" :
                               (det.first == CTOF_DETECTOR) ? "CTOF" : "BAND";
                break;
            }
        }

        if (beta >= -0.1) { // Allow negative β for neutrons
            if (detector_used == "BAND" && beta < 0) assigned_pid = 0; // BAND unassignment
            else if (hasHit(p, ECAL_DETECTOR)) {
                assigned_pid = (beta < neutronMaxBeta) ? 2112 : 22;
            } else if (hasHit(p, CND_DETECTOR) || hasHit(p, CTOF_DETECTOR)) {
                assigned_pid = (beta < cndNeutronMaxBeta) ? 2112 : 22;
            } else if (hasHit(p, BAND_DETECTOR)) {
                assigned_pid = (beta < 0.9f) ? 2112 : 0;
            } else {
                assigned_pid = (beta < neutronMaxBeta) ? 2112 : 22; // Default
            }
        }

        cout << "Neutral: pindex=" << p.pid << ", beta=" << beta 
             << ", detector=" << detector_used << ", PID=" << assigned_pid 
             << ", status=" << p.status << ", startTime=" << startTime 
             << ", REC::Particle.beta=" << p.beta << endl;
        return assigned_pid;
    } else {
        // Charged logic unchanged
        set<int> hypotheses = (p.charge > 0) ? PID_POSITIVE : PID_NEGATIVE;
        float min_dt = numeric_limits<float>::max();
        int best_pid = 0;
        for (const auto& det : chargedBetaDetectors) {
            for (int layer : det.second) {
                if (hasHit(p, det.first, layer)) {
                    for (int pid : hypotheses) {
                        if (abs(pid) == 11) continue;
                        float mass = PDG_MASS[pid];
                        float beta_theory = p.p / sqrt(p.p * p.p + mass * mass);
                        for (const auto& hit : p.hits) {
                            if (hit.detector == det.first && hit.layer == layer) {
                                float vt = hit.time - hit.path / (SPEED_OF_LIGHT * beta_theory);
                                float dt = fabs(vt - startTime);
                                if (dt < min_dt) {
                                    min_dt = dt;
                                    best_pid = pid;
                                }
                            }
                        }
                    }
                    if (best_pid != 0) {
                        cout << "Charged: PID=" << best_pid << ", beta_theory=" 
                             << p.p / sqrt(p.p * p.p + PDG_MASS[best_pid] * PDG_MASS[best_pid]) 
                             << ", startTime=" << startTime << endl;
                        return best_pid;
                    }
                }
            }
        }
        cout << "Charged: PID=" << best_pid << ", startTime=" << startTime << endl;
        return best_pid;
    }
}

void assignPids(vector<ParticleData>& particles, float eventStartTime) {
    for (auto& p : particles) {
        if (p.is_trigger && p.pid == 11) {
            p.assigned_pid = 11;
            continue;
        }
        p.assigned_pid = 0;
        bool is_elec = isSimpleElectron(p);
        int pid_from_timing = bestPidFromTiming(p, eventStartTime);
        bool htcc_signal = p.nphe_htcc > htccNpheCut;
        bool ltcc_signal = p.nphe_ltcc > ltccNpheCut;
        bool htcc_pion = p.p > HTCC_PION_THRESHOLD;
        bool ltcc_pion = p.p > LTCC_PION_THRESHOLD;

        set<int> hypotheses = (p.charge > 0) ? PID_POSITIVE : (p.charge < 0) ? PID_NEGATIVE : PID_NEUTRAL;
        for (int pid : hypotheses) {
            float theory_beta = getTheoryBeta(p.p, PDG_MASS[pid]);
            bool timing_check = (pid == pid_from_timing && theory_beta > 0);
            switch (abs(pid)) {
                case 11:
                    if (is_elec) p.assigned_pid = (p.charge > 0) ? -11 : 11;
                    break;
                case 211:
                    if (timing_check && !is_elec) p.assigned_pid = pid;
                    else if (!is_elec && htcc_signal && htcc_pion) p.assigned_pid = pid;
                    break;
                case 321:
                    if (timing_check && !is_elec) {
                        if (ltcc_signal && ltcc_pion) p.assigned_pid = (pid > 0) ? 211 : -211;
                        else p.assigned_pid = pid;
                    }
                    break;
                case 2212:
                    if (timing_check && !is_elec) {
                        if (ltcc_signal && ltcc_pion) p.assigned_pid = (pid > 0) ? 211 : -211;
                        else p.assigned_pid = pid;
                    }
                    break;
                case 45:
                    if (timing_check && !is_elec) p.assigned_pid = pid;
                    break;
                case 2112:
                    if (timing_check) p.assigned_pid = pid;
                    break;
                case 22:
                    if (timing_check) p.assigned_pid = pid; // Relaxed PCAL
                    break;
            }
            if (p.assigned_pid != 0) break;
        }
        cout << "Assigned PID=" << p.assigned_pid << " for status=" << p.status 
             << ", REC::Particle.pid=" << p.pid << endl;
    }
}

int main() {
    TStopwatch timer;
    timer.Start();

    loadCCDBParams();

    ifstream inputFile("directories.txt");
    if (!inputFile.is_open()) { cerr << "Error: Could not open directories.txt" << endl; return 1; }
    vector<string> directories;
    string dir;
    while (getline(inputFile, dir)) {
        if (!dir.empty()) directories.push_back(dir);
    }
    inputFile.close();
    if (directories.empty()) { cerr << "Error: No directories found" << endl; return 1; }

    map<int, int> pid_counts;
    map<int, string> pid_names = {
        {11, "Electron"}, {-11, "Positron"}, {211, "Pi+"}, {-211, "Pi-"},
        {321, "K+"}, {-321, "K-"}, {2212, "Proton"}, {-2212, "Anti-Proton"},
        {45, "Deuteron"}, {22, "Photon"}, {2112, "Neutron"}, {0, "Unidentified"}
    };

    int total_events = 0;
    const int maxEvents = 100000;

    for (const auto& dir : directories) {
        vector<string> hipoFiles;
        for (const auto& entry : fs::directory_iterator(dir)) {
            if (entry.path().extension() == ".hipo") hipoFiles.push_back(entry.path().string());
        }
        if (hipoFiles.empty()) { cout << "  No .hipo files in " << dir << endl; continue; }

        for (const auto& file : hipoFiles) {
            cout << "  Opening " << file << endl;
            hipo::reader reader;
            reader.open(file.c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);

            while (reader.next() && total_events < maxEvents) {
                total_events++;

                hipo::event event;
                reader.read(event);

                hipo::bank EVENT(factory.getSchema("REC::Event"));
                hipo::bank PART(factory.getSchema("REC::Particle"));
                hipo::bank CHER(factory.getSchema("REC::Cherenkov"));
                hipo::bank CALO(factory.getSchema("REC::Calorimeter"));
                hipo::bank SCIN(factory.getSchema("REC::Scintillator"));
                event.getStructure(EVENT);
                event.getStructure(PART);
                event.getStructure(CHER);
                event.getStructure(CALO);
                event.getStructure(SCIN);

                IndexMap cherMap = loadMapByIndex(CHER, "pindex");
                IndexMap caloMap = loadMapByIndex(CALO, "pindex");
                IndexMap scinMap = loadMapByIndex(SCIN, "pindex");

                float event_start_time = EVENT.getFloat("startTime", 0);
                vector<ParticleData> particles(PART.getRows());
                int valid_particles = 0;
                for (int i = 0; i < PART.getRows(); i++) {
                    ParticleData pd = getParticleData(i, PART, CHER, CALO, SCIN, cherMap, caloMap, scinMap);
                    if (pd.status == 0) continue;
                    particles[valid_particles++] = pd;
                }
                particles.resize(valid_particles);

                assignPids(particles, event_start_time);

                for (const auto& p : particles) {
                    pid_counts[p.assigned_pid]++;
                }
            }
            if (total_events >= maxEvents) break;
        }
        if (total_events >= maxEvents) break;
    }

    cout << "Total events: " << total_events << endl;
    for (const auto& [pid, name] : pid_names) {
        cout << name << ": " << pid_counts[pid] << endl;
    }
    cout << "Time elapsed: " << timer.RealTime() << " seconds" << endl;

    return 0;
}