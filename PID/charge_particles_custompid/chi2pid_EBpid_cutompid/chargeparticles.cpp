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
#include "CCDB/Calibration.h"
#include "CCDB/CalibrationGenerator.h"

using namespace std;
namespace fs = std::filesystem;

// Define IndexMap as a map of int to vector<int>
using IndexMap = std::map<int, std::vector<int>>;

// Constants
const float MIN_PCAL_ENERGY = 0.06f;
const float NSIGMA_CUT = 5.0f;
const int HTCC_DETECTOR = 15;
const int LTCC_DETECTOR = 16;
const int ECAL_DETECTOR = 7;
const int FTOF_DETECTOR = 12;
const int CTOF_DETECTOR = 4;
const set<int> PID_POSITIVE = {-11, 211, 321, 2212, 45};
const set<int> PID_NEGATIVE = {11, -211, -321, -2212};
const float SPEED_OF_LIGHT = 29.9792458; // cm/ns
const float HTCC_PION_THRESHOLD = 4.9f;
const float LTCC_PION_THRESHOLD = 3.0f;

// Particle masses (in GeV)
map<int, float> PDG_MASS = {
    {11, 0.0005}, {-11, 0.0005}, {211, 0.13957018}, {-211, 0.13957018},
    {321, 0.49367716}, {-321, 0.49367716}, {2212, 0.938272046}, {-2212, 0.938272046},
    {45, 1.87705}
};

// Detector priorities (only for charged particles)
const vector<pair<int, vector<int>>> chargedBetaDetectors = {
    {FTOF_DETECTOR, {2, 1, 3}}, {CTOF_DETECTOR, {1}}, {ECAL_DETECTOR, {1, 4, 7}}
};

// Structures
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
    float beta;
    vector<DetectorHit> hits;
    bool is_trigger;
    int assigned_pid;
    float start_time;
    int status;
};

// Sampling Fraction Parameters
struct SamplingFractionParams {
    float sf[4];
    float sfs[4];
};
std::vector<SamplingFractionParams> sfParams(7);
float htccNpheCut = 2.0f;
float ltccNpheCut = 2.0f;

void loadCCDBParams() {
    ccdb::Calibration *calib = ccdb::CalibrationGenerator::CreateCalibration("mysql://clas12reader@clasdb.jlab.org/clas12", 18536, "default");
    if (!calib) { 
        cerr << "Failed to create CCDB Calibration object! Using defaults." << endl;
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
        htccNpheCut = 2.0f;
        ltccNpheCut = 2.0f;
        return;
    }

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
    delete calib;
}

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

int bestPidFromTiming(const ParticleData& p, float eventStartTime) {
    float startTime = eventStartTime;
    if (p.charge != 0 && p.vt != -9999.0f) startTime = p.vt;

    if (p.charge == 0) {
        return 0; // Skip neutrals
    } else {
        set<int> hypotheses = (p.charge > 0) ? PID_POSITIVE : PID_NEGATIVE;
        float min_dt = numeric_limits<float>::max();
        int best_pid = 0;
        for (const auto& det : chargedBetaDetectors) {
            for (int layer : det.second) {
                if (hasHit(p, det.first, layer)) {
                    for (int pid : hypotheses) {
                        if (abs(pid) == 11) continue; // Skip electrons/positrons for timing
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

void assignPids(vector<ParticleData>& particles, float eventStartTime, int eventNum) {
    for (size_t i = 0; i < particles.size(); ++i) {
        ParticleData& p = particles[i];
        p.assigned_pid = 0;
        bool is_elec = isSimpleElectron(p);
        int pid_from_timing = bestPidFromTiming(p, eventStartTime);
        bool htcc_signal = p.nphe_htcc > htccNpheCut;
        bool ltcc_signal = p.nphe_ltcc > ltccNpheCut;
        bool htcc_pion = p.p > HTCC_PION_THRESHOLD;
        bool ltcc_pion = p.p > LTCC_PION_THRESHOLD;

        set<int> hypotheses = (p.charge > 0) ? PID_POSITIVE : (p.charge < 0) ? PID_NEGATIVE : set<int>();
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
            }
            if (p.assigned_pid != 0) break;
        }
        // Print if assigned PID differs from REC::Particle PID or is unidentified
        if (p.assigned_pid != p.pid ) {
            cout << "Event " << eventNum << ", Particle " << i << ": REC::PID=" << p.pid 
                 << ", Assigned PID=" << p.assigned_pid << ", Charge=" << (int)p.charge 
                 << ", Momentum=" << p.p << ", NPHE_HTCC=" << p.nphe_htcc 
                 << ", NPHE_LTCC=" << p.nphe_ltcc << ", Status=" << p.status << endl;
        }
    }
}

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

int main() {
    TStopwatch timer;
    timer.Start();

    // Read directories from a file
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

    loadCCDBParams();

    int total_events_processed = 0;
    const int maxEvents = 100000;
    int event_count = 0;

    int total_electrons = 0;
    int total_positrons = 0;
    int total_positive_pions = 0;
    int total_negative_pions = 0;
    int total_positive_kaons = 0;
    int total_negative_kaons = 0;
    int total_protons = 0;
    int total_antiprotons = 0;
    int total_deuterons = 0;
    int total_unidentified = 0;

    int count_electrons = 0;
    int count_positrons = 0;
    int count_positive_pions = 0;
    int count_negative_pions = 0;
    int count_positive_kaons = 0;
    int count_negative_kaons = 0;
    int count_protons = 0;
    int count_antiprotons = 0;
    int count_deuterons = 0;
    int count_pidzero = 0;

    int directory_count = 0;
    int hipo_file_count = 0;

    // Structure to store mismatch details
    struct Mismatch {
        int eventNum;
        size_t particleIdx;
        int recPid;
        int assignedPid;
        int charge;
        float momentum;
        float npheHtcc;
        float npheLtcc;
        int status;
    };
    vector<Mismatch> mismatches;

    for (const auto& dir : directories) {
        directory_count++;
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

        for (const auto& file : hipoFiles) {
            hipo_file_count++;
            cout << "  Opening file: " << file << endl;
            hipo::reader reader;
            reader.open(file.c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);

            while (reader.next() && event_count < maxEvents) {
                event_count++;
                total_events_processed++;

                hipo::event event;
                reader.read(event);

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

                IndexMap cherMap = loadMapByIndex(CHER, "pindex");
                IndexMap caloMap = loadMapByIndex(CALO, "pindex");
                IndexMap scinMap = loadMapByIndex(SCIN, "pindex");

                float startTime = EVENT.getFloat("startTime", 0);

                // Process particles for PID
                vector<ParticleData> particles(PART.getRows());
                int valid_particles = 0;
                for (int i = 0; i < PART.getRows(); i++) {
                    ParticleData pd = getParticleData(i, PART, CHER, CALO, SCIN, cherMap, caloMap, scinMap);
                    if (pd.status == 0 || pd.charge == 0) continue; // Skip neutrals
                    particles[valid_particles++] = pd;
                }
                particles.resize(valid_particles);

                assignPids(particles, startTime, event_count);

                // Count all particles and collect mismatches
                for (size_t j = 0; j < particles.size(); ++j) {
                    const auto& p = particles[j];
                    int pid = p.pid;
                    int assigned_pid = p.assigned_pid;

                    // REC::Particle PID counts
                    if (pid == 11) total_electrons++;
                    else if (pid == -11) total_positrons++;
                    else if (pid == 211) total_positive_pions++;
                    else if (pid == -211) total_negative_pions++;
                    else if (pid == 321) total_positive_kaons++;
                    else if (pid == -321) total_negative_kaons++;
                    else if (pid == 2212) total_protons++;
                    else if (pid == -2212) total_antiprotons++;
                    else if (pid == 45) total_deuterons++;
                    else if (pid == 0) total_unidentified++;

                    // Assigned PID counts (all particles)
                    if (assigned_pid == 11) count_electrons++;
                    else if (assigned_pid == -11) count_positrons++;
                    else if (assigned_pid == 211) count_positive_pions++;
                    else if (assigned_pid == -211) count_negative_pions++;
                    else if (assigned_pid == 321) count_positive_kaons++;
                    else if (assigned_pid == -321) count_negative_kaons++;
                    else if (assigned_pid == 2212) count_protons++;
                    else if (assigned_pid == -2212) count_antiprotons++;
                    else if (assigned_pid == 45) count_deuterons++;
                    else if (assigned_pid == 0) count_pidzero++;

                    // Collect mismatch details
                    if (assigned_pid != pid ) {
                        mismatches.push_back({event_count, j, pid, assigned_pid, (int)p.charge, p.p, p.nphe_htcc, p.nphe_ltcc, p.status});
                    }
                }
            }
            if (event_count >= maxEvents) break;
        }
        if (event_count >= maxEvents) break;
    }

    timer.Stop();

    // Print results
    cout << "****************************************************************************************" << endl;
    cout << "Total number of events: " << event_count << endl;
    cout << "REC::Particle PID Counts:" << endl;
    cout << "Total electrons (11): " << total_electrons << endl;
    cout << "Total positrons (-11): " << total_positrons << endl;
    cout << "Total positive pions (211): " << total_positive_pions << endl;
    cout << "Total negative pions (-211): " << total_negative_pions << endl;
    cout << "Total positive kaons (321): " << total_positive_kaons << endl;
    cout << "Total negative kaons (-321): " << total_negative_kaons << endl;
    cout << "Total protons (2212): " << total_protons << endl;
    cout << "Total antiprotons (-2212): " << total_antiprotons << endl;
    cout << "Total deuterons (45): " << total_deuterons << endl;
    cout << "Total unidentified (0): " << total_unidentified << endl;
    cout << "Assigned PID Counts (after custom logic, all particles):" << endl;
    cout << "Total electrons (11): " << count_electrons << endl;
    cout << "Total positrons (-11): " << count_positrons << endl;
    cout << "Total positive pions (211): " << count_positive_pions << endl;
    cout << "Total negative pions (-211): " << count_negative_pions << endl;
    cout << "Total positive kaons (321): " << count_positive_kaons << endl;
    cout << "Total negative kaons (-321): " << count_negative_kaons << endl;
    cout << "Total protons (2212): " << count_protons << endl;
    cout << "Total antiprotons (-2212): " << count_antiprotons << endl;
    cout << "Total deuterons (45): " << count_deuterons << endl;
    cout << "Total unidentified (0): " << count_pidzero << endl;
    cout << "Total events processed: " << total_events_processed << endl;
    cout << "Real time: " << timer.RealTime() << " s, CPU time: " << timer.CpuTime() << " s" << endl;

    // Print mismatches separately
    if (!mismatches.empty()) {
        cout << "****************************************************************************************" << endl;
        cout << "Mismatch Details:" << endl;
        for (const auto& mismatch : mismatches) {
            cout << "Event " << mismatch.eventNum << ", Particle Index " << mismatch.particleIdx 
                 << ": REC::PID=" << mismatch.recPid << ", Assigned PID=" << mismatch.assignedPid 
                 << ", Charge=" << mismatch.charge << ", Momentum=" << mismatch.momentum 
                 << ", NPHE_HTCC=" << mismatch.npheHtcc << ", NPHE_LTCC=" << mismatch.npheLtcc 
                 << ", Status=" << mismatch.status << endl;
        }
        cout << "Total mismatches: " << mismatches.size() << endl;
        cout << "****************************************************************************************" << endl;
    } else {
        cout << "No mismatches found." << endl;
    }

    cout << "****************************************************************************************" << endl;

    return 0;
}