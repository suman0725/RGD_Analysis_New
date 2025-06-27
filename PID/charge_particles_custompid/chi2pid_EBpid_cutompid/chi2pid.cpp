#include <cstdlib>
#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <fstream>
#include <string>
#include <filesystem>
#include <set>
#include <tuple>
#include <cstddef>
#include <limits>
#include "reader.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include <TStopwatch.h>
#include <TLegend.h>
#include <TStyle.h>
#include "CCDB/Calibration.h"
#include "CCDB/CalibrationGenerator.h"

using namespace std;
namespace fs = std::filesystem;

// Define IndexMap for mapping particle indices to detector hits
using IndexMap = std::map<int, std::vector<int>>;

// Cache for FTOF timing resolutions
std::map<std::tuple<int, int, int>, float> tresCache;

// Constants for detector IDs and thresholds
const float MIN_PCAL_ENERGY = 0.06f; // Minimum PCAL energy for electron identification
const float NSIGMA_CUT = 5.0f;       // Sampling fraction cut for electrons
const int HTCC_DETECTOR = 15;        // High Threshold Cherenkov Counter
const int LTCC_DETECTOR = 16;        // Low Threshold Cherenkov Counter
const int ECAL_DETECTOR = 7;         // Electromagnetic Calorimeter
const int FTOF_DETECTOR = 12;        // Forward Time-of-Flight
const int CTOF_DETECTOR = 4;         // Central Time-of-Flight
const std::set<int> PID_POSITIVE = {-11, 211, 321, 2212, 45}; // Positive charged particle hypotheses
const std::set<int> PID_NEGATIVE = {11, -211, -321, -2212};    // Negative charged particle hypotheses
const float SPEED_OF_LIGHT = 29.9792458; // cm/ns
const float HTCC_PION_THRESHOLD = 4.9f;  // Momentum threshold for pions in HTCC (GeV)
const float LTCC_PION_THRESHOLD = 3.0f;  // Momentum threshold for pions in LTCC (GeV)

// Particle masses (in GeV) for theoretical beta calculation
std::map<int, float> PDG_MASS = {
    {11, 0.0005}, {-11, 0.0005}, {211, 0.13957018}, {-211, 0.13957018},
    {321, 0.49367716}, {-321, 0.49367716}, {2212, 0.938272046}, {-2212, 0.938272046},
    {45, 1.87705}
};

// Detector priorities for timing-based PID (only for charged particles)
const std::vector<std::pair<int, std::vector<int>>> chargedBetaDetectors = {
    {FTOF_DETECTOR, {2, 1, 3}}, {CTOF_DETECTOR, {1}}, {ECAL_DETECTOR, {1, 4, 7}}
};

// Structures
struct DetectorHit {
    int detector;
    int layer;
    int sector;    // Added for tres lookup
    int component; // Added for tres lookup
    float time;
    float path;
};

struct ParticleData {
    int charge;              // Particle charge from REC::Particle
    float p;                 // Momentum magnitude (GeV)
    int sector;              // Detector sector (from ECAL)
    float nphe_htcc;         // Number of photoelectrons in HTCC
    float nphe_ltcc;         // Number of photoelectrons in LTCC
    float energy_pcal;       // Energy in PCAL layer (GeV)
    float energy_total;      // Total energy in ECAL (GeV)
    int pid;                 // PID from REC::Particle
    float vt;                // Vertex time (ns) from REC::Particle
    float beta;              // Beta from REC::Particle (measured by TOF)
    std::vector<DetectorHit> hits; // Detector hits (FTOF, CTOF, ECAL)
    bool is_trigger;         // True if particle is a trigger particle (status < 0)
    int assigned_pid;        // Reassigned PID by custom logic
    float start_time;        // Start time for timing calculations (vt for charged particles)
    int status;              // Status from REC::Particle (detector collection passed)
    float chi2pid;           // Chi2pid from REC::Particle (set by COATJAVA)
    float chi2pid_custom;    // Computed PID quality for assigned PID
    float vz_tele;           // Vertex z-position from EVENT bank

    ParticleData() : charge(0), p(0.0f), sector(-1), nphe_htcc(0.0f), nphe_ltcc(0.0f),
                     energy_pcal(0.0f), energy_total(0.0f), pid(0), vt(0.0f), beta(-99.0f),
                     hits(), is_trigger(false), assigned_pid(0), start_time(0.0f), status(0),
                     chi2pid(0.0f), chi2pid_custom(99999.0f), vz_tele(0.0f) {}
};

// Sampling Fraction Parameters for electron identification
struct SamplingFractionParams {
    float sf[4];  // Mean sampling fraction parameters
    float sfs[4]; // Sigma sampling fraction parameters
};
std::vector<SamplingFractionParams> sfParams(7);
float htccNpheCut = 2.0f;
float ltccNpheCut = 2.0f;

// Enum to define regions
enum class Region { Forward, Central, Band, Unknown };

struct ParticleStats {
    int total_particles_before_filtering = 0;
    int total_unidentified_before_filtering = 0;
    int total_filtered_out = 0;
    int total_filtered_out_pid_zero = 0;
    int total_invalid_chi2pid = 0;
    int total_invalid_chi2pid_nonzero_pid = 0;
    int total_invalid_chi2pid_custom = 0;
    int total_invalid_chi2pid_custom_nonzero_pid = 0;
    int valid_chi2pid_custom = 0;
    int valid_chi2pid_custom_pid_zero = 0;
    int invalid_chi2pid_custom_pid_zero = 0;
    int valid_chi2pid_custom_nonzero_pid = 0;
    int valid_chi2pid = 0;

    // Maps for PID and assigned PID counts per region
    std::map<int, int> pid_counts;                        // Overall PID counts
    std::map<int, int> assigned_pid_counts;               // Overall assigned PID counts
    std::map<Region, std::map<int, int>> pid_counts_by_region;    // PID counts by region
    std::map<Region, std::map<int, int>> assigned_pid_counts_by_region; // Assigned PID counts by region
};

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
    float chi2pid;
    float chi2pidCustom;
    bool hasFtof;
    bool hasCtof;
};

// Function declarations
void loadCCDBParams();
ParticleData getParticleData(int partIdx, hipo::bank& PART, hipo::bank& CHER, hipo::bank& CALO, hipo::bank& SCIN,
                            IndexMap& cherMap, IndexMap& caloMap, IndexMap& scinMap);
IndexMap loadMapByIndex(hipo::bank& fromBank, const char* idxVarName);
bool isSimpleElectron(const ParticleData& pd);
int bestPidFromTiming(const ParticleData& p, float eventStartTime, int eventNum, size_t particleIdx);
void assignPids(std::vector<ParticleData>& particles, float eventStartTime, int eventNum);
float PIDQuality(const ParticleData& p, int pid);
float computeDeltaT(const ParticleData& p, int pid);
float getSamplingFractionMean(int sector, float measuredEnergy);
float getSamplingFractionSigma(int sector, float measuredEnergy);
float getSamplingFractionNSigma(float samplingFraction, float mean, float sigma);
float getTheoryBeta(float p, float mass);
bool hasHit(const ParticleData& p, int detector, int layer = -1);
Region determineRegion(int status);
void validateAndCompare(std::vector<ParticleData>& particles, int eventNum, ParticleStats& statsBefore, ParticleStats& statsAfter,
                       std::vector<Mismatch>& mismatchesBefore, std::vector<Mismatch>& mismatchesAfter,
                       std::vector<float>& chi2pidListBefore, std::vector<float>& chi2pidListAfter,
                       std::vector<float>& chi2pidCustomListBefore, std::vector<float>& chi2pidCustomListAfter,
                       std::map<int, TH1F*>& hMomentumRecBefore, std::map<int, TH1F*>& hMomentumCustomBefore,
                       std::map<int, TH1F*>& hMomentumRecAfter, std::map<int, TH1F*>& hMomentumCustomAfter,
                       std::map<int, TH1F*>& hChi2pidRecBefore, std::map<int, TH1F*>& hChi2pidCustomBefore,
                       std::map<int, TH1F*>& hChi2pidRecAfter, std::map<int, TH1F*>& hChi2pidCustomAfter,
                       std::map<int, TH1F*>& hDtRecBefore, std::map<int, TH1F*>& hDtCustomBefore,
                       std::map<int, TH1F*>& hDtRecAfter, std::map<int, TH1F*>& hDtCustomAfter,
                       std::map<int, TH2F*>& hDtVsPRecBefore, std::map<int, TH2F*>& hDtVsPCustomBefore,
                       std::map<int, TH2F*>& hDtVsPRecAfter, std::map<int, TH2F*>& hDtVsPCustomAfter,
                       std::map<int, TH2F*>& hChi2pidVsPRecBefore, std::map<int, TH2F*>& hChi2pidVsPCustomBefore,
                       std::map<int, TH2F*>& hChi2pidVsPRecAfter, std::map<int, TH2F*>& hChi2pidVsPCustomAfter,
                       std::map<int, TH2F*>& hBetaVsPRecBefore, std::map<int, TH2F*>& hBetaVsPCustomBefore,
                       std::map<int, TH2F*>& hBetaVsPRecAfter, std::map<int, TH2F*>& hBetaVsPCustomAfter,
                       bool hasTriggerElectronRec, bool hasTriggerElectronCustom);
void saveStatsToCSV(const ParticleStats& stats, const string& filename, const string& prefix);
void savePlots(const string& dir, const string& suffix,
               std::vector<float>& chi2pidList, std::vector<float>& chi2pidCustomList,
               std::map<int, TH1F*>& hMomentumRec, std::map<int, TH1F*>& hMomentumCustom,
               std::map<int, TH1F*>& hChi2pidRec, std::map<int, TH1F*>& hChi2pidCustom,
               std::map<int, TH1F*>& hDtRec, std::map<int, TH1F*>& hDtCustom,
               std::map<int, TH2F*>& hDtVsPRec, std::map<int, TH2F*>& hDtVsPCustom,
               std::map<int, TH2F*>& hChi2pidVsPRec, std::map<int, TH2F*>& hChi2pidVsPCustom,
               std::map<int, TH2F*>& hBetaVsPRec, std::map<int, TH2F*>& hBetaVsPCustom);

// Function implementations
void loadCCDBParams() {
    ccdb::Calibration *calib = ccdb::CalibrationGenerator::CreateCalibration(
        "mysql://clas12reader@clasdb.jlab.org/clas12", 18536, "default");
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
    if (sigma == 0) return 99999.0f;
    return (samplingFraction - mean) / sigma;
}

float getTheoryBeta(float p, float mass) {
    return p / sqrt(p * p + mass * mass);
}

bool hasHit(const ParticleData& p, int detector, int layer) {
    for (const auto& hit : p.hits) {
        if (hit.detector == detector && (layer == -1 || hit.layer == layer)) {
            return true;
        }
    }
    return false;
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
    pd.is_trigger = (pd.status < 0);
    pd.start_time = pd.vt;
    pd.vz_tele = PART.getFloat("vz", partIdx);

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

int bestPidFromTiming(const ParticleData& p, float eventStartTime, int eventNum, size_t particleIdx) {
    float startTime = eventStartTime;
    if (p.charge != 0 && p.vt != -9999.0f) startTime = p.vt;

    if (p.charge == 0) {
        return 0;
    } else {
        std::set<int> hypotheses = (p.charge > 0) ? PID_POSITIVE : PID_NEGATIVE;
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
                                float dt = vt - startTime;
                                if (abs(dt) < min_dt) {
                                    min_dt = abs(dt);
                                    best_pid = pid;
                                    cout << "Timing: Event " << eventNum << ", Particle " << particleIdx 
                                         << ", Hypothesis PID=" << pid << ", dt=" << dt 
                                         << ", beta_theory=" << beta_theory << endl;
                                }
                            }
                        }
                    }
                    if (best_pid != 0) return best_pid;
                }
            }
        }
        return best_pid;
    }
}

void assignPids(std::vector<ParticleData>& particles, float eventStartTime, int eventNum) {
    for (size_t i = 0; i < particles.size(); ++i) {
        ParticleData& p = particles[i];
        p.assigned_pid = 0;
        if (p.charge == 0) continue;

        bool is_elec = isSimpleElectron(p);
        int pid_from_timing = bestPidFromTiming(p, eventStartTime, eventNum, i);
        bool htcc_signal = p.nphe_htcc > htccNpheCut;
        bool ltcc_signal = p.nphe_ltcc > ltccNpheCut;
        bool htcc_pion = p.p > HTCC_PION_THRESHOLD;
        bool ltcc_pion = p.p > LTCC_PION_THRESHOLD;

        std::set<int> hypotheses = (p.charge > 0) ? PID_POSITIVE : PID_NEGATIVE;
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

        if (p.assigned_pid != 0) {
            p.chi2pid_custom = PIDQuality(p, p.assigned_pid);
        }
    }
}

float PIDQuality(const ParticleData& p, int pid) {
    float q = 99999.0f;

    if (abs(pid) == 11) {
        float sf = p.energy_total / p.p;
        float mean = getSamplingFractionMean(p.sector, p.energy_total);
        float sigma = getSamplingFractionSigma(p.sector, p.energy_total);
        q = getSamplingFractionNSigma(sf, mean, sigma);
    } else if (p.charge != 0) {
        float sigma = -1.0f;
        float delta_t = 99999.0f;
        bool found = false;
        for (const auto& det : chargedBetaDetectors) {
            for (int layer : det.second) {
                if (hasHit(p, det.first, layer)) {
                    for (const auto& hit : p.hits) {
                        if (hit.detector == det.first && hit.layer == layer) {
                            if (hit.detector == FTOF_DETECTOR) {
                                if (tresCache.count(make_tuple(hit.sector, hit.layer, hit.component))) {
                                    sigma = tresCache[make_tuple(hit.sector, hit.layer, hit.component)];
                                } else {
                                    return 99999.0f;
                                }
                            } else if (hit.detector == CTOF_DETECTOR) {
                                sigma = 0.065f;
                            } else {
                                return 99999.0f;
                            }
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
        if (sigma > 0) q = delta_t / sigma;
        else q = 99999.0f;
    }
    return q;
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

Region determineRegion(int status) {
    int abs_status = std::abs(status);
    if (abs_status / 2000 == 1) return Region::Forward;
    else if (abs_status / 4000 == 1) return Region::Central;
    else if (abs_status / 8000 == 1) return Region::Band;
    else {
        std::cout << "Unknown status: " << abs_status << std::endl;
        return Region::Unknown;
    }
}

void validateAndCompare(std::vector<ParticleData>& particles, int eventNum, ParticleStats& statsBefore, ParticleStats& statsAfter,
                       std::vector<Mismatch>& mismatchesBefore, std::vector<Mismatch>& mismatchesAfter,
                       std::vector<float>& chi2pidListBefore, std::vector<float>& chi2pidListAfter,
                       std::vector<float>& chi2pidCustomListBefore, std::vector<float>& chi2pidCustomListAfter,
                       std::map<int, TH1F*>& hMomentumRecBefore, std::map<int, TH1F*>& hMomentumCustomBefore,
                       std::map<int, TH1F*>& hMomentumRecAfter, std::map<int, TH1F*>& hMomentumCustomAfter,
                       std::map<int, TH1F*>& hChi2pidRecBefore, std::map<int, TH1F*>& hChi2pidCustomBefore,
                       std::map<int, TH1F*>& hChi2pidRecAfter, std::map<int, TH1F*>& hChi2pidCustomAfter,
                       std::map<int, TH1F*>& hDtRecBefore, std::map<int, TH1F*>& hDtCustomBefore,
                       std::map<int, TH1F*>& hDtRecAfter, std::map<int, TH1F*>& hDtCustomAfter,
                       std::map<int, TH2F*>& hDtVsPRecBefore, std::map<int, TH2F*>& hDtVsPCustomBefore,
                       std::map<int, TH2F*>& hDtVsPRecAfter, std::map<int, TH2F*>& hDtVsPCustomAfter,
                       std::map<int, TH2F*>& hChi2pidVsPRecBefore, std::map<int, TH2F*>& hChi2pidVsPCustomBefore,
                       std::map<int, TH2F*>& hChi2pidVsPRecAfter, std::map<int, TH2F*>& hChi2pidVsPCustomAfter,
                       std::map<int, TH2F*>& hBetaVsPRecBefore, std::map<int, TH2F*>& hBetaVsPCustomBefore,
                       std::map<int, TH2F*>& hBetaVsPRecAfter, std::map<int, TH2F*>& hBetaVsPCustomAfter,
                       bool hasTriggerElectronRec, bool hasTriggerElectronCustom) {
    int valid_chi2pid_count_before = 0, valid_chi2pid_count_after = 0;
    int valid_chi2pid_custom_count_before = 0, valid_chi2pid_custom_count_after = 0;
    int charged_particles_before = 0, charged_particles_after = 0;

    for (size_t j = 0; j < particles.size(); ++j) {
        const auto& p = particles[j];

        // Skip neutral particles for charged-particle-specific histograms
        if (p.charge == 0) continue;

        // For "before_trigger", process all particles regardless of trigger electron presence
        {
            ParticleStats& currentStats = statsBefore;
            std::vector<Mismatch>& currentMismatches = mismatchesBefore;
            std::vector<float>& currentChi2pidList = chi2pidListBefore;
            std::vector<float>& currentChi2pidCustomList = chi2pidCustomListBefore;
            std::map<int, TH1F*>& currentMomentumRec = hMomentumRecBefore;
            std::map<int, TH1F*>& currentMomentumCustom = hMomentumCustomBefore;
            std::map<int, TH1F*>& currentChi2pidRec = hChi2pidRecBefore;
            std::map<int, TH1F*>& currentChi2pidCustom = hChi2pidCustomBefore;
            std::map<int, TH1F*>& currentDtRec = hDtRecBefore;
            std::map<int, TH1F*>& currentDtCustom = hDtCustomBefore;
            std::map<int, TH2F*>& currentDtVsPRec = hDtVsPRecBefore;
            std::map<int, TH2F*>& currentDtVsPCustom = hDtVsPCustomBefore;
            std::map<int, TH2F*>& currentChi2pidVsPRec = hChi2pidVsPRecBefore;
            std::map<int, TH2F*>& currentChi2pidVsPCustom = hChi2pidVsPCustomBefore;
            std::map<int, TH2F*>& currentBetaVsPRec = hBetaVsPRecBefore;
            std::map<int, TH2F*>& currentBetaVsPCustom = hBetaVsPCustomBefore;

            charged_particles_before++;

            Region region = determineRegion(p.status);
            currentStats.pid_counts[p.pid]++;
            currentStats.pid_counts_by_region[region][p.pid]++;
            currentStats.assigned_pid_counts[p.assigned_pid]++;
            currentStats.assigned_pid_counts_by_region[region][p.assigned_pid]++;

            if (currentMomentumRec.count(p.pid)) {
                currentMomentumRec[p.pid]->Fill(p.p);
            }
            if (currentMomentumCustom.count(p.assigned_pid)) {
                currentMomentumCustom[p.assigned_pid]->Fill(p.p);
            }

            cout << "Event " << eventNum << ", Particle " << j 
                 << " (Before Trigger): PID=" << p.pid << ", Assigned PID=" << p.assigned_pid 
                 << ", Momentum=" << p.p << ", Measured Beta=" << p.beta 
                 << ", Status=" << p.status << endl;

            if (p.p > 0 && p.beta >= 0 && p.beta <= 1.2) {
                if (currentBetaVsPRec.count(p.pid)) {
                    currentBetaVsPRec[p.pid]->Fill(p.p, p.beta);
                }
                if (p.assigned_pid != 0 && currentBetaVsPCustom.count(p.assigned_pid)) {
                    currentBetaVsPCustom[p.assigned_pid]->Fill(p.p, p.beta);
                }
            }

            if (p.chi2pid != 9999.0f) {
                valid_chi2pid_count_before++;
                statsBefore.valid_chi2pid++;
                currentChi2pidList.push_back(p.chi2pid);
                if (currentChi2pidRec.count(p.pid)) {
                    currentChi2pidRec[p.pid]->Fill(p.chi2pid);
                }
                if (currentChi2pidVsPRec.count(p.pid)) {
                    currentChi2pidVsPRec[p.pid]->Fill(p.p, p.chi2pid);
                }
            } else {
                currentStats.total_invalid_chi2pid++;
                if (p.pid != 0) currentStats.total_invalid_chi2pid_nonzero_pid++;
            }

            if (p.chi2pid_custom != 99999.0f) {
                valid_chi2pid_custom_count_before++;
                statsBefore.valid_chi2pid_custom++;
                if (p.assigned_pid == 0) statsBefore.valid_chi2pid_custom_pid_zero++;
                else statsBefore.valid_chi2pid_custom_nonzero_pid++;
                currentChi2pidCustomList.push_back(p.chi2pid_custom);
                if (currentChi2pidCustom.count(p.assigned_pid)) {
                    currentChi2pidCustom[p.assigned_pid]->Fill(p.chi2pid_custom);
                }
                if (currentChi2pidVsPCustom.count(p.assigned_pid)) {
                    currentChi2pidVsPCustom[p.assigned_pid]->Fill(p.p, p.chi2pid_custom);
                }
            } else {
                currentStats.total_invalid_chi2pid_custom++;
                if (p.assigned_pid != 0) {
                    currentStats.total_invalid_chi2pid_custom_nonzero_pid++;
                } else {
                    currentStats.invalid_chi2pid_custom_pid_zero++;
                }
            }

            float dt_rec = computeDeltaT(p, p.pid);
            if (dt_rec != 99999.0f && currentDtRec.count(p.pid)) {
                currentDtRec[p.pid]->Fill(dt_rec);
                if (currentDtVsPRec.count(p.pid)) {
                    currentDtVsPRec[p.pid]->Fill(p.p, dt_rec);
                }
            }

            float dt_custom = computeDeltaT(p, p.assigned_pid);
            if (dt_custom != 99999.0f && currentDtCustom.count(p.assigned_pid)) {
                currentDtCustom[p.assigned_pid]->Fill(dt_custom);
                if (currentDtVsPCustom.count(p.assigned_pid)) {
                    currentDtVsPCustom[p.assigned_pid]->Fill(p.p, dt_custom);
                }
            }

            if (p.assigned_pid != p.pid) {
                bool has_ftof = hasHit(p, FTOF_DETECTOR);
                bool has_ctof = hasHit(p, CTOF_DETECTOR);
                currentMismatches.push_back({eventNum, j, p.pid, p.assigned_pid, (int)p.charge, p.p,
                                             p.nphe_htcc, p.nphe_ltcc, p.status, p.chi2pid,
                                             p.chi2pid_custom, has_ftof, has_ctof});
            }
        }

        // For "after_trigger", only process if the event has a trigger electron
        if (hasTriggerElectronRec || hasTriggerElectronCustom) {
            ParticleStats& currentStats = statsAfter;
            std::vector<Mismatch>& currentMismatches = mismatchesAfter;
            std::vector<float>& currentChi2pidList = chi2pidListAfter;
            std::vector<float>& currentChi2pidCustomList = chi2pidCustomListAfter;
            std::map<int, TH1F*>& currentMomentumRec = hMomentumRecAfter;
            std::map<int, TH1F*>& currentMomentumCustom = hMomentumCustomAfter;
            std::map<int, TH1F*>& currentChi2pidRec = hChi2pidRecAfter;
            std::map<int, TH1F*>& currentChi2pidCustom = hChi2pidCustomAfter;
            std::map<int, TH1F*>& currentDtRec = hDtRecAfter;
            std::map<int, TH1F*>& currentDtCustom = hDtCustomAfter;
            std::map<int, TH2F*>& currentDtVsPRec = hDtVsPRecAfter;
            std::map<int, TH2F*>& currentDtVsPCustom = hDtVsPCustomAfter;
            std::map<int, TH2F*>& currentChi2pidVsPRec = hChi2pidVsPRecAfter;
            std::map<int, TH2F*>& currentChi2pidVsPCustom = hChi2pidVsPCustomAfter;
            std::map<int, TH2F*>& currentBetaVsPRec = hBetaVsPRecAfter;
            std::map<int, TH2F*>& currentBetaVsPCustom = hBetaVsPCustomAfter;

            charged_particles_after++;

            Region region = determineRegion(p.status);
            currentStats.pid_counts[p.pid]++;
            currentStats.pid_counts_by_region[region][p.pid]++;
            currentStats.assigned_pid_counts[p.assigned_pid]++;
            currentStats.assigned_pid_counts_by_region[region][p.assigned_pid]++;

            if (currentMomentumRec.count(p.pid)) {
                currentMomentumRec[p.pid]->Fill(p.p);
            }
            if (currentMomentumCustom.count(p.assigned_pid)) {
                currentMomentumCustom[p.assigned_pid]->Fill(p.p);
            }

            cout << "Event " << eventNum << ", Particle " << j 
                 << " (After Trigger): PID=" << p.pid << ", Assigned PID=" << p.assigned_pid 
                 << ", Momentum=" << p.p << ", Measured Beta=" << p.beta 
                 << ", Status=" << p.status << endl;

            if (p.p > 0 && p.beta >= 0 && p.beta <= 1.2) {
                if (currentBetaVsPRec.count(p.pid)) {
                    currentBetaVsPRec[p.pid]->Fill(p.p, p.beta);
                }
                if (p.assigned_pid != 0 && currentBetaVsPCustom.count(p.assigned_pid)) {
                    currentBetaVsPCustom[p.assigned_pid]->Fill(p.p, p.beta);
                }
            }

            if (p.chi2pid != 9999.0f) {
                valid_chi2pid_count_after++;
                statsAfter.valid_chi2pid++;
                currentChi2pidList.push_back(p.chi2pid);
                if (currentChi2pidRec.count(p.pid)) {
                    currentChi2pidRec[p.pid]->Fill(p.chi2pid);
                }
                if (currentChi2pidVsPRec.count(p.pid)) {
                    currentChi2pidVsPRec[p.pid]->Fill(p.p, p.chi2pid);
                }
            } else {
                currentStats.total_invalid_chi2pid++;
                if (p.pid != 0) currentStats.total_invalid_chi2pid_nonzero_pid++;
            }

            if (p.chi2pid_custom != 99999.0f) {
                valid_chi2pid_custom_count_after++;
                statsAfter.valid_chi2pid_custom++;
                if (p.assigned_pid == 0) statsAfter.valid_chi2pid_custom_pid_zero++;
                else statsAfter.valid_chi2pid_custom_nonzero_pid++;
                currentChi2pidCustomList.push_back(p.chi2pid_custom);
                if (currentChi2pidCustom.count(p.assigned_pid)) {
                    currentChi2pidCustom[p.assigned_pid]->Fill(p.chi2pid_custom);
                }
                if (currentChi2pidVsPCustom.count(p.assigned_pid)) {
                    currentChi2pidVsPCustom[p.assigned_pid]->Fill(p.p, p.chi2pid_custom);
                }
            } else {
                currentStats.total_invalid_chi2pid_custom++;
                if (p.assigned_pid != 0) {
                    currentStats.total_invalid_chi2pid_custom_nonzero_pid++;
                } else {
                    currentStats.invalid_chi2pid_custom_pid_zero++;
                }
            }

            float dt_rec = computeDeltaT(p, p.pid);
            if (dt_rec != 99999.0f && currentDtRec.count(p.pid)) {
                currentDtRec[p.pid]->Fill(dt_rec);
                if (currentDtVsPRec.count(p.pid)) {
                    currentDtVsPRec[p.pid]->Fill(p.p, dt_rec);
                }
            }

            float dt_custom = computeDeltaT(p, p.assigned_pid);
            if (dt_custom != 99999.0f && currentDtCustom.count(p.assigned_pid)) {
                currentDtCustom[p.assigned_pid]->Fill(dt_custom);
                if (currentDtVsPCustom.count(p.assigned_pid)) {
                    currentDtVsPCustom[p.assigned_pid]->Fill(p.p, dt_custom);
                }
            }

            if (p.assigned_pid != p.pid) {
                bool has_ftof = hasHit(p, FTOF_DETECTOR);
                bool has_ctof = hasHit(p, CTOF_DETECTOR);
                currentMismatches.push_back({eventNum, j, p.pid, p.assigned_pid, (int)p.charge, p.p,
                                             p.nphe_htcc, p.nphe_ltcc, p.status, p.chi2pid,
                                             p.chi2pid_custom, has_ftof, has_ctof});
            }
        }
    }

    cout << "Event " << eventNum << ": Charged particles (Before Trigger)=" << charged_particles_before
         << ", Charged particles (After Trigger)=" << charged_particles_after
         << ", Valid chi2pid (Before)=" << valid_chi2pid_count_before
         << ", Valid chi2pid (After)=" << valid_chi2pid_count_after
         << ", Valid chi2pid_custom (Before)=" << valid_chi2pid_custom_count_before
         << ", Valid chi2pid_custom (After)=" << valid_chi2pid_custom_count_after << endl;
}

void saveStatsToCSV(const ParticleStats& stats, const string& filename, const string& prefix) {
    ofstream csvFile(filename);
    if (!csvFile.is_open()) {
        cerr << "Error: Could not open " << filename << " for writing!" << endl;
        return;
    }

    csvFile << "Statistic,Value\n";
    csvFile << prefix << "Total particles before filtering," << stats.total_particles_before_filtering << "\n";
    csvFile << prefix << "Total unidentified before filtering (pid == 0)," << stats.total_unidentified_before_filtering << "\n";
    csvFile << prefix << "Total particles filtered out (status == 0 || charge == 0)," << stats.total_filtered_out << "\n";
    csvFile << prefix << "Total particles filtered out with pid == 0," << stats.total_filtered_out_pid_zero << "\n";
    csvFile << prefix << "Total invalid chi2pid," << stats.total_invalid_chi2pid << "\n";
    csvFile << prefix << "Total invalid chi2pid nonzero pid," << stats.total_invalid_chi2pid_nonzero_pid << "\n";
    csvFile << prefix << "Total invalid chi2pid custom," << stats.total_invalid_chi2pid_custom << "\n";
    csvFile << prefix << "Total invalid chi2pid custom nonzero pid," << stats.total_invalid_chi2pid_custom_nonzero_pid << "\n";
    csvFile << prefix << "Valid chi2pid custom," << stats.valid_chi2pid_custom << "\n";
    csvFile << prefix << "Valid chi2pid custom pid zero," << stats.valid_chi2pid_custom_pid_zero << "\n";
    csvFile << prefix << "Invalid chi2pid custom pid zero," << stats.invalid_chi2pid_custom_pid_zero << "\n";
    csvFile << prefix << "Valid chi2pid custom nonzero pid," << stats.valid_chi2pid_custom_nonzero_pid << "\n";
    csvFile << prefix << "Valid chi2pid," << stats.valid_chi2pid << "\n";

    csvFile << "\nPID Counts (Overall)\n";
    csvFile << "pid,Pid from REC,Custom Pid\n";
    vector<int> particle_types = {11, -11, 211, -211, 321, -321, 2212, -2212, 45};
    for (int pid : particle_types) {
        int rec_count = stats.pid_counts.count(pid) ? stats.pid_counts.at(pid) : 0;
        int custom_count = stats.assigned_pid_counts.count(pid) ? stats.assigned_pid_counts.at(pid) : 0;
        csvFile << pid << "," << rec_count << "," << custom_count << "\n";
    }
    int rec_count_pid0 = stats.pid_counts.count(0) ? stats.pid_counts.at(0) : 0;
    int custom_count_pid0 = stats.assigned_pid_counts.count(0) ? stats.assigned_pid_counts.at(0) : 0;
    csvFile << "0," << rec_count_pid0 << "," << custom_count_pid0 << "\n";

    const std::map<Region, std::string> region_names = {
        {Region::Forward, "Forward"}, {Region::Central, "Central"}, {Region::Band, "Band"}, {Region::Unknown, "Unknown"}
    };
    for (const auto& region_pair : region_names) {
        csvFile << "\nPID Counts (" << region_pair.second << ")\n";
        csvFile << "pid,Pid from REC,Custom Pid\n";
        for (int pid : particle_types) {
            int rec_region_count = stats.pid_counts_by_region.count(region_pair.first) && 
                                  stats.pid_counts_by_region.at(region_pair.first).count(pid) ?
                                  stats.pid_counts_by_region.at(region_pair.first).at(pid) : 0;
            int custom_region_count = stats.assigned_pid_counts_by_region.count(region_pair.first) && 
                                     stats.assigned_pid_counts_by_region.at(region_pair.first).count(pid) ?
                                     stats.assigned_pid_counts_by_region.at(region_pair.first).at(pid) : 0;
            csvFile << pid << "," << rec_region_count << "," << custom_region_count << "\n";
        }
        int rec_region_count_pid0 = stats.pid_counts_by_region.count(region_pair.first) && 
                                   stats.pid_counts_by_region.at(region_pair.first).count(0) ?
                                   stats.pid_counts_by_region.at(region_pair.first).at(0) : 0;
        int custom_region_count_pid0 = stats.assigned_pid_counts_by_region.count(region_pair.first) && 
                                      stats.assigned_pid_counts_by_region.at(region_pair.first).count(0) ?
                                      stats.assigned_pid_counts_by_region.at(region_pair.first).at(0) : 0;
        csvFile << "0," << rec_region_count_pid0 << "," << custom_region_count_pid0 << "\n";
    }

    csvFile.close();
    cout << "Saved statistics to " << filename << endl;
}

void savePlots(const string& dir, const string& suffix,
               std::vector<float>& chi2pidList, std::vector<float>& chi2pidCustomList,
               std::map<int, TH1F*>& hMomentumRec, std::map<int, TH1F*>& hMomentumCustom,
               std::map<int, TH1F*>& hChi2pidRec, std::map<int, TH1F*>& hChi2pidCustom,
               std::map<int, TH1F*>& hDtRec, std::map<int, TH1F*>& hDtCustom,
               std::map<int, TH2F*>& hDtVsPRec, std::map<int, TH2F*>& hDtVsPCustom,
               std::map<int, TH2F*>& hChi2pidVsPRec, std::map<int, TH2F*>& hChi2pidVsPCustom,
               std::map<int, TH2F*>& hBetaVsPRec, std::map<int, TH2F*>& hBetaVsPCustom) {
    vector<int> particle_types = {11, -11, 211, -211, 321, -321, 2212, -2212, 45};

    // Chi2pid Comparison Plot (Combined for all charged particles)
    TCanvas *canvasChi2pid = new TCanvas(("canvasChi2pid_" + suffix).c_str(), ("Chi2pid Comparison for Charged Particles " + suffix).c_str(), 1200, 600);
    canvasChi2pid->Divide(3, 1);
    TH1F *hChi2pidAll = new TH1F(("hChi2pidAll_" + suffix).c_str(), ("Chi2pid (REC::Particle) " + suffix + ";Chi2pid;Counts").c_str(), 200, -15, 15);
    TH1F *hChi2pidCustomAll = new TH1F(("hChi2pidCustomAll_" + suffix).c_str(), ("Chi2pid_custom (Assigned) " + suffix + ";Chi2pid;Counts").c_str(), 200, -15, 15);
    TH1F *hChi2pidDiff = new TH1F(("hChi2pidDiff_" + suffix).c_str(), ("Chi2pid_custom - Chi2pid (Charged Particles) " + suffix + ";Difference;Counts").c_str(), 200, -15, 15);
    hChi2pidAll->SetLineColor(kBlue);
    hChi2pidAll->SetFillColor(kBlue);
    hChi2pidAll->SetFillStyle(3004);
    hChi2pidCustomAll->SetLineColor(kRed);
    hChi2pidCustomAll->SetFillColor(kRed);
    hChi2pidCustomAll->SetFillStyle(3005);
    hChi2pidDiff->SetLineColor(kGreen);
    hChi2pidDiff->SetFillColor(kGreen);
    hChi2pidDiff->SetFillStyle(3006);

    for (size_t i = 0; i < chi2pidList.size() && i < chi2pidCustomList.size(); ++i) {
        hChi2pidAll->Fill(chi2pidList[i]);
        hChi2pidCustomAll->Fill(chi2pidCustomList[i]);
        hChi2pidDiff->Fill(chi2pidCustomList[i] - chi2pidList[i]);
    }

    canvasChi2pid->cd(1);
    hChi2pidAll->Draw("HIST");
    canvasChi2pid->cd(2);
    hChi2pidCustomAll->Draw("HIST");
    canvasChi2pid->cd(3);
    hChi2pidAll->Draw("HIST");
    hChi2pidCustomAll->Draw("HIST SAME");
    hChi2pidDiff->Draw("HIST SAME");
    TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg1->AddEntry(hChi2pidAll, "REC::Particle Chi2pid", "f");
    leg1->AddEntry(hChi2pidCustomAll, "Chi2pid_custom (Assigned PID)", "f");
    leg1->AddEntry(hChi2pidDiff, "Chi2pid_custom - Chi2pid", "f");
    leg1->Draw();

    string outputChi2pid = dir + "/output_" + suffix + ".pdf";
    canvasChi2pid->Print(outputChi2pid.c_str());
    cout << "Chi2pid histograms saved to " << outputChi2pid << endl;

    // Momentum Plots (Combined)
    TCanvas *canvasMomentum = new TCanvas(("canvasMomentum_" + suffix).c_str(), ("Momentum Comparison for Charged Particles " + suffix).c_str(), 1200, 1200);
    canvasMomentum->Divide(3, 3);
    int pad = 1;
    for (int pid : particle_types) {
        canvasMomentum->cd(pad++);
        hMomentumRec[pid]->Draw("HIST");
        hMomentumCustom[pid]->Draw("HIST SAME");
        TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        string rec_label = "REC::PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hMomentumRec[pid]->GetEntries())) + ")";
        string custom_label = "Assigned PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hMomentumCustom[pid]->GetEntries())) + ")";
        leg->AddEntry(hMomentumRec[pid], rec_label.c_str(), "f");
        leg->AddEntry(hMomentumCustom[pid], custom_label.c_str(), "f");
        leg->Draw();
    }

    string outputMomentum = dir + "/momentum_plots_" + suffix + ".pdf";
    canvasMomentum->Print(outputMomentum.c_str());
    cout << "Momentum plots saved to " << outputMomentum << endl;

    // Chi2pid Plots (Combined)
    TCanvas *canvasChi2pidPerPid = new TCanvas(("canvasChi2pidPerPid_" + suffix).c_str(), ("Chi2pid Comparison Per Particle Type " + suffix).c_str(), 1200, 1200);
    canvasChi2pidPerPid->Divide(3, 3);
    pad = 1;
    for (int pid : particle_types) {
        canvasChi2pidPerPid->cd(pad++);
        hChi2pidRec[pid]->Draw("HIST");
        hChi2pidCustom[pid]->Draw("HIST SAME");
        TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        string rec_label = "REC::PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hChi2pidRec[pid]->GetEntries())) + ")";
        string custom_label = "Assigned PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hChi2pidCustom[pid]->GetEntries())) + ")";
        leg->AddEntry(hChi2pidRec[pid], rec_label.c_str(), "f");
        leg->AddEntry(hChi2pidCustom[pid], custom_label.c_str(), "f");
        leg->Draw();
    }

    string outputChi2pidPerPid = dir + "/chi2pid_plots_" + suffix + ".pdf";
    canvasChi2pidPerPid->Print(outputChi2pidPerPid.c_str());
    cout << "Chi2pid per particle plots saved to " << outputChi2pidPerPid << endl;

    // Delta T Plots (Combined)
    TCanvas *canvasDt = new TCanvas(("canvasDt_" + suffix).c_str(), ("Delta T Comparison for Charged Particles " + suffix).c_str(), 1200, 1200);
    canvasDt->Divide(3, 3);
    pad = 1;
    for (int pid : particle_types) {
        canvasDt->cd(pad++);
        hDtRec[pid]->Draw("HIST");
        hDtCustom[pid]->Draw("HIST SAME");
        TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        string rec_label = "REC::PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hDtRec[pid]->GetEntries())) + ")";
        string custom_label = "Assigned PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hDtCustom[pid]->GetEntries())) + ")";
        leg->AddEntry(hDtRec[pid], rec_label.c_str(), "f");
        leg->AddEntry(hDtCustom[pid], custom_label.c_str(), "f");
        leg->Draw();
    }

    string outputDt = dir + "/delta_t_plots_" + suffix + ".pdf";
    canvasDt->Print(outputDt.c_str());
    cout << "Delta T plots saved to " << outputDt << endl;

    // Delta T vs Momentum Plots (Separated)
    TCanvas *canvasDtVsPRec = new TCanvas(("canvasDtVsPRec_" + suffix).c_str(), ("Delta T vs Momentum (REC::PID) " + suffix).c_str(), 1200, 1200);
    canvasDtVsPRec->Divide(3, 3);
    pad = 1;
    for (int pid : particle_types) {
        canvasDtVsPRec->cd(pad++);
        hDtVsPRec[pid]->Draw("COLZ");
        hDtVsPRec[pid]->SetStats(0);
        TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        string label = "REC::PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hDtVsPRec[pid]->GetEntries())) + ")";
        leg->AddEntry(hDtVsPRec[pid], label.c_str(), "p");
        leg->Draw();
    }
    string outputDtVsPRec = dir + "/delta_t_vs_p_rec_" + suffix + ".pdf";
    canvasDtVsPRec->Print(outputDtVsPRec.c_str());
    cout << "Delta T vs Momentum plots for REC::PID saved to " << outputDtVsPRec << endl;

    TCanvas *canvasDtVsPCustom = new TCanvas(("canvasDtVsPCustom_" + suffix).c_str(), ("Delta T vs Momentum (Assigned PID) " + suffix).c_str(), 1200, 1200);
    canvasDtVsPCustom->Divide(3, 3);
    pad = 1;
    for (int pid : particle_types) {
        canvasDtVsPCustom->cd(pad++);
        hDtVsPCustom[pid]->Draw("COLZ");
        hDtVsPCustom[pid]->SetStats(0);
        TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        string label = "Assigned PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hDtVsPCustom[pid]->GetEntries())) + ")";
        leg->AddEntry(hDtVsPCustom[pid], label.c_str(), "p");
        leg->Draw();
    }
    string outputDtVsPCustom = dir + "/delta_t_vs_p_custom_" + suffix + ".pdf";
    canvasDtVsPCustom->Print(outputDtVsPCustom.c_str());
    cout << "Delta T vs Momentum plots for Assigned PID saved to " << outputDtVsPCustom << endl;

    // Chi2pid vs Momentum Plots (Separated)
    TCanvas *canvasChi2pidVsPRec = new TCanvas(("canvasChi2pidVsPRec_" + suffix).c_str(), ("Chi2pid vs Momentum (REC::PID) " + suffix).c_str(), 1200, 1200);
    canvasChi2pidVsPRec->Divide(3, 3);
    pad = 1;
    for (int pid : particle_types) {
        canvasChi2pidVsPRec->cd(pad++);
        hChi2pidVsPRec[pid]->Draw("COLZ");
        hChi2pidVsPRec[pid]->SetStats(0);
        TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        string label = "REC::PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hChi2pidVsPRec[pid]->GetEntries())) + ")";
        leg->AddEntry(hChi2pidVsPRec[pid], label.c_str(), "p");
        leg->Draw();
    }
    string outputChi2pidVsPRec = dir + "/chi2pid_vs_p_rec_" + suffix + ".pdf";
    canvasChi2pidVsPRec->Print(outputChi2pidVsPRec.c_str());
    cout << "Chi2pid vs Momentum plots for REC::PID saved to " << outputChi2pidVsPRec << endl;

    TCanvas *canvasChi2pidVsPCustom = new TCanvas(("canvasChi2pidVsPCustom_" + suffix).c_str(), ("Chi2pid vs Momentum (Assigned PID) " + suffix).c_str(), 1200, 1200);
    canvasChi2pidVsPCustom->Divide(3, 3);
    pad = 1;
    for (int pid : particle_types) {
        canvasChi2pidVsPCustom->cd(pad++);
        hChi2pidVsPCustom[pid]->Draw("COLZ");
        hChi2pidVsPCustom[pid]->SetStats(0);
        TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        string label = "Assigned PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hChi2pidVsPCustom[pid]->GetEntries())) + ")";
        leg->AddEntry(hChi2pidVsPCustom[pid], label.c_str(), "p");
        leg->Draw();
    }
    string outputChi2pidVsPCustom = dir + "/chi2pid_vs_p_custom_" + suffix + ".pdf";
    canvasChi2pidVsPCustom->Print(outputChi2pidVsPCustom.c_str());
    cout << "Chi2pid vs Momentum plots for Assigned PID saved to " << outputChi2pidVsPCustom << endl;

    // Beta vs Momentum Plots (Separated)
    TCanvas *canvasBetaVsPRec = new TCanvas(("canvasBetaVsPRec_" + suffix).c_str(), ("Beta vs Momentum (REC::PID) " + suffix).c_str(), 1200, 1200);
    canvasBetaVsPRec->Divide(3, 3);
    pad = 1;
    for (int pid : particle_types) {
        canvasBetaVsPRec->cd(pad++);
        hBetaVsPRec[pid]->Draw("COLZ");
        hBetaVsPRec[pid]->SetStats(0);
        TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        string label = "REC::PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hBetaVsPRec[pid]->GetEntries())) + ")";
        leg->AddEntry(hBetaVsPRec[pid], label.c_str(), "p");
        leg->Draw();
    }
    string outputBetaVsPRec = dir + "/beta_vs_p_rec_" + suffix + ".pdf";
    canvasBetaVsPRec->Print(outputBetaVsPRec.c_str());
    cout << "Beta vs Momentum plots for REC::PID saved to " << outputBetaVsPRec << endl;

    TCanvas *canvasBetaVsPCustom = new TCanvas(("canvasBetaVsPCustom_" + suffix).c_str(), ("Beta vs Momentum (Assigned PID) " + suffix).c_str(), 1200, 1200);
    canvasBetaVsPCustom->Divide(3, 3);
    pad = 1;
    for (int pid : particle_types) {
        canvasBetaVsPCustom->cd(pad++);
        hBetaVsPCustom[pid]->Draw("COLZ");
        hBetaVsPCustom[pid]->SetStats(0);
        TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        string label = "Assigned PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hBetaVsPCustom[pid]->GetEntries())) + ")";
        leg->AddEntry(hBetaVsPCustom[pid], label.c_str(), "p");
        leg->Draw();
    }
    string outputBetaVsPCustom = dir + "/beta_vs_p_custom_" + suffix + ".pdf";
    canvasBetaVsPCustom->Print(outputBetaVsPCustom.c_str());
    cout << "Beta vs Momentum plots for Assigned PID saved to " << outputBetaVsPCustom << endl;

    // Clean up canvases and histograms
    delete hChi2pidAll;
    delete hChi2pidCustomAll;
    delete hChi2pidDiff;
    delete canvasChi2pid;
    delete canvasMomentum;
    delete canvasChi2pidPerPid;
    delete canvasDt;
    delete canvasDtVsPRec;
    delete canvasDtVsPCustom;
    delete canvasChi2pidVsPRec;
    delete canvasChi2pidVsPCustom;
    delete canvasBetaVsPRec;
    delete canvasBetaVsPCustom;
}

int main() {
    TStopwatch timer;
    timer.Start();

    // Create directories for "before" and "after" trigger outputs
    fs::create_directory("before_trigger");
    fs::create_directory("after_trigger");

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
        cerr << "Error: No directories found" << endl;
        return 1;
    }

    loadCCDBParams();

    const int maxEvents = 100000;
    int event_count = 0;
    int total_events_processed = 0;
    int directory_count = 0;
    int hipo_file_count = 0;
    set<int> unique_events;
    ParticleStats statsBefore, statsAfter;
    vector<Mismatch> mismatchesBefore, mismatchesAfter;
    vector<float> chi2pidListBefore, chi2pidListAfter;
    vector<float> chi2pidCustomListBefore, chi2pidCustomListAfter;

    // Initialize histograms for "before" trigger
    map<int, TH1F*> hMomentumRecBefore, hMomentumCustomBefore;
    map<int, TH1F*> hChi2pidRecBefore, hChi2pidCustomBefore;
    map<int, TH1F*> hDtRecBefore, hDtCustomBefore;
    map<int, TH2F*> hDtVsPRecBefore, hDtVsPCustomBefore;
    map<int, TH2F*> hChi2pidVsPRecBefore, hChi2pidVsPCustomBefore;
    map<int, TH2F*> hBetaVsPRecBefore, hBetaVsPCustomBefore;

    // Initialize histograms for "after" trigger
    map<int, TH1F*> hMomentumRecAfter, hMomentumCustomAfter;
    map<int, TH1F*> hChi2pidRecAfter, hChi2pidCustomAfter;
    map<int, TH1F*> hDtRecAfter, hDtCustomAfter;
    map<int, TH2F*> hDtVsPRecAfter, hDtVsPCustomAfter;
    map<int, TH2F*> hChi2pidVsPRecAfter, hChi2pidVsPCustomAfter;
    map<int, TH2F*> hBetaVsPRecAfter, hBetaVsPCustomAfter;

    vector<int> particle_types = {11, -11, 211, -211, 321, -321, 2212, -2212, 45};

    // Initialize momentum histograms
    for (int pid : particle_types) {
        string name_rec_before = "hMomentumRecBefore_" + to_string(pid);
        string name_custom_before = "hMomentumCustomBefore_" + to_string(pid);
        string name_rec_after = "hMomentumRecAfter_" + to_string(pid);
        string name_custom_after = "hMomentumCustomAfter_" + to_string(pid);
        string title_rec_before = "Momentum (REC::PID Before Trigger) " + to_string(pid) + ";Momentum (GeV);Counts";
        string title_custom_before = "Momentum (Assigned PID Before Trigger) " + to_string(pid) + ";Momentum (GeV);Counts";
        string title_rec_after = "Momentum (REC::PID After Trigger) " + to_string(pid) + ";Momentum (GeV);Counts";
        string title_custom_after = "Momentum (Assigned PID After Trigger) " + to_string(pid) + ";Momentum (GeV);Counts";
        hMomentumRecBefore[pid] = new TH1F(name_rec_before.c_str(), title_rec_before.c_str(), 200, 0, 10);
        hMomentumCustomBefore[pid] = new TH1F(name_custom_before.c_str(), title_custom_before.c_str(), 200, 0, 10);
        hMomentumRecAfter[pid] = new TH1F(name_rec_after.c_str(), title_rec_after.c_str(), 200, 0, 10);
        hMomentumCustomAfter[pid] = new TH1F(name_custom_after.c_str(), title_custom_after.c_str(), 200, 0, 10);
        hMomentumRecBefore[pid]->SetLineColor(kBlue);
        hMomentumRecBefore[pid]->SetFillColor(kBlue);
        hMomentumRecBefore[pid]->SetFillStyle(3004);
        hMomentumCustomBefore[pid]->SetLineColor(kRed);
        hMomentumCustomBefore[pid]->SetFillColor(kRed);
        hMomentumCustomBefore[pid]->SetFillStyle(3005);
        hMomentumRecAfter[pid]->SetLineColor(kBlue);
        hMomentumRecAfter[pid]->SetFillColor(kBlue);
        hMomentumRecAfter[pid]->SetFillStyle(3004);
        hMomentumCustomAfter[pid]->SetLineColor(kRed);
        hMomentumCustomAfter[pid]->SetFillColor(kRed);
        hMomentumCustomAfter[pid]->SetFillStyle(3005);
    }

    // Initialize chi2pid histograms
    for (int pid : particle_types) {
        string name_rec_before = "hChi2pidRecBefore_" + to_string(pid);
        string name_custom_before = "hChi2pidCustomBefore_" + to_string(pid);
        string name_rec_after = "hChi2pidRecAfter_" + to_string(pid);
        string name_custom_after = "hChi2pidCustomAfter_" + to_string(pid);
        string title_rec_before = "Chi2pid (REC::PID Before Trigger) " + to_string(pid) + ";Chi2pid;Counts";
        string title_custom_before = "Chi2pid_custom (Assigned PID Before Trigger) " + to_string(pid) + ";Chi2pid;Counts";
        string title_rec_after = "Chi2pid (REC::PID After Trigger) " + to_string(pid) + ";Chi2pid;Counts";
        string title_custom_after = "Chi2pid_custom (Assigned PID After Trigger) " + to_string(pid) + ";Chi2pid;Counts";
        hChi2pidRecBefore[pid] = new TH1F(name_rec_before.c_str(), title_rec_before.c_str(), 200, -15, 15);
        hChi2pidCustomBefore[pid] = new TH1F(name_custom_before.c_str(), title_custom_before.c_str(), 200, -15, 15);
        hChi2pidRecAfter[pid] = new TH1F(name_rec_after.c_str(), title_rec_after.c_str(), 200, -15, 15);
        hChi2pidCustomAfter[pid] = new TH1F(name_custom_after.c_str(), title_custom_after.c_str(), 200, -15, 15);
        hChi2pidRecBefore[pid]->SetLineColor(kBlue);
        hChi2pidRecBefore[pid]->SetFillColor(kBlue);
        hChi2pidRecBefore[pid]->SetFillStyle(3004);
        hChi2pidCustomBefore[pid]->SetLineColor(kRed);
        hChi2pidCustomBefore[pid]->SetFillColor(kRed);
        hChi2pidCustomBefore[pid]->SetFillStyle(3005);
        hChi2pidRecAfter[pid]->SetLineColor(kBlue);
        hChi2pidRecAfter[pid]->SetFillColor(kBlue);
        hChi2pidRecAfter[pid]->SetFillStyle(3004);
        hChi2pidCustomAfter[pid]->SetLineColor(kRed);
        hChi2pidCustomAfter[pid]->SetFillColor(kRed);
        hChi2pidCustomAfter[pid]->SetFillStyle(3005);
    }

    // Initialize Delta T histograms
    for (int pid : particle_types) {
        string name_rec_before = "hDtRecBefore_" + to_string(pid);
        string name_custom_before = "hDtCustomBefore_" + to_string(pid);
        string name_rec_after = "hDtRecAfter_" + to_string(pid);
        string name_custom_after = "hDtCustomAfter_" + to_string(pid);
        string title_rec_before = "Delta T (REC::PID Before Trigger) " + to_string(pid) + ";Delta T (ns);Counts";
        string title_custom_before = "Delta T (Assigned PID Before Trigger) " + to_string(pid) + ";Delta T (ns);Counts";
        string title_rec_after = "Delta T (REC::PID After Trigger) " + to_string(pid) + ";Delta T (ns);Counts";
        string title_custom_after = "Delta T (Assigned PID After Trigger) " + to_string(pid) + ";Delta T (ns);Counts";
        hDtRecBefore[pid] = new TH1F(name_rec_before.c_str(), title_rec_before.c_str(), 200, -5, 5);
        hDtCustomBefore[pid] = new TH1F(name_custom_before.c_str(), title_custom_before.c_str(), 200, -5, 5);
        hDtRecAfter[pid] = new TH1F(name_rec_after.c_str(), title_rec_after.c_str(), 200, -5, 5);
        hDtCustomAfter[pid] = new TH1F(name_custom_after.c_str(), title_custom_after.c_str(), 200, -5, 5);
        hDtRecBefore[pid]->SetLineColor(kBlue);
        hDtRecBefore[pid]->SetFillColor(kBlue);
        hDtRecBefore[pid]->SetFillStyle(3004);
        hDtCustomBefore[pid]->SetLineColor(kRed);
        hDtCustomBefore[pid]->SetFillColor(kRed);
        hDtCustomBefore[pid]->SetFillStyle(3005);
        hDtRecAfter[pid]->SetLineColor(kBlue);
        hDtRecAfter[pid]->SetFillColor(kBlue);
        hDtRecAfter[pid]->SetFillStyle(3004);
        hDtCustomAfter[pid]->SetLineColor(kRed);
        hDtCustomAfter[pid]->SetFillColor(kRed);
        hDtCustomAfter[pid]->SetFillStyle(3005);
    }

    // Initialize 2D histograms for Delta T vs Momentum
    for (int pid : particle_types) {
        string name_rec_before = "hDtVsPRecBefore_" + to_string(pid);
        string name_custom_before = "hDtVsPCustomBefore_" + to_string(pid);
        string name_rec_after = "hDtVsPRecAfter_" + to_string(pid);
        string name_custom_after = "hDtVsPCustomAfter_" + to_string(pid);
        string title_rec_before = "Delta T vs Momentum (REC::PID Before Trigger) " + to_string(pid) + ";Momentum (GeV);Delta T (ns)";
        string title_custom_before = "Delta T vs Momentum (Assigned PID Before Trigger) " + to_string(pid) + ";Momentum (GeV);Delta T (ns)";
        string title_rec_after = "Delta T vs Momentum (REC::PID After Trigger) " + to_string(pid) + ";Momentum (GeV);Delta T (ns)";
        string title_custom_after = "Delta T vs Momentum (Assigned PID After Trigger) " + to_string(pid) + ";Momentum (GeV);Delta T (ns)";
        hDtVsPRecBefore[pid] = new TH2F(name_rec_before.c_str(), title_rec_before.c_str(), 200, 0, 10, 200, -5, 5);
        hDtVsPCustomBefore[pid] = new TH2F(name_custom_before.c_str(), title_custom_before.c_str(), 200, 0, 10, 200, -5, 5);
        hDtVsPRecAfter[pid] = new TH2F(name_rec_after.c_str(), title_rec_after.c_str(), 200, 0, 10, 200, -5, 5);
        hDtVsPCustomAfter[pid] = new TH2F(name_custom_after.c_str(), title_custom_after.c_str(), 200, 0, 10, 200, -5, 5);
        hDtVsPRecBefore[pid]->SetMarkerColor(kBlue);
        hDtVsPCustomBefore[pid]->SetMarkerColor(kRed);
        hDtVsPRecAfter[pid]->SetMarkerColor(kBlue);
        hDtVsPCustomAfter[pid]->SetMarkerColor(kRed);
    }

    // Initialize 2D histograms for Chi2pid vs Momentum
    for (int pid : particle_types) {
        string name_rec_before = "hChi2pidVsPRecBefore_" + to_string(pid);
        string name_custom_before = "hChi2pidVsPCustomBefore_" + to_string(pid);
        string name_rec_after = "hChi2pidVsPRecAfter_" + to_string(pid);
        string name_custom_after = "hChi2pidVsPCustomAfter_" + to_string(pid);
        string title_rec_before = "Chi2pid vs Momentum (REC::PID Before Trigger) " + to_string(pid) + ";Momentum (GeV);Chi2pid";
        string title_custom_before = "Chi2pid vs Momentum (Assigned PID Before Trigger) " + to_string(pid) + ";Momentum (GeV);Chi2pid";
        string title_rec_after = "Chi2pid vs Momentum (REC::PID After Trigger) " + to_string(pid) + ";Momentum (GeV);Chi2pid";
        string title_custom_after = "Chi2pid vs Momentum (Assigned PID After Trigger) " + to_string(pid) + ";Momentum (GeV);Chi2pid";
        hChi2pidVsPRecBefore[pid] = new TH2F(name_rec_before.c_str(), title_rec_before.c_str(), 200, 0, 10, 200, -15, 15);
        hChi2pidVsPCustomBefore[pid] = new TH2F(name_custom_before.c_str(), title_custom_before.c_str(), 200, 0, 10, 200, -15, 15);
        hChi2pidVsPRecAfter[pid] = new TH2F(name_rec_after.c_str(), title_rec_after.c_str(), 200, 0, 10, 200, -15, 15);
        hChi2pidVsPCustomAfter[pid] = new TH2F(name_custom_after.c_str(), title_custom_after.c_str(), 200, 0, 10, 200, -15, 15);
        hChi2pidVsPRecBefore[pid]->SetMarkerColor(kBlue);
        hChi2pidVsPCustomBefore[pid]->SetMarkerColor(kRed);
        hChi2pidVsPRecAfter[pid]->SetMarkerColor(kBlue);
        hChi2pidVsPCustomAfter[pid]->SetMarkerColor(kRed);
    }

    // Initialize beta vs. momentum histograms
    for (int pid : particle_types) {
        string name_rec_before = "hBetaVsPRecBefore_" + to_string(pid);
        string name_custom_before = "hBetaVsPCustomBefore_" + to_string(pid);
        string name_rec_after = "hBetaVsPRecAfter_" + to_string(pid);
        string name_custom_after = "hBetaVsPCustomAfter_" + to_string(pid);
        string title_rec_before = "Beta vs Momentum (REC::PID Before Trigger) " + to_string(pid) + ";Momentum (GeV);Beta";
        string title_custom_before = "Beta vs Momentum (Assigned PID Before Trigger) " + to_string(pid) + ";Momentum (GeV);Beta";
        string title_rec_after = "Beta vs Momentum (REC::PID After Trigger) " + to_string(pid) + ";Momentum (GeV);Beta";
        string title_custom_after = "Beta vs Momentum (Assigned PID After Trigger) " + to_string(pid) + ";Momentum (GeV);Beta";
        hBetaVsPRecBefore[pid] = new TH2F(name_rec_before.c_str(), title_rec_before.c_str(), 200, 0, 10, 200, 0, 1.2);
        hBetaVsPCustomBefore[pid] = new TH2F(name_custom_before.c_str(), title_custom_before.c_str(), 200, 0, 10, 200, 0, 1.2);
        hBetaVsPRecAfter[pid] = new TH2F(name_rec_after.c_str(), title_rec_after.c_str(), 200, 0, 10, 200, 0, 1.2);
        hBetaVsPCustomAfter[pid] = new TH2F(name_custom_after.c_str(), title_custom_after.c_str(), 200, 0, 10, 200, 0, 1.2);
        int color = (pid == 11 || pid == -11) ? kBlue : (pid == 211 || pid == -211) ? kRed : (pid == 321 || pid == -321) ? kGreen+2 : (pid == 2212 || pid == -2212) ? kMagenta : kCyan;
        hBetaVsPRecBefore[pid]->SetMarkerColor(color);
        hBetaVsPCustomBefore[pid]->SetMarkerColor(color);
        hBetaVsPRecAfter[pid]->SetMarkerColor(color);
        hBetaVsPCustomAfter[pid]->SetMarkerColor(color);
        hBetaVsPRecBefore[pid]->SetMarkerStyle(20);
        hBetaVsPCustomBefore[pid]->SetMarkerStyle(20);
        hBetaVsPRecAfter[pid]->SetMarkerStyle(20);
        hBetaVsPCustomAfter[pid]->SetMarkerStyle(20);
        hBetaVsPRecBefore[pid]->SetMarkerSize(0.5);
        hBetaVsPCustomBefore[pid]->SetMarkerSize(0.5);
        hBetaVsPRecAfter[pid]->SetMarkerSize(0.5);
        hBetaVsPCustomAfter[pid]->SetMarkerSize(0.5);
    }

    int events_with_trigger_rec = 0;
int events_with_trigger_custom = 0;
int total_trigger_electrons_rec = 0;  // Total number of REC::Particle trigger electrons
int total_trigger_electrons_custom = 0;  // Total number of custom-assigned trigger electrons
int total_electrons_rec = 0;  // Total number of REC::Particle electrons (pid == 11)
int total_electrons_custom = 0;  // Total number of custom-assigned electrons (assigned_pid == 11)
for (const auto& dir : directories) {
    directory_count++;
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
        hipo_file_count++;
        cout << "Processing file: " << file << endl;
        hipo::reader reader;
        reader.open(file.c_str());
        hipo::dictionary factory;
        reader.readDictionary(factory);

        while (reader.next() && event_count < maxEvents) {
            event_count++;
            unique_events.insert(event_count);
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

            set<tuple<int, int, int, float, float, float>> unique_particles_per_event;
            for (int i = 0; i < PART.getRows(); i++) {
                int pid = PART.getInt("pid", i);
                int charge = PART.getByte("charge", i);
                int status = PART.getShort("status", i);
                float px = PART.getFloat("px", i);
                float py = PART.getFloat("py", i);
                float pz = PART.getFloat("pz", i);
                unique_particles_per_event.insert(make_tuple(pid, charge, status, px, py, pz));
            }
            if (unique_particles_per_event.size() != PART.getRows()) {
                cout << "Event " << event_count << ": Found duplicates! Unique particles = "
                     << unique_particles_per_event.size() << ", Total rows = " << PART.getRows() << endl;
            }

            IndexMap cherMap = loadMapByIndex(CHER, "pindex");
            IndexMap caloMap = loadMapByIndex(CALO, "pindex");
            IndexMap scinMap = loadMapByIndex(SCIN, "pindex");

            float startTime = EVENT.getFloat("startTime", 0);

            map<tuple<int, int, int, float, float, float>, int> particle_indices;
            vector<int> unique_indices;
            for (int i = 0; i < PART.getRows(); i++) {
                int pid = PART.getInt("pid", i);
                int charge = PART.getByte("charge", i);
                int status = PART.getShort("status", i);
                float px = PART.getFloat("px", i);
                float py = PART.getFloat("py", i);
                float pz = PART.getFloat("pz", i);
                auto key = make_tuple(pid, charge, status, px, py, pz);
                if (particle_indices.find(key) == particle_indices.end()) {
                    particle_indices[key] = i;
                    unique_indices.push_back(i);
                }
            }

           vector<ParticleData> particles(unique_indices.size());
int valid_particles = 0;
bool has_trigger_electron_rec = false;
bool has_trigger_electron_custom = false;
int trigger_electrons_rec_in_event = 0;
int trigger_electrons_custom_in_event = 0;

// First pass: Collect particles and count REC::Particle electrons and trigger electrons
for (int idx : unique_indices) {
    ParticleData pd = getParticleData(idx, PART, CHER, CALO, SCIN, cherMap, caloMap, scinMap);
    statsBefore.total_particles_before_filtering++;
    if (pd.pid == 0) {
        statsBefore.total_unidentified_before_filtering++;
    }
    if (pd.status == 0 || pd.charge == 0) {
        statsBefore.total_filtered_out++;
        if (pd.pid == 0) {
            statsBefore.total_filtered_out_pid_zero++;
        }
        continue;
    }
    particles[valid_particles++] = pd;

    // Count REC::Particle electrons and trigger electrons
    if (pd.pid == 11) {
        total_electrons_rec++;
        // Apply stricter trigger electron condition
        if (pd.status < 0 /* && pd.chi2pid > -5 */ && pd.chi2pid < 5 && 
            pd.vz_tele >= -20 && pd.vz_tele <= 5 && abs(pd.status)/1000 == 2) {
            has_trigger_electron_rec = true;
            trigger_electrons_rec_in_event++;
            total_trigger_electrons_rec++;
        }
    }
}
particles.resize(valid_particles);

// Assign PIDs to determine custom electrons and trigger electrons
assignPids(particles, startTime, event_count);

// Second pass: Count custom-assigned electrons and trigger electrons
for (const auto& pd : particles) {
    if (pd.assigned_pid == 11) {
        total_electrons_custom++;
        // Apply stricter trigger electron condition
        if (pd.status < 0 /* && pd.chi2pid_custom > -5 */ && pd.chi2pid_custom < 5 && 
            pd.vz_tele >= -20 && pd.vz_tele <= 5 && abs(pd.status)/1000 == 2) {
            has_trigger_electron_custom = true;
            trigger_electrons_custom_in_event++;
            total_trigger_electrons_custom++;
        }
    }
}

// Debug logging
cout << "Event " << event_count << ": has_trigger_electron_rec=" << has_trigger_electron_rec
     << ", has_trigger_electron_custom=" << has_trigger_electron_custom
     << ", trigger_electrons_rec=" << trigger_electrons_rec_in_event
     << ", trigger_electrons_custom=" << trigger_electrons_custom_in_event << endl;

// Update event counts
if (has_trigger_electron_rec) events_with_trigger_rec++;
if (has_trigger_electron_custom) events_with_trigger_custom++;

// Process the event for both cases
validateAndCompare(particles, event_count, statsBefore, statsAfter, mismatchesBefore, mismatchesAfter,
                   chi2pidListBefore, chi2pidListAfter, chi2pidCustomListBefore, chi2pidCustomListAfter,
                   hMomentumRecBefore, hMomentumCustomBefore, hMomentumRecAfter, hMomentumCustomAfter,
                   hChi2pidRecBefore, hChi2pidCustomBefore, hChi2pidRecAfter, hChi2pidCustomAfter,
                   hDtRecBefore, hDtCustomBefore, hDtRecAfter, hDtCustomAfter,
                   hDtVsPRecBefore, hDtVsPCustomBefore, hDtVsPRecAfter, hDtVsPCustomAfter,
                   hChi2pidVsPRecBefore, hChi2pidVsPCustomBefore, hChi2pidVsPRecAfter, hChi2pidVsPCustomAfter,
                   hBetaVsPRecBefore, hBetaVsPCustomBefore, hBetaVsPRecAfter, hBetaVsPCustomAfter,
                   has_trigger_electron_rec, has_trigger_electron_custom);
        }
        if (event_count >= maxEvents) break;
    }
    if (event_count >= maxEvents) break;
}

    timer.Stop();

    // Print Summary
    cout << "****************************************************************************************" << endl;
    cout << "Total number of events processed: " << event_count << endl;
    cout << "Total unique events: " << unique_events.size() << endl;
    cout << "Total electrons (REC::Particle): " << total_electrons_rec << endl;
    cout << "Total electrons (Assigned PID): " << total_electrons_custom << endl;
    cout << "Events with trigger electrons (REC::Particle): " << events_with_trigger_rec << endl;
    cout << "Events with trigger electrons (Assigned PID): " << events_with_trigger_custom << endl;

    // Save statistics to CSV
    saveStatsToCSV(statsBefore, "before_trigger/stats.csv", "Before Trigger: ");
    saveStatsToCSV(statsAfter, "after_trigger/stats.csv", "After Trigger: ");

    // Save plots
    savePlots("before_trigger", "before_trigger", chi2pidListBefore, chi2pidCustomListBefore,
              hMomentumRecBefore, hMomentumCustomBefore, hChi2pidRecBefore, hChi2pidCustomBefore,
              hDtRecBefore, hDtCustomBefore, hDtVsPRecBefore, hDtVsPCustomBefore,
              hChi2pidVsPRecBefore, hChi2pidVsPCustomBefore, hBetaVsPRecBefore, hBetaVsPCustomBefore);
    savePlots("after_trigger", "after_trigger", chi2pidListAfter, chi2pidCustomListAfter,
              hMomentumRecAfter, hMomentumCustomAfter, hChi2pidRecAfter, hChi2pidCustomAfter,
              hDtRecAfter, hDtCustomAfter, hDtVsPRecAfter, hDtVsPCustomAfter,
              hChi2pidVsPRecAfter, hChi2pidVsPCustomAfter, hBetaVsPRecAfter, hBetaVsPCustomAfter);

    // Print Mismatch Details
    cout << "\nMismatches (Before Trigger):\n";
    if (!mismatchesBefore.empty()) {
        for (const auto& mismatch : mismatchesBefore) {
            cout << "Event " << mismatch.eventNum << ", Particle Index " << mismatch.particleIdx
                 << ": REC::PID=" << mismatch.recPid << ", Assigned PID=" << mismatch.assignedPid
                 << ", Charge=" << mismatch.charge << ", Momentum=" << mismatch.momentum
                 << ", NPHE_HTCC=" << mismatch.npheHtcc << ", NPHE_LTCC=" << mismatch.npheLtcc
                 << ", Status=" << mismatch.status
                 << ", Chi2pid (REC)=" << mismatch.chi2pid
                 << ", Chi2pid_custom=" << mismatch.chi2pidCustom
                 << ", Diff=" << (mismatch.chi2pidCustom - mismatch.chi2pid)
                 << ", Has_FTOF=" << (mismatch.hasFtof ? "Yes" : "No")
                 << ", Has_CTOF=" << (mismatch.hasCtof ? "Yes" : "No") << endl;
        }
        cout << "Total mismatches (Before Trigger): " << mismatchesBefore.size() << endl;
    } else {
        cout << "No mismatches found (Before Trigger).\n";
    }

    cout << "\nMismatches (After Trigger):\n";
    if (!mismatchesAfter.empty()) {
        for (const auto& mismatch : mismatchesAfter) {
            cout << "Event " << mismatch.eventNum << ", Particle Index " << mismatch.particleIdx
                 << ": REC::PID=" << mismatch.recPid << ", Assigned PID=" << mismatch.assignedPid
                 << ", Charge=" << mismatch.charge << ", Momentum=" << mismatch.momentum
                 << ", NPHE_HTCC=" << mismatch.npheHtcc << ", NPHE_LTCC=" << mismatch.npheLtcc
                 << ", Status=" << mismatch.status
                 << ", Chi2pid (REC)=" << mismatch.chi2pid
                 << ", Chi2pid_custom=" << mismatch.chi2pidCustom
                 << ", Diff=" << (mismatch.chi2pidCustom - mismatch.chi2pid)
                 << ", Has_FTOF=" << (mismatch.hasFtof ? "Yes" : "No")
                 << ", Has_CTOF=" << (mismatch.hasCtof ? "Yes" : "No") << endl;
        }
        cout << "Total mismatches (After Trigger): " << mismatchesAfter.size() << endl;
    } else {
        cout << "No mismatches found (After Trigger).\n";
    }

    cout << "\nTotal events processed: " << total_events_processed << endl;

    // Print Summary
    cout << "Real time: " << timer.RealTime() << " s, CPU time: " << timer.CpuTime() << " s" << endl;

    // Clean up histograms
    for (auto& h : hMomentumRecBefore) delete h.second;
    for (auto& h : hMomentumCustomBefore) delete h.second;
    for (auto& h : hChi2pidRecBefore) delete h.second;
    for (auto& h : hChi2pidCustomBefore) delete h.second;
    for (auto& h : hDtRecBefore) delete h.second;
    for (auto& h : hDtCustomBefore) delete h.second;
    for (auto& h : hDtVsPRecBefore) delete h.second;
    for (auto& h : hDtVsPCustomBefore) delete h.second;
    for (auto& h : hChi2pidVsPRecBefore) delete h.second;
    for (auto& h : hChi2pidVsPCustomBefore) delete h.second;
    for (auto& h : hBetaVsPRecBefore) delete h.second;
    for (auto& h : hBetaVsPCustomBefore) delete h.second;

    for (auto& h : hMomentumRecAfter) delete h.second;
    for (auto& h : hMomentumCustomAfter) delete h.second;
    for (auto& h : hChi2pidRecAfter) delete h.second;
    for (auto& h : hChi2pidCustomAfter) delete h.second;
    for (auto& h : hDtRecAfter) delete h.second;
    for (auto& h : hDtCustomAfter) delete h.second;
    for (auto& h : hDtVsPRecAfter) delete h.second;
    for (auto& h : hDtVsPCustomAfter) delete h.second;
    for (auto& h : hChi2pidVsPRecAfter) delete h.second;
    for (auto& h : hChi2pidVsPCustomAfter) delete h.second;
    for (auto& h : hBetaVsPRecAfter) delete h.second;
    for (auto& h : hBetaVsPCustomAfter) delete h.second;

    return 0;
}