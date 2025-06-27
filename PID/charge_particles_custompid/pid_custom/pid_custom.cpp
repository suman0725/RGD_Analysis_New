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
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include <TStopwatch.h>
#include <TLegend.h>
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
    float vt;                // Vertex time (ns) from REC::Particle
    float beta;              // Beta from REC::Particle (measured by TOF)
    std::vector<DetectorHit> hits; // Detector hits (FTOF, CTOF, ECAL)
    bool is_trigger;         // True if particle is a trigger particle (status < 0)
    int assigned_pid;        // Custom-assigned PID
    float start_time;        // Start time for timing calculations (vt for charged particles)
    int status;              // Status from REC::Particle (detector collection passed)
    float chi2pid_custom;    // Computed PID quality for assigned PID
    float vz_tele;           // Vertex z-position from EVENT bank

    ParticleData() : charge(0), p(0.0f), sector(-1), nphe_htcc(0.0f), nphe_ltcc(0.0f),
                     energy_pcal(0.0f), energy_total(0.0f), vt(0.0f), beta(-99.0f),
                     hits(), is_trigger(false), assigned_pid(0), start_time(0.0f), status(0),
                     chi2pid_custom(99999.0f), vz_tele(0.0f) {}
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

// Statistics for custom PID
struct ParticleStats {
    int total_particles_before_filtering = 0;
    int total_filtered_out = 0;
    int total_assigned_pid_zero = 0;
    int valid_chi2pid_custom = 0;
    int invalid_chi2pid_custom = 0;
    std::map<int, int> assigned_pid_counts;               // Overall assigned PID counts
    std::map<Region, std::map<int, int>> assigned_pid_counts_by_region; // Assigned PID counts by region
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
void processEvent(std::vector<ParticleData>& particles, int eventNum, ParticleStats& stats,
                  std::vector<float>& chi2pidCustomList,
                  std::map<int, TH1F*>& hMomentumCustom,
                  std::map<int, TH1F*>& hChi2pidCustom,
                  std::map<int, TH1F*>& hDtCustom,
                  std::map<int, TH2F*>& hDtVsPCustom,
                  std::map<int, TH2F*>& hChi2pidVsPCustom,
                  std::map<int, TH2F*>& hBetaVsPCustom,
                  bool hasTriggerElectronCustom);
void saveStatsToCSV(const ParticleStats& stats, const string& filename);
void savePlots(const string& dir,
               std::vector<float>& chi2pidCustomList,
               std::map<int, TH1F*>& hMomentumCustom,
               std::map<int, TH1F*>& hChi2pidCustom,
               std::map<int, TH1F*>& hDtCustom,
               std::map<int, TH2F*>& hDtVsPCustom,
               std::map<int, TH2F*>& hChi2pidVsPCustom,
               std::map<int, TH2F*>& hBetaVsPCustom);

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
    float px = PART.getFloat("px", partIdx);
    float py = PART.getFloat("py", partIdx);
    float pz = PART.getFloat("pz", partIdx);
    pd.p = std::sqrt(px * px + py * py + pz * pz);
    pd.vt = PART.getFloat("vt", partIdx);
    pd.beta = PART.getFloat("beta", partIdx);
    pd.status = PART.getShort("status", partIdx);
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
                        if (abs(pid) == 11) continue; // Electrons are identified via ECAL, not timing
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

void processEvent(std::vector<ParticleData>& particles, int eventNum, ParticleStats& stats,
                 std::vector<float>& chi2pidCustomList,
                 std::map<int, TH1F*>& hMomentumCustom,
                 std::map<int, TH1F*>& hChi2pidCustom,
                 std::map<int, TH1F*>& hDtCustom,
                 std::map<int, TH2F*>& hDtVsPCustom,
                 std::map<int, TH2F*>& hChi2pidVsPCustom,
                 std::map<int, TH2F*>& hBetaVsPCustom,
                 bool hasTriggerElectronCustom) {
    if (!hasTriggerElectronCustom) return; // Only process events with a trigger electron

    for (size_t j = 0; j < particles.size(); ++j) {
        const auto& p = particles[j];
        if (p.charge == 0) continue; // Skip neutral particles

        Region region = determineRegion(p.status);
        stats.assigned_pid_counts[p.assigned_pid]++;
        stats.assigned_pid_counts_by_region[region][p.assigned_pid]++;

        if (p.assigned_pid == 0) stats.total_assigned_pid_zero++;

        if (hMomentumCustom.count(p.assigned_pid)) {
            hMomentumCustom[p.assigned_pid]->Fill(p.p);
        }

        cout << "Event " << eventNum << ", Particle " << j 
             << ": Assigned PID=" << p.assigned_pid 
             << ", Momentum=" << p.p << ", Measured Beta=" << p.beta 
             << ", Status=" << p.status << endl;

        if (p.p > 0 && p.beta >= 0 && p.beta <= 1.2 && p.assigned_pid != 0) {
            if (hBetaVsPCustom.count(p.assigned_pid)) {
                hBetaVsPCustom[p.assigned_pid]->Fill(p.p, p.beta);
            }
        }

        if (p.chi2pid_custom != 99999.0f) {
            stats.valid_chi2pid_custom++;
            chi2pidCustomList.push_back(p.chi2pid_custom);
            if (hChi2pidCustom.count(p.assigned_pid)) {
                hChi2pidCustom[p.assigned_pid]->Fill(p.chi2pid_custom);
            }
            if (hChi2pidVsPCustom.count(p.assigned_pid)) {
                hChi2pidVsPCustom[p.assigned_pid]->Fill(p.p, p.chi2pid_custom);
            }
        } else {
            stats.invalid_chi2pid_custom++;
        }

        float dt_custom = computeDeltaT(p, p.assigned_pid);
        if (dt_custom != 99999.0f && hDtCustom.count(p.assigned_pid)) {
            hDtCustom[p.assigned_pid]->Fill(dt_custom);
            if (hDtVsPCustom.count(p.assigned_pid)) {
                hDtVsPCustom[p.assigned_pid]->Fill(p.p, dt_custom);
            }
        }
    }
}

void saveStatsToCSV(const ParticleStats& stats, const string& filename) {
    ofstream csvFile(filename);
    if (!csvFile.is_open()) {
        cerr << "Error: Could not open " << filename << " for writing!" << endl;
        return;
    }

    csvFile << "Statistic,Value\n";
    csvFile << "Total particles before filtering," << stats.total_particles_before_filtering << "\n";
    csvFile << "Total particles filtered out (status == 0 || charge == 0)," << stats.total_filtered_out << "\n";
    csvFile << "Total particles with assigned PID == 0," << stats.total_assigned_pid_zero << "\n";
    csvFile << "Valid chi2pid_custom," << stats.valid_chi2pid_custom << "\n";
    csvFile << "Invalid chi2pid_custom," << stats.invalid_chi2pid_custom << "\n";

    csvFile << "\nAssigned PID Counts (Overall)\n";
    csvFile << "pid,Count\n";
    vector<int> particle_types = {11, -11, 211, -211, 321, -321, 2212, -2212, 45};
    for (int pid : particle_types) {
        int count = stats.assigned_pid_counts.count(pid) ? stats.assigned_pid_counts.at(pid) : 0;
        csvFile << pid << "," << count << "\n";
    }
    int count_pid0 = stats.assigned_pid_counts.count(0) ? stats.assigned_pid_counts.at(0) : 0;
    csvFile << "0," << count_pid0 << "\n";

    const std::map<Region, std::string> region_names = {
        {Region::Forward, "Forward"}, {Region::Central, "Central"}, {Region::Band, "Band"}, {Region::Unknown, "Unknown"}
    };
    for (const auto& region_pair : region_names) {
        csvFile << "\nAssigned PID Counts (" << region_pair.second << ")\n";
        csvFile << "pid,Count\n";
        for (int pid : particle_types) {
            int region_count = stats.assigned_pid_counts_by_region.count(region_pair.first) && 
                              stats.assigned_pid_counts_by_region.at(region_pair.first).count(pid) ?
                              stats.assigned_pid_counts_by_region.at(region_pair.first).at(pid) : 0;
            csvFile << pid << "," << region_count << "\n";
        }
        int region_count_pid0 = stats.assigned_pid_counts_by_region.count(region_pair.first) && 
                               stats.assigned_pid_counts_by_region.at(region_pair.first).count(0) ?
                               stats.assigned_pid_counts_by_region.at(region_pair.first).at(0) : 0;
        csvFile << "0," << region_count_pid0 << "\n";
    }

    csvFile.close();
    cout << "Saved statistics to " << filename << endl;
}

void savePlots(const string& dir,
               std::vector<float>& chi2pidCustomList,
               std::map<int, TH1F*>& hMomentumCustom,
               std::map<int, TH1F*>& hChi2pidCustom,
               std::map<int, TH1F*>& hDtCustom,
               std::map<int, TH2F*>& hDtVsPCustom,
               std::map<int, TH2F*>& hChi2pidVsPCustom,
               std::map<int, TH2F*>& hBetaVsPCustom) {
    vector<int> particle_types = {11, -11, 211, -211, 321, -321, 2212, -2212, 45};

    // Chi2pid_custom Plot (Combined for all charged particles)
    TCanvas *canvasChi2pid = new TCanvas("canvasChi2pid", "Chi2pid_custom for Charged Particles", 800, 600);
    TH1F *hChi2pidCustomAll = new TH1F("hChi2pidCustomAll", "Chi2pid_custom (Assigned);Chi2pid;Counts", 200, -15, 15);
    hChi2pidCustomAll->SetLineColor(kRed);
    hChi2pidCustomAll->SetFillColor(kRed);
    hChi2pidCustomAll->SetFillStyle(3005);

    for (float chi2 : chi2pidCustomList) {
        hChi2pidCustomAll->Fill(chi2);
    }

    hChi2pidCustomAll->Draw("HIST");
    string outputChi2pid = dir + "/chi2pid_custom.pdf";
    canvasChi2pid->Print(outputChi2pid.c_str());
    cout << "Chi2pid_custom histogram saved to " << outputChi2pid << endl;

    // Momentum Plots (Combined)
    TCanvas *canvasMomentum = new TCanvas("canvasMomentum", "Momentum for Charged Particles", 1200, 1200);
    canvasMomentum->Divide(3, 3);
    int pad = 1;
    for (int pid : particle_types) {
        canvasMomentum->cd(pad++);
        hMomentumCustom[pid]->Draw("HIST");
        TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        string label = "Assigned PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hMomentumCustom[pid]->GetEntries())) + ")";
        leg->AddEntry(hMomentumCustom[pid], label.c_str(), "f");
        leg->Draw();
    }

    string outputMomentum = dir + "/momentum_plots.pdf";
    canvasMomentum->Print(outputMomentum.c_str());
    cout << "Momentum plots saved to " << outputMomentum << endl;

    // Chi2pid_custom Plots (Combined)
    TCanvas *canvasChi2pidPerPid = new TCanvas("canvasChi2pidPerPid", "Chi2pid_custom Per Particle Type", 1200, 1200);
    canvasChi2pidPerPid->Divide(3, 3);
    pad = 1;
    for (int pid : particle_types) {
        canvasChi2pidPerPid->cd(pad++);
        hChi2pidCustom[pid]->Draw("HIST");
        TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        string label = "Assigned PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hChi2pidCustom[pid]->GetEntries())) + ")";
        leg->AddEntry(hChi2pidCustom[pid], label.c_str(), "f");
        leg->Draw();
    }

    string outputChi2pidPerPid = dir + "/chi2pid_custom_plots.pdf";
    canvasChi2pidPerPid->Print(outputChi2pidPerPid.c_str());
    cout << "Chi2pid_custom per particle plots saved to " << outputChi2pidPerPid << endl;

    // Delta T Plots (Combined)
    TCanvas *canvasDt = new TCanvas("canvasDt", "Delta T for Charged Particles", 1200, 1200);
    canvasDt->Divide(3, 3);
    pad = 1;
    for (int pid : particle_types) {
        canvasDt->cd(pad++);
        hDtCustom[pid]->Draw("HIST");
        TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        string label = "Assigned PID " + to_string(pid) + " (Entries: " + to_string(static_cast<int>(hDtCustom[pid]->GetEntries())) + ")";
        leg->AddEntry(hDtCustom[pid], label.c_str(), "f");
        leg->Draw();
    }

    string outputDt = dir + "/delta_t_plots.pdf";
    canvasDt->Print(outputDt.c_str());
    cout << "Delta T plots saved to " << outputDt << endl;

    // Delta T vs Momentum Plots
    TCanvas *canvasDtVsPCustom = new TCanvas("canvasDtVsPCustom", "Delta T vs Momentum (Assigned PID)", 1200, 1200);
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
    string outputDtVsPCustom = dir + "/delta_t_vs_p_custom.pdf";
    canvasDtVsPCustom->Print(outputDtVsPCustom.c_str());
    cout << "Delta T vs Momentum plots saved to " << outputDtVsPCustom << endl;

    // Chi2pid_custom vs Momentum Plots
    TCanvas *canvasChi2pidVsPCustom = new TCanvas("canvasChi2pidVsPCustom", "Chi2pid_custom vs Momentum (Assigned PID)", 1200, 1200);
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
    string outputChi2pidVsPCustom = dir + "/chi2pid_vs_p_custom.pdf";
    canvasChi2pidVsPCustom->Print(outputChi2pidVsPCustom.c_str());
    cout << "Chi2pid_custom vs Momentum plots saved to " << outputChi2pidVsPCustom << endl;

    // Beta vs Momentum Plots
    TCanvas *canvasBetaVsPCustom = new TCanvas("canvasBetaVsPCustom", "Beta vs Momentum (Assigned PID)", 1200, 1200);
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
    string outputBetaVsPCustom = dir + "/beta_vs_p_custom.pdf";
    canvasBetaVsPCustom->Print(outputBetaVsPCustom.c_str());
    cout << "Beta vs Momentum plots saved to " << outputBetaVsPCustom << endl;

    // Clean up
    delete hChi2pidCustomAll;
    delete canvasChi2pid;
    delete canvasMomentum;
    delete canvasChi2pidPerPid;
    delete canvasDt;
    delete canvasDtVsPCustom;
    delete canvasChi2pidVsPCustom;
    delete canvasBetaVsPCustom;
}

int main() {
    TStopwatch timer;
    timer.Start();

    // Create output directory
    fs::create_directory("output");

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
        cerr << "Error: No directories found" << endl;
        return 1;
    }

    loadCCDBParams();

    const int maxEvents = 100000;
    int event_count = 0;
    int total_events_processed = 0;
    set<int> unique_events;
    ParticleStats stats;
    vector<float> chi2pidCustomList;

    // Initialize histograms for custom PID
    map<int, TH1F*> hMomentumCustom;
    map<int, TH1F*> hChi2pidCustom;
    map<int, TH1F*> hDtCustom;
    map<int, TH2F*> hDtVsPCustom;
    map<int, TH2F*> hChi2pidVsPCustom;
    map<int, TH2F*> hBetaVsPCustom;

    vector<int> particle_types = {11, -11, 211, -211, 321, -321, 2212, -2212, 45};

    // Initialize momentum histograms
    for (int pid : particle_types) {
        string name = "hMomentumCustom_" + to_string(pid);
        string title = "Momentum (Assigned PID) " + to_string(pid) + ";Momentum (GeV);Counts";
        hMomentumCustom[pid] = new TH1F(name.c_str(), title.c_str(), 200, 0, 10);
        hMomentumCustom[pid]->SetLineColor(kRed);
        hMomentumCustom[pid]->SetFillColor(kRed);
        hMomentumCustom[pid]->SetFillStyle(3005);
    }

    // Initialize chi2pid_custom histograms
    for (int pid : particle_types) {
        string name = "hChi2pidCustom_" + to_string(pid);
        string title = "Chi2pid_custom (Assigned PID) " + to_string(pid) + ";Chi2pid;Counts";
        hChi2pidCustom[pid] = new TH1F(name.c_str(), title.c_str(), 200, -15, 15);
        hChi2pidCustom[pid]->SetLineColor(kRed);
        hChi2pidCustom[pid]->SetFillColor(kRed);
        hChi2pidCustom[pid]->SetFillStyle(3005);
    }

    // Initialize Delta T histograms
    for (int pid : particle_types) {
        string name = "hDtCustom_" + to_string(pid);
        string title = "Delta T (Assigned PID) " + to_string(pid) + ";Delta T (ns);Counts";
        hDtCustom[pid] = new TH1F(name.c_str(), title.c_str(), 200, -5, 5);
        hDtCustom[pid]->SetLineColor(kRed);
        hDtCustom[pid]->SetFillColor(kRed);
        hDtCustom[pid]->SetFillStyle(3005);
    }

    // Initialize 2D histograms for Delta T vs Momentum
    for (int pid : particle_types) {
        string name = "hDtVsPCustom_" + to_string(pid);
        string title = "Delta T vs Momentum (Assigned PID) " + to_string(pid) + ";Momentum (GeV);Delta T (ns)";
        hDtVsPCustom[pid] = new TH2F(name.c_str(), title.c_str(), 200, 0, 10, 200, -5, 5);
        hDtVsPCustom[pid]->SetMarkerColor(kRed);
    }

    // Initialize 2D histograms for Chi2pid_custom vs Momentum
    for (int pid : particle_types) {
        string name = "hChi2pidVsPCustom_" + to_string(pid);
        string title = "Chi2pid_custom vs Momentum (Assigned PID) " + to_string(pid) + ";Momentum (GeV);Chi2pid";
        hChi2pidVsPCustom[pid] = new TH2F(name.c_str(), title.c_str(), 200, 0, 10, 200, -15, 15);
        hChi2pidVsPCustom[pid]->SetMarkerColor(kRed);
    }

    // Initialize beta vs. momentum histograms
    for (int pid : particle_types) {
        string name = "hBetaVsPCustom_" + to_string(pid);
        string title = "Beta vs Momentum (Assigned PID) " + to_string(pid) + ";Momentum (GeV);Beta";
        hBetaVsPCustom[pid] = new TH2F(name.c_str(), title.c_str(), 200, 0, 10, 200, 0, 1.2);
        int color = (pid == 11 || pid == -11) ? kBlue : (pid == 211 || pid == -211) ? kRed : (pid == 321 || pid == -321) ? kGreen+2 : (pid == 2212 || pid == -2212) ? kMagenta : kCyan;
        hBetaVsPCustom[pid]->SetMarkerColor(color);
        hBetaVsPCustom[pid]->SetMarkerStyle(20);
        hBetaVsPCustom[pid]->SetMarkerSize(0.5);
    }

    int events_with_trigger_custom = 0;
    int total_trigger_electrons_custom = 0;
    int total_electrons_custom = 0;

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
            hipo::reader reader;
            reader.open(file.c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);

            if (!factory.hasSchema("REC::Particle") || !factory.hasSchema("REC::Event") ||
                !factory.hasSchema("REC::Scintillator") || !factory.hasSchema("REC::Cherenkov") ||
                !factory.hasSchema("REC::Calorimeter")) {
                cerr << "Error: Required schemas not found in file: " << file << endl;
                continue;
            }

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

                // Check for duplicate particles
                set<tuple<int, int, float, float, float>> unique_particles_per_event;
                for (int i = 0; i < PART.getRows(); i++) {
                    int charge = PART.getByte("charge", i);
                    int status = PART.getShort("status", i);
                    float px = PART.getFloat("px", i);
                    float py = PART.getFloat("py", i);
                    float pz = PART.getFloat("pz", i);
                    unique_particles_per_event.insert(make_tuple(charge, status, px, py, pz));
                }
                if (unique_particles_per_event.size() != PART.getRows()) {
                    cout << "Event " << event_count << ": Found duplicates! Unique particles = "
                         << unique_particles_per_event.size() << ", Total rows = " << PART.getRows() << endl;
                }

                IndexMap cherMap = loadMapByIndex(CHER, "pindex");
                IndexMap caloMap = loadMapByIndex(CALO, "pindex");
                IndexMap scinMap = loadMapByIndex(SCIN, "pindex");

                float startTime = EVENT.getFloat("startTime", 0);

                // Remove duplicates
                map<tuple<int, int, float, float, float>, int> particle_indices;
                vector<int> unique_indices;
                for (int i = 0; i < PART.getRows(); i++) {
                    int charge = PART.getByte("charge", i);
                    int status = PART.getShort("status", i);
                    float px = PART.getFloat("px", i);
                    float py = PART.getFloat("py", i);
                    float pz = PART.getFloat("pz", i);
                    auto key = make_tuple(charge, status, px, py, pz);
                    if (particle_indices.find(key) == particle_indices.end()) {
                        particle_indices[key] = i;
                        unique_indices.push_back(i);
                    }
                }

                vector<ParticleData> particles(unique_indices.size());
                int valid_particles = 0;
                bool has_trigger_electron_custom = false;
                int trigger_electrons_custom_in_event = 0;

                // Collect particles and filter
                for (int idx : unique_indices) {
                    ParticleData pd = getParticleData(idx, PART, CHER, CALO, SCIN, cherMap, caloMap, scinMap);
                    stats.total_particles_before_filtering++;
                    if (pd.status == 0 || pd.charge == 0) {
                        stats.total_filtered_out++;
                        continue;
                    }
                    particles[valid_particles++] = pd;
                }
                particles.resize(valid_particles);

                // Assign custom PIDs
                assignPids(particles, startTime, event_count);

                // Check for trigger electrons (custom PID)
                for (const auto& pd : particles) {
                    if (pd.assigned_pid == 11) {
                        total_electrons_custom++;
                        if (pd.status < 0 && pd.chi2pid_custom < 5 && 
                            pd.vz_tele >= -20 && pd.vz_tele <= 5 && abs(pd.status)/1000 == 2) {
                            has_trigger_electron_custom = true;
                            trigger_electrons_custom_in_event++;
                            total_trigger_electrons_custom++;
                        }
                    }
                }

                if (has_trigger_electron_custom) {
                    events_with_trigger_custom++;
                }

                // Process the event
                processEvent(particles, event_count, stats, chi2pidCustomList,
                            hMomentumCustom, hChi2pidCustom, hDtCustom,
                            hDtVsPCustom, hChi2pidVsPCustom, hBetaVsPCustom,
                            has_trigger_electron_custom);

                if (event_count % 1000 == 0) {
                    cout << "Processed " << event_count << " events" << endl;
                }
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
    cout << "Total electrons (Assigned PID): " << total_electrons_custom << endl;
    cout << "Events with trigger electrons (Assigned PID): " << events_with_trigger_custom << endl;
    cout << "Total trigger electrons (Assigned PID): " << total_trigger_electrons_custom << endl;

    // Save statistics to CSV
    saveStatsToCSV(stats, "output/stats.csv");

    // Save plots
    savePlots("output", chi2pidCustomList, hMomentumCustom, hChi2pidCustom, hDtCustom,
              hDtVsPCustom, hChi2pidVsPCustom, hBetaVsPCustom);

    cout << "\nTotal events processed: " << total_events_processed << endl;
    cout << "Real time: " << timer.RealTime() << " s, CPU time: " << timer.CpuTime() << " s" << endl;

    // Clean up histograms
    for (auto& h : hMomentumCustom) delete h.second;
    for (auto& h : hChi2pidCustom) delete h.second;
    for (auto& h : hDtCustom) delete h.second;
    for (auto& h : hDtVsPCustom) delete h.second;
    for (auto& h : hChi2pidVsPCustom) delete h.second;
    for (auto& h : hBetaVsPCustom) delete h.second;

    return 0;
}