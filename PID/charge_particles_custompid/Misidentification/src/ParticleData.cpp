#include "ParticleData.h"
#include "CCDB/Calibration.h"
#include "CCDB/CalibrationGenerator.h"
#include <cmath>

// Global variables
std::map<std::tuple<int, int, int>, float> tresCache;
std::vector<SamplingFractionParams> sfParams(7);
float htccNpheCut = 2.0f;
const std::vector<std::pair<int, std::vector<int>>> chargedBetaDetectors = {
    {FTOF_DETECTOR, {2, 1, 3}}, {CTOF_DETECTOR, {1}} , {ECAL_DETECTOR, {1, 4, 7}}
};

// Load CCDB parameters
void loadCCDBParams() {
    ccdb::Calibration *calib = ccdb::CalibrationGenerator::CreateCalibration(
        "mysql://clas12reader@clasdb.jlab.org/clas12", 18536, "default");
    if (!calib) { 
        std::cerr << "Failed to create CCDB Calibration object! Using defaults." << std::endl;
        return;
    }

    // Load sampling fraction parameters
    std::vector<std::vector<double>> sfValues;
    if (!calib->GetCalib(sfValues, "/calibration/eb/electron_sf")) {
        cerr << "Failed to get electron_sf! Using defaults." << endl;
    } else {
        for (const auto& row : sfValues) {
            int sector = static_cast<int>(row[0]);
            if (sector < 1 || sector > 6 || row.size() < 11) continue;
            sfParams[sector].sf[0] = row[3]; sfParams[sector].sf[1] = row[4];
            sfParams[sector].sf[2] = row[5]; sfParams[sector].sf[3] = row[6];
            sfParams[sector].sfs[0] = row[7]; sfParams[sector].sfs[1] = row[8];
            sfParams[sector].sfs[2] = row[9]; sfParams[sector].sfs[3] = row[10];
        }
        std::cout << "Loaded sampling fraction parameters for " << sfValues.size() << " sectors from CCDB" << std::endl;
    }

    // Load HTCC nphe cut
    std::vector<std::vector<double>> htccValues;
    if (!calib->GetCalib(htccValues, "/calibration/eb/htcc_matching") || htccValues.empty() || htccValues[0].size() < 7) {
        std::cerr << "Failed to get HTCC nphe cut! Using 2.0" << std::endl;
        htccNpheCut = 2.0f;
    } else {
        htccNpheCut = htccValues[0][6];
    }

    // Load timing resolutions
    std::vector<std::vector<double>> tresValues;
    if (calib->GetCalib(tresValues, "/calibration/ftof/tres")) {
        for (const auto& row : tresValues) {
            int sector = static_cast<int>(row[0]);
            int layer = static_cast<int>(row[1]);
            int component = static_cast<int>(row[2]);
            float tres = row[3];
            tresCache[std::make_tuple(sector, layer, component)] = tres;
        }
        std::cout << "Cached " << tresCache.size() << " tres values from CCDB" << std::endl;
    } else {
        std::cerr << "Failed to load /calibration/ftof/tres from CCDB! No timing resolutions available." << std::endl;
    }

    delete calib;
}

// Sampling fraction helper functions
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



// Load mapping from detector banks to particle indices
IndexMap loadMapByIndex(hipo::bank& fromBank, const char* idxVarName) {
    IndexMap map;
    for (int i = 0; i < fromBank.getRows(); ++i) {
        int iTo = fromBank.getInt(idxVarName, i);
        map[iTo].push_back(i);
    }
    return map;
}

// Extract particle data from HIPO banks
ParticleData getParticleData(int partIdx, hipo::bank& PART, hipo::bank& CHER, hipo::bank& CALO, hipo::bank& SCIN,
                             IndexMap& cherMap, IndexMap& caloMap, IndexMap& scinMap) {
    ParticleData pd;

    pd.charge = PART.getByte("charge", partIdx);
    pd.px = PART.getFloat("px", partIdx);
    pd.py = PART.getFloat("py", partIdx);
    pd.pz = PART.getFloat("pz", partIdx);
    pd.p = std::sqrt(pd.px * pd.px + pd.py * pd.py + pd.pz * pd.pz);
    pd.theta = (pd.p > 0) ? std::acos(pd.pz / pd.p) * 180.0 / M_PI : 0.0f;
    pd.phi = std::atan2(pd.py, pd.px) * 180.0 / M_PI;
    pd.vx = PART.getFloat("vx", partIdx);
    pd.vy = PART.getFloat("vy", partIdx);
    pd.vz = PART.getFloat("vz", partIdx);
    pd.vt = PART.getFloat("vt", partIdx);
    pd.beta = PART.getFloat("beta", partIdx);
    pd.status = PART.getShort("status", partIdx);
    pd.is_trigger = (pd.status < 0);
    pd.start_time = pd.vt;
    pd.pid = PART.getInt("pid", partIdx);
    pd.chi2pid = PART.getFloat("chi2pid", partIdx);

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

// Electron/positron identification
bool isSimpleElectron(const ParticleData& pd) {
    if (pd.charge != -1 && pd.charge != 1) return false;
    if (pd.sector < 1 || pd.sector > 6) return false;
    if (pd.nphe_htcc < htccNpheCut) return false;
    if (pd.energy_pcal < MIN_PCAL_ENERGY) return false;
    float sf = pd.energy_total / pd.p;
    float mean = getSamplingFractionMean(pd.sector, pd.energy_total);
    float sigma = getSamplingFractionSigma(pd.sector, pd.energy_total);
    float nSigma = getSamplingFractionNSigma(sf, mean, sigma);
    return std::abs(nSigma) <= NSIGMA_CUT;
}

// Trigger electron identification
bool isTriggerElectron(const ParticleData& pd) {
    if (pd.status >= 0) return false; // Trigger particles have negative status
    if (!isSimpleElectron(pd) || pd.charge != -1) return false; // Only electrons (not positrons)
    return true;
}

float getTheoryBeta(const ParticleData& p, const float mass){

    return p.p / std::sqrt(p.p * p.p + mass * mass); 
}

// To evaluate the calculated mass for the reconstruted particle 
float getCalculatedMass(const ParticleData& p){
    //float beta = getTheoryBeta(p); 
    float beta = p.beta;
    if (beta <= 0 || beta >= 1) return 99999.0f;
    float mass = (p.p / beta) * std::sqrt(1 - beta * beta); 
    return mass; 
}

// Compute delta T for timing
float computeDeltaT(const ParticleData& p) {
    float delta_t = 99999.0f;
    for (const auto& det : chargedBetaDetectors) {
        for (int layer : det.second) {
            for (const auto& hit : p.hits) {
                if (hit.detector == det.first && hit.layer == layer) {
                    float beta_theory = p.p / std::sqrt(p.p * p.p + PION_MASS * PION_MASS);
                    if (beta_theory > 0) {
                        float vt = hit.time - hit.path / (SPEED_OF_LIGHT * beta_theory);
                        delta_t = vt - p.start_time;
                        return delta_t;
                    }
                }
            }
        }
    }
    return delta_t;
}

// Compute chi2pid for pion hypothesis
float computeChi2pid(const ParticleData& p, const float massHypothesis) {
    float q = 99999.0f;
    float sigma = -1.0f;
    float delta_t = 99999.0f;
    bool found = false;

    for (const auto& det : chargedBetaDetectors) {
        for (int layer : det.second) {
            for (const auto& hit : p.hits) {
                if (hit.detector == det.first && hit.layer == layer) {
                    if (hit.detector == FTOF_DETECTOR) {
                        auto key = make_tuple(hit.sector, hit.layer, hit.component);
                        if (tresCache.count(key)) {
                            sigma = tresCache[key];
                        } else {
                            return 99999.0f; // Skip if resolution not found in CCDB
                        }
                    } else if (hit.detector == CTOF_DETECTOR) {
                        sigma = 0.065f;
                    }else if (hit.detector == ECAL_DETECTOR) {
                        return 99999.0f; // Align with EB's getDetTimingResolution returning 0
                    } else {
                        return 99999.0f;
                    }
                    float beta_theory = p.p / sqrt(p.p * p.p + massHypothesis * massHypothesis);
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
    if (sigma > 0 && found) q = delta_t / sigma;
    return q;
}