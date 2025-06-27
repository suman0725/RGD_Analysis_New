#ifndef PARTICLE_DATA_H
#define PARTICLE_DATA_H

#include <vector>
#include <map>
#include <tuple>
#include "reader.h"

using namespace std;  // Add this to avoid std:: qualifiers

// Constants for detector IDs and cuts
const float MIN_PCAL_ENERGY = 0.06f;
const float NSIGMA_CUT = 5.0f;
const int HTCC_DETECTOR = 15;
const int LTCC_DETECTOR = 16;
const int ECAL_DETECTOR = 7;
const int FTOF_DETECTOR = 12;
const int CTOF_DETECTOR = 4;
const float SPEED_OF_LIGHT = 29.9792458; // cm/ns
const float PION_MASS = 0.13957018; // GeV, for π⁺ hypothesis
const float KAON_MASS = 0.493677; 
const float ELECTRON_MASS = 0.00051099895; // GeV
const float PROTON_MASS = 0.93827208816; // GeV
const float BEAM_ENERGY = 10.6; // GeV, adjust as needed

// Type alias for mapping particle indices to detector bank indices
using IndexMap = map<int, vector<int>>;

// Cache for timing resolutions from CCDB
extern map<tuple<int, int, int>, float> tresCache;

// Structure to hold detector hit information
struct DetectorHit {
    int detector, layer, sector, component;
    float time, path;
};

// Structure to hold particle data
struct ParticleData {
    int charge, sector;
    float px, py, pz, p, theta, phi;
    float nphe_htcc, nphe_ltcc, energy_pcal, energy_total, vt, beta;
    float vx, vy, vz;
    vector<DetectorHit> hits;
    float start_time;
    int status;
    bool is_trigger;
    int pid;
    float chi2pid; 

    ParticleData() : charge(0), sector(-1), px(0.0f), py(0.0f), pz(0.0f), p(0.0f), theta(0.0f), phi(0.0f),
                     nphe_htcc(0.0f), nphe_ltcc(0.0f), energy_pcal(0.0f), energy_total(0.0f), vt(0.0f),
                     beta(-99.0f), vx(0.0f), vy(0.0f), vz(0.0f), start_time(0.0f), status(0),
                     is_trigger(false), pid(0), chi2pid(0.0f) {}
};

// Structure for sampling fraction parameters
struct SamplingFractionParams {
    float sf[4], sfs[4];
};
extern vector<SamplingFractionParams> sfParams;
extern float htccNpheCut;

// List of detectors for beta calculation
extern const vector<pair<int, vector<int>>> chargedBetaDetectors;

// Function declarations
void loadCCDBParams();
float getSamplingFractionMean(int sector, float measuredEnergy);
float getSamplingFractionSigma(int sector, float measuredEnergy);
float getSamplingFractionNSigma(float samplingFraction, float mean, float sigma);
IndexMap loadMapByIndex(hipo::bank& fromBank, const char* idxVarName);
ParticleData getParticleData(int partIdx, hipo::bank& PART, hipo::bank& CHER, hipo::bank& CALO, hipo::bank& SCIN,
                             IndexMap& cherMap, IndexMap& caloMap, IndexMap& scinMap);
bool isSimpleElectron(const ParticleData& pd);
bool isTriggerElectron(const ParticleData& pd);
float computeDeltaT(const ParticleData& p);
float computeChi2pid(const ParticleData& p, const float massHypothesis);
float getTheoryBeta(const ParticleData& p, const float mass);
float getCalculatedMass(const ParticleData& p);

#endif // PARTICLE_DATA_H