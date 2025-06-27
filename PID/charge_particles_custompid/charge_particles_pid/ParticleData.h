#ifndef PARTICLE_DATA_H
#define PARTICLE_DATA_H

#include <map>
#include <vector>
#include <set>
#include <tuple>

namespace fs = std::filesystem;

// Define IndexMap for mapping particle indices to detector hits
using IndexMap = std::map<int, std::vector<int>>;

// Cache for FTOF timing resolutions
extern std::map<std::tuple<int, int, int>, float> tresCache;

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
extern std::map<int, float> PDG_MASS;

// Detector priorities for timing-based PID (only for charged particles)
extern const std::vector<std::pair<int, std::vector<int>>> chargedBetaDetectors;

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

    ParticleData() : charge(0), p(0.0f), sector(-1), nphe_htcc(0.0f), nphe_ltcc(0.0f),
                     energy_pcal(0.0f), energy_total(0.0f), pid(0), vt(0.0f), beta(-99.0f),
                     hits(), is_trigger(false), assigned_pid(0), start_time(0.0f), status(0),
                     chi2pid(0.0f), chi2pid_custom(99999.0f) {}
};

// Sampling Fraction Parameters for electron identification
struct SamplingFractionParams {
    float sf[4];  // Mean sampling fraction parameters
    float sfs[4]; // Sigma sampling fraction parameters
};
extern std::vector<SamplingFractionParams> sfParams;
extern float htccNpheCut;
extern float ltccNpheCut;

#endif