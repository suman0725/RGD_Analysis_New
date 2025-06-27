#ifndef VALIDATION_H
#define VALIDATION_H

#include "ParticleData.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"


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
    /* std::map<int, int> pid_counts;
    std::map<int, int> assigned_pid_counts; */

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

void validateAndCompare(std::vector<ParticleData>& particles, int eventNum, ParticleStats& stats,
                       std::vector<Mismatch>& mismatches, std::vector<float>& chi2pidList,
                       std::vector<float>& chi2pidCustomList,
                       std::map<int, TH1F*>& hMomentumRec, std::map<int, TH1F*>& hMomentumCustom,
                       std::map<int, TH1F*>& hChi2pidRec, std::map<int, TH1F*>& hChi2pidCustom,
                       std::map<int, TH1F*>& hDtRec, std::map<int, TH1F*>& hDtCustom,
                       std::map<int, TH2F*>& hDtVsPRec, std::map<int, TH2F*>& hDtVsPCustom,
                       std::map<int, TH2F*>& hChi2pidVsPRec, std::map<int, TH2F*>& hChi2pidVsPCustom);

// Helper function to determine region from status
Region determineRegion(int status);

#endif