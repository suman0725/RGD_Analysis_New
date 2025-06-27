#include "Validation.h"
#include <iostream>
#include "TLegend.h"
#include "TStyle.h"

Region determineRegion(int status) {
    int abs_status = std::abs(status);
    if (abs_status / 2000 == 1) return Region::Forward;
    else if (abs_status / 4000 == 1) return Region::Central;
    else if (abs_status / 8000 == 1) return Region::Band;
    else return Region::Unknown;
}

void validateAndCompare(vector<ParticleData>& particles, int eventNum, ParticleStats& stats,
                       vector<Mismatch>& mismatches, vector<float>& chi2pidList,
                       vector<float>& chi2pidCustomList,
                       map<int, TH1F*>& hMomentumRec, map<int, TH1F*>& hMomentumCustom,
                       map<int, TH1F*>& hChi2pidRec, map<int, TH1F*>& hChi2pidCustom,
                       map<int, TH1F*>& hDtRec, map<int, TH1F*>& hDtCustom,
                       map<int, TH2F*>& hDtVsPRec, map<int, TH2F*>& hDtVsPCustom,
                       map<int, TH2F*>& hChi2pidVsPRec, map<int, TH2F*>& hChi2pidVsPCustom) {
    int valid_chi2pid_count = 0;
    int valid_chi2pid_custom_count = 0;
    int charged_particles = 0;

    for (size_t j = 0; j < particles.size(); ++j) {
        const auto& p = particles[j];

        // Skip neutral particles
        if (p.charge == 0) continue;

        charged_particles++;

        // Determine the region based on status
        Region region = determineRegion(p.status);

        // Update overall and region-specific PID counts
        stats.pid_counts[p.pid]++;
        stats.pid_counts_by_region[region][p.pid]++;
        stats.assigned_pid_counts[p.assigned_pid]++;
        stats.assigned_pid_counts_by_region[region][p.assigned_pid]++;

        // Fill momentum histograms
        if (hMomentumRec.count(p.pid)) {
            hMomentumRec[p.pid]->Fill(p.p);
        }
        if (hMomentumCustom.count(p.assigned_pid)) {
            hMomentumCustom[p.assigned_pid]->Fill(p.p);
        }

        // Fill chi2pid histograms
        if (p.chi2pid != 9999.0f) {
            valid_chi2pid_count++;
            chi2pidList.push_back(p.chi2pid);
            if (hChi2pidRec.count(p.pid)) {
                hChi2pidRec[p.pid]->Fill(p.chi2pid);
            }
            if (hChi2pidVsPRec.count(p.pid)) {
                hChi2pidVsPRec[p.pid]->Fill(p.p, p.chi2pid);
            }
        } else {
            stats.total_invalid_chi2pid++;
            if (p.pid != 0) stats.total_invalid_chi2pid_nonzero_pid++;
        }

        if (p.chi2pid_custom != 99999.0f) {
            valid_chi2pid_custom_count++;
            stats.valid_chi2pid_custom++;
            chi2pidCustomList.push_back(p.chi2pid_custom);
            if (hChi2pidCustom.count(p.assigned_pid)) {
                hChi2pidCustom[p.assigned_pid]->Fill(p.chi2pid_custom);
            }
            if (hChi2pidVsPCustom.count(p.assigned_pid)) {
                hChi2pidVsPCustom[p.assigned_pid]->Fill(p.p, p.chi2pid_custom);
            }
            if (p.assigned_pid == 0) {
                stats.valid_chi2pid_custom_pid_zero++;
            } else {
                stats.valid_chi2pid_custom_nonzero_pid++;
            }
        } else {
            stats.total_invalid_chi2pid_custom++;
            if (p.assigned_pid != 0) {
                stats.total_invalid_chi2pid_custom_nonzero_pid++;
            } else {
                stats.invalid_chi2pid_custom_pid_zero++;
            }
        }

        // Compute and fill Delta T histograms
        float dt_rec = computeDeltaT(p, p.pid);
        if (dt_rec != 99999.0f && hDtRec.count(p.pid)) {
            hDtRec[p.pid]->Fill(dt_rec);
            if (hDtVsPRec.count(p.pid)) {
                hDtVsPRec[p.pid]->Fill(p.p, dt_rec);
            }
        }

        float dt_custom = computeDeltaT(p, p.assigned_pid);
        if (dt_custom != 99999.0f && hDtCustom.count(p.assigned_pid)) {
            hDtCustom[p.assigned_pid]->Fill(dt_custom);
            if (hDtVsPCustom.count(p.assigned_pid)) {
                hDtVsPCustom[p.assigned_pid]->Fill(p.p, dt_custom);
            }
        }

        if (p.assigned_pid != p.pid) {
            bool has_ftof = hasHit(p, FTOF_DETECTOR);
            bool has_ctof = hasHit(p, CTOF_DETECTOR);
            mismatches.push_back({eventNum, j, p.pid, p.assigned_pid, (int)p.charge, p.p,
                                 p.nphe_htcc, p.nphe_ltcc, p.status, p.chi2pid,
                                 p.chi2pid_custom, has_ftof, has_ctof});
        }
    }

    stats.valid_chi2pid += valid_chi2pid_count;

    cout << "Event " << eventNum << ": Total particles=" << particles.size()
         << ", Charged particles=" << charged_particles
         << ", Valid chi2pid (REC)=" << valid_chi2pid_count
         << ", Valid chi2pid_custom (Computed for Assigned)=" << valid_chi2pid_custom_count << endl;
}