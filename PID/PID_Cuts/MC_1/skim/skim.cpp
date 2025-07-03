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

using IndexMap = std::map<int, std::vector<int>>;

const double C = 29.9792458; // Speed of light in cm/ns
const int FTOF_DETECTOR = 12;
const float MASS_PION = 0.139570;
const double DEG_PER_RAD = 57.2957795131; // 180 / π

map<tuple<int, int, int>, float> tresCache;

struct ParticleData {
    int pid, status, charge;
    float p, beta, chi2pid, vz, vt, px, py, pz;
    vector<tuple<int, int, int>> hits;
    int hit_sector, hit_layer, hit_component;
    ParticleData() : hit_sector(-1), hit_layer(-1), hit_component(-1) {}
};

ParticleData getParticleData(int partIdx, hipo::bank& PART, hipo::bank& SCIN, IndexMap& scinMap) {
    ParticleData pd;
    pd.pid = PART.getInt("pid", partIdx);
    pd.px = PART.getFloat("px", partIdx);
    pd.py = PART.getFloat("py", partIdx);
    pd.pz = PART.getFloat("pz", partIdx);
    pd.p = sqrt(pd.px * pd.px + pd.py * pd.py + pd.pz * pd.pz);
    pd.beta = PART.getFloat("beta", partIdx);
    pd.chi2pid = PART.getFloat("chi2pid", partIdx);
    pd.vz = PART.getFloat("vz", partIdx);
    pd.vt = PART.getFloat("vt", partIdx);
    pd.status = PART.getShort("status", partIdx);
    pd.charge = PART.getByte("charge", partIdx);

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

void loadCCDBParams(int runNumber) {
    ccdb::Calibration *calib = ccdb::CalibrationGenerator::CreateCalibration(
        "mysql://clas12reader@clasdb.jlab.org/clas12", runNumber, "default");
    if (!calib) {
        cerr << "Failed to connect to CCDB!" << endl;
        return;
    }

    std::vector<std::vector<double>> tresValues;
    if (calib->GetCalib(tresValues, "/calibration/ftof/tres")) {
        tresCache.clear();
        for (const auto& row : tresValues) {
            int sector = static_cast<int>(row[0]);
            int layer = static_cast<int>(row[1]);
            int component = static_cast<int>(row[2]);
            float tres = row[3];
            tresCache[make_tuple(sector, layer, component)] = tres;
        }
        cout << "Cached " << tresCache.size() << " tres values for run " << runNumber << endl;
    } else {
        cerr << "Failed to load /calibration/ftof/tres for run " << runNumber << endl;
    }
    delete calib;
}

IndexMap loadMapByIndex(hipo::bank& fromBank, const char* idxVarName) {
    IndexMap map;
    for (int i = 0; i < fromBank.getRows(); ++i) {
        int pindex = fromBank.getInt(idxVarName, i);
        if (pindex >= 0) map[pindex].push_back(i);
    }
    return map;
}

float computeDeltaT(ParticleData& p, hipo::bank& SCIN) {
    float delta_t = 99999.0f;
    float mass = MASS_PION;
    float beta_theory = p.p / sqrt(p.p * p.p + mass * mass);
    if (beta_theory <= 0) return delta_t;

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

float computeChi2PidPion(const ParticleData& p, float delta_t, float& sigma_out) {
    if (delta_t >= 99999.0f || p.hit_sector < 1 || p.hit_layer < 1 || p.hit_component < 1) {
        sigma_out = -1.0f;
        return 99999.0f;
    }
    if (!tresCache.count(make_tuple(p.hit_sector, p.hit_layer, p.hit_component))) {
        sigma_out = -1.0f;
        return 99999.0f;
    }
    float sigma = tresCache[make_tuple(p.hit_sector, p.hit_layer, p.hit_component)];
    if (sigma <= 0) {
        sigma_out = -1.0f;
        return 99999.0f;
    }
    float q = delta_t / sigma;
    sigma_out = sigma;
    return q;
}

int findTruthMatch(const ParticleData& pd, hipo::bank& MCPART) {
    float theta_reco = acos(pd.pz / pd.p) * DEG_PER_RAD;
    float phi_reco = atan2(pd.py, pd.px) * DEG_PER_RAD;
    float min_dist = 99999.0f;
    int best_match_idx = -1;

    for (int i = 0; i < MCPART.getRows(); ++i) {
        float px_gen = MCPART.getFloat("px", i);
        float py_gen = MCPART.getFloat("py", i);
        float pz_gen = MCPART.getFloat("pz", i);
        float p_gen = sqrt(px_gen * px_gen + py_gen * py_gen + pz_gen * pz_gen);
        if (p_gen <= 0) continue;

        float theta_gen = acos(pz_gen / p_gen) * DEG_PER_RAD;
        float phi_gen = atan2(py_gen, px_gen) * DEG_PER_RAD;

        float delta_theta = fabs(theta_reco - theta_gen);
        float delta_phi = fabs(phi_reco - phi_gen);
        if (delta_phi > 180.0) delta_phi = 360.0 - delta_phi;

        if (delta_theta < 1.0 && delta_phi < 3.0) {
            float dist = sqrt(delta_theta * delta_theta + delta_phi * delta_phi);
            if (dist < min_dist) {
                min_dist = dist;
                best_match_idx = i;
            }
        }
    }
    return best_match_idx;
}

int main() {
    auto startTime = high_resolution_clock::now();

    ofstream processedFilesStream("processed_mc_hipo_files.txt");
    if (!processedFilesStream.is_open()) {
        cerr << "Error: Could not open processed_mc_hipo_files.txt" << endl;
        return 1;
    }

    TFile* outputFile = new TFile("/w/hallb-scshelf2102/clas12/suman/new_RGD_Analysis/PID/PID_Cuts/MC_1/mc_pion_candidates.root", "RECREATE");
    TTree* treeEBPosPionAssumed = new TTree("EB_pos_pion_assumed_mc", "Positive charge particles assumed as pions (MC)");

    float momentum, beta, vz, dt, recomputed_chi2pid;
    int true_pid, is_matched;
    treeEBPosPionAssumed->Branch("p", &momentum, "p/F");
    treeEBPosPionAssumed->Branch("beta", &beta, "beta/F");
    treeEBPosPionAssumed->Branch("vz", &vz, "vz/F");
    treeEBPosPionAssumed->Branch("dt", &dt, "dt/F");
    treeEBPosPionAssumed->Branch("recomputed_chi2pid", &recomputed_chi2pid, "recomputed_chi2pid/F");
    treeEBPosPionAssumed->Branch("true_pid", &true_pid, "true_pid/I");
    treeEBPosPionAssumed->Branch("is_matched", &is_matched, "is_matched/I");

    string hipoFile = "/w/hallb-scshelf2102/clas12/suman/new_RGD_Analysis/PID/PID_Cuts/MC_1/out.hipo";
    hipo::reader reader;
    reader.open(hipoFile.c_str());
    hipo::dictionary factory;
    reader.readDictionary(factory);
    hipo::event event;
    hipo::bank RUN(factory.getSchema("RUN::config"));
    hipo::bank PART(factory.getSchema("REC::Particle"));
    hipo::bank SCIN(factory.getSchema("REC::Scintillator"));
    hipo::bank MCPART(factory.getSchema("MC::Particle"));

    int runNumber = -1;
    if (reader.gotoEvent(0)) {
        reader.read(event);
        event.getStructure(RUN);
        runNumber = RUN.getInt("run", 0);
        loadCCDBParams(runNumber);
        processedFilesStream << hipoFile << endl;
        cout << "Processing MC file: " << hipoFile << endl;
    } else {
        cout << "Error: Empty MC file: " << hipoFile << endl;
        return 1;
    }

    reader.rewind();

    int event_count = 0;
    int events_with_trigger = 0;
    int maxEvents = 300000; // Updated to match 300,000 events in log
    int electron_count = 0;
    int positron_count = 0;
    int trigger_electron_count = 0;
    int particle_count_charge = 0;
    auto lastUpdateTime = startTime;

    while (reader.next() /* && event_count < maxEvents */) {
        event_count++;
        reader.read(event);
        event.getStructure(PART);
        event.getStructure(SCIN);
        event.getStructure(MCPART);
        IndexMap scinMap = loadMapByIndex(SCIN, "pindex");

        vector<ParticleData> particles(PART.getRows());
        for (int i = 0; i < PART.getRows(); ++i) {
            particles[i] = getParticleData(i, PART, SCIN, scinMap);
        }

        bool hasTriggerElectron = false;
        for (const auto& pd : particles) {
            if (pd.pid == 11) {
                electron_count++;
                if (abs(pd.chi2pid) < 5 && pd.vz >= -10.56 && pd.vz <= 5 && abs(pd.status) / 1000 == 2 && pd.status < 0) {
                    trigger_electron_count++;
                    hasTriggerElectron = true;
                }
            } else if (pd.pid == -11) {
                positron_count++;
            }
        }

        if (hasTriggerElectron) {
            events_with_trigger++;
            for (auto& pd : particles) {
                if (abs(pd.status) / 1000 != 2) continue;

                momentum = pd.p;
                beta = pd.beta;
                vz = pd.vz;
                dt = computeDeltaT(pd, SCIN);
                float sigma;
                recomputed_chi2pid = computeChi2PidPion(pd, dt, sigma);

                if (pd.charge > 0) {
                    int mc_idx = findTruthMatch(pd, MCPART);
                    is_matched = (mc_idx >= 0) ? 1 : 0;
                    true_pid = (mc_idx >= 0) ? MCPART.getInt("pid", mc_idx) : -999;

                    treeEBPosPionAssumed->Fill();
                    particle_count_charge++;
                }
            }
        }

       /*  auto currentTime = high_resolution_clock::now();
        auto elapsedSinceLastUpdate = duration_cast<seconds>(currentTime - lastUpdateTime).count();
        if (elapsedSinceLastUpdate >= 10) { // Update every 10 seconds
            auto elapsed = duration_cast<seconds>(currentTime - startTime).count();
            double eventsPerSecond = event_count / (elapsed + 1e-6);
            int remainingEvents = maxEvents - event_count;
            double remainingSeconds = (remainingEvents > 0) ? remainingEvents / eventsPerSecond : 0;
            int hours = remainingSeconds / 3600;
            int minutes = (remainingSeconds - hours * 3600) / 60;
            int seconds = remainingSeconds - hours * 3600 - minutes * 60;
            cout << "Estimated time remaining: " << hours << "h " << minutes << "m " << seconds << "s" << endl;
            lastUpdateTime = currentTime;
        } */
    }

    auto endTime = high_resolution_clock::now();
    auto totalElapsed = duration_cast<seconds>(endTime - startTime).count();
    cout << "Total events processed: " << event_count << endl;
    cout << "Total particles saved in charge-based tree: " << particle_count_charge << endl;
    cout << "Total electron count: " << electron_count << endl;
    cout << "Total positron count: " << positron_count << endl;
    cout << "Total trigger electron count: " << trigger_electron_count << endl;
    cout << "Events with trigger electrons: " << events_with_trigger << endl;
    cout << "Total processing time: " << totalElapsed / 3600 << "h "
         << (totalElapsed % 3600) / 60 << "m " << (totalElapsed % 60) << "s" << endl;

    outputFile->Write();
    outputFile->Close();
    delete outputFile;
    processedFilesStream.close();

    return 0;
}