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
#include <limits>
#include "reader.h"
#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include <TStopwatch.h>
#include "CCDB/Calibration.h"
#include "CCDB/CalibrationGenerator.h"
#include "TH1F.h"

using namespace std;
namespace fs = std::filesystem;

using IndexMap = std::map<int, std::vector<int>>;
std::map<std::tuple<int, int, int>, float> tresCache;

const float MIN_PCAL_ENERGY = 0.06f;
const float NSIGMA_CUT = 5.0f;
const int HTCC_DETECTOR = 15;
const int LTCC_DETECTOR = 18; // Added for LTCC handling
const int ECAL_DETECTOR = 7;
const int FTOF_DETECTOR = 12;
const int CTOF_DETECTOR = 4;
const float SPEED_OF_LIGHT = 29.9792458; // cm/ns
const float PION_MASS = 0.13957018; // GeV, for π⁺ hypothesis
const std::vector<std::pair<int, std::vector<int>>> chargedBetaDetectors = {
    {FTOF_DETECTOR, {2, 1, 3}}, {CTOF_DETECTOR, {1}}
};

struct DetectorHit {
    int detector, layer, sector, component;
    float time, path;
};

struct ParticleData {
    int charge, sector;
    float p, nphe_htcc, nphe_ltcc, energy_pcal, energy_total, vt, beta;
    std::vector<DetectorHit> hits;
    float start_time;
    int status;
    float vz_tele;
    bool is_trigger;

    ParticleData() : charge(0), p(0.0f), sector(-1), nphe_htcc(0.0f), nphe_ltcc(0.0f), energy_pcal(0.0f),
                     energy_total(0.0f), vt(0.0f), beta(-99.0f), start_time(0.0f), status(0), vz_tele(0.0f), is_trigger(false) {}
};

struct SamplingFractionParams {
    float sf[4], sfs[4];
};
std::vector<SamplingFractionParams> sfParams(7);
float htccNpheCut = 2.0f;

void loadCCDBParams() {
    ccdb::Calibration *calib = ccdb::CalibrationGenerator::CreateCalibration(
        "mysql://clas12reader@clasdb.jlab.org/clas12", 18536, "default");
    if (!calib) { 
        cerr << "Failed to create CCDB Calibration object! Using defaults." << endl;
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
        cerr << "Failed to load /calibration/ftof/tres from CCDB! No timing resolutions available." << endl;
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
    for (int i = 0; i < fromBank.getRows(); ++i) {
        int iTo = fromBank.getInt(idxVarName, i);
        map[iTo].push_back(i);
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
    float sigma = getSamplingFractionSigma(pd.sector,  pd.energy_total);
    float nSigma = getSamplingFractionNSigma(sf, mean, sigma);
    return abs(nSigma) <= NSIGMA_CUT;
}

bool isTriggerElectron(const ParticleData& pd) {
    if (pd.status >= 0) return false; // Trigger particles have negative status
    if (!isSimpleElectron(pd) || pd.charge != -1) return false; // Only electrons (not positrons)
    float sf = pd.energy_total / pd.p;
    float mean = getSamplingFractionMean(pd.sector, pd.energy_total);
    float sigma = getSamplingFractionSigma(pd.sector,  pd.energy_total);
    float chi2pid = getSamplingFractionNSigma(sf, mean, sigma);; // As per your example
    //return (chi2pid < 5 && pd.vz_tele >= -20 && pd.vz_tele <= 5 && abs(pd.status) / 1000 == 2);
}

float computeDeltaT(const ParticleData& p) {
    float delta_t = 99999.0f;
    for (const auto& det : chargedBetaDetectors) {
        for (int layer : det.second) {
            for (const auto& hit : p.hits) {
                if (hit.detector == det.first && hit.layer == layer) {
                    float beta_theory = p.p / sqrt(p.p * p.p + PION_MASS * PION_MASS);
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

float computeChi2pid(const ParticleData& p) {
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
                    } else {
                        return 99999.0f;
                    }
                    float beta_theory = p.p / sqrt(p.p * p.p + PION_MASS * PION_MASS);
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

void processEvent(std::vector<ParticleData>& particles, TH2F* hDtVsP, TH2F* hChi2pidVsP,
                  std::vector<TH1F*>& hDtSlices, std::vector<TH1F*>& hChi2pidSlices) {
    for (const auto& p : particles) {
        if (p.charge <= 0) continue; // Only positive hadrons
        if (isSimpleElectron(p)) continue; // Skip electrons/positrons
        if ((abs(p.status) / 2000) != 1) continue;

        float dt = computeDeltaT(p);
        float chi2pid = computeChi2pid(p);
        if (dt != 99999.0f) hDtVsP->Fill(p.p, dt);
        if (chi2pid != 99999.0f) hChi2pidVsP->Fill(p.p, chi2pid);

        int sliceIdx = static_cast<int>(p.p / 0.3f); // 0.3 GeV bins
        if (sliceIdx >= 0 && sliceIdx < hDtSlices.size()) {
            if (dt != 99999.0f) hDtSlices[sliceIdx]->Fill(dt);
            if (chi2pid != 99999.0f) hChi2pidSlices[sliceIdx]->Fill(chi2pid);
        }
    }
}
int main() {
    TStopwatch timer;
    timer.Start();

    fs::create_directory("output_copy");
    loadCCDBParams();

    TH2F* hDtVsP = new TH2F("hDtVsP", "#DeltaT vs P (Pion Hypothesis);Momentum (GeV);#DeltaT (ns)",
                            200, 0, 10, 200, -5, 5);
    TH2F* hChi2pidVsP = new TH2F("hChi2pidVsP", "chi2pid vs P (Pion Hypothesis);Momentum (GeV);chi2pid",
                                 200, 0, 10, 200, -15, 15);

    // Define momentum slices: 0 to 10 GeV with 0.3 GeV bins
    const int nSlices = static_cast<int>(10.0 / 0.3) + 1; // 34 slices: 0-0.3, 0.3-0.6, ..., 9.9-10.0
    std::vector<TH1F*> hDtSlices(nSlices);
    std::vector<TH1F*> hChi2pidSlices(nSlices);
    for (int i = 0; i < nSlices; ++i) {
        float pMin = i * 0.3f;
        float pMax = (i + 1) * 0.3f;
        string dtName = "hDt_p" + to_string(pMin) + "_" + to_string(pMax);
        string chiName = "hChi2pid_p" + to_string(pMin) + "_" + to_string(pMax);
        hDtSlices[i] = new TH1F(dtName.c_str(), ("#DeltaT (p: " + to_string(pMin) + "-" + to_string(pMax) + " GeV);#DeltaT (ns);Counts").c_str(),
                                100, -5, 5);
        hChi2pidSlices[i] = new TH1F(chiName.c_str(), ("chi2pid (p: " + to_string(pMin) + "-" + to_string(pMax) + " GeV);chi2pid;Counts").c_str(),
                                     100, -15, 15);
    }

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

    int event_count = 0;
    int events_with_trigger = 0;
    int maxEvents = 1000000;
    int electron_count = 0;
    int positron_count = 0; 
    int trigger_electron_count = 0; 

    for (const auto& dir : directories) {
        for (const auto& entry : fs::directory_iterator(dir)) {
            if (entry.path().extension() != ".hipo") continue;
            hipo::reader reader(entry.path().c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);

            while (reader.next() && event_count < maxEvents) {
                event_count++;
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

                vector<ParticleData> particles(PART.getRows());
                for (int i = 0; i < PART.getRows(); ++i) {
                    particles[i] = getParticleData(i, PART, CHER, CALO, SCIN, cherMap, caloMap, scinMap);
                }

                bool hadTrigger = false;
               for (const auto& p : particles) {
                    if (isSimpleElectron(p)) {
                        if (p.charge == -1) electron_count++;
                        else if (p.charge == 1) positron_count++;
                    }
                    if (isTriggerElectron(p)) {
                        // Apply additional cuts for trigger electrons
                        float sf = p.energy_total / p.p;
                        float mean = getSamplingFractionMean(p.sector, p.p);
                        float sigma = getSamplingFractionSigma(p.sector, p.p);
                        float chi2pid = getSamplingFractionNSigma(sf, mean, sigma);
                    
                        if (abs(chi2pid) < 5 && p.vz_tele >= -20 && p.vz_tele <= 5 && abs(p.status) / 1000 == 2) {
                            trigger_electron_count++;
                            hadTrigger = true;
                        }
                    }
                }
                if (hadTrigger) {
                    processEvent(particles, hDtVsP, hChi2pidVsP, hDtSlices, hChi2pidSlices);
                }
            }
            if (event_count >= maxEvents) break;
        }
        if (event_count >= maxEvents) break;
    }

           // Original 2D plots
    TCanvas* c1 = new TCanvas("c1", "Delta T vs P", 800, 600);
    hDtVsP->SetStats(0);
    hDtVsP->Draw("COLZ");
    gPad->SetLogz();
    c1->SaveAs("output_copy/delta_t_vs_p_pion_FR.pdf");

    TCanvas* c2 = new TCanvas("c2", "Chi2 PID vs P", 800, 600);
    hChi2pidVsP->SetStats(0);
    hChi2pidVsP->Draw("COLZ");
    gPad->SetLogz();
    c2->SaveAs("output_copy/chi2pid_vs_p_pion_FR.pdf");

    // Single PDF for all ΔT slices, one canvas per bin
    string dtPdfName = "output_copy/delta_t_all_bins_pion.pdf";
    TCanvas* cDt = new TCanvas("cDt", "Delta T", 800, 600);
    cDt->Print((dtPdfName + "[").c_str()); // Open the PDF

    for (int i = 0; i < nSlices; ++i) {
        if (hDtSlices[i]->GetEntries() == 0) continue; // Skip empty histograms

        cDt->Clear(); // Clear the canvas for the new histogram 
        hDtSlices[i]->Scale(1.0 / hDtSlices[i]->Integral()); 
        hDtSlices[i]->SetStats(0);
        hDtSlices[i]->SetLineColor(kBlue);
        hChi2pidSlices[i]->SetMarkerStyle(20);  // Circular marker
        hChi2pidSlices[i]->SetMarkerSize(1.0); 
        hDtSlices[i]->Draw("p");
        //gPad->SetLogy(); // Optional: log scale for better visibility
        cDt->Print(dtPdfName.c_str()); // Save the current canvas as a new page
    }

    cDt->Print((dtPdfName + "]").c_str()); // Close the PDF

    // Single PDF for all χ²_PID slices, one canvas per bin
    string chi2PdfName = "output_copy/chi2pid_all_bins_pion.pdf";
    TCanvas* cChi2 = new TCanvas("cChi2", "Chi2 PID", 800, 600);
    cChi2->Print((chi2PdfName + "[").c_str()); // Open the PDF

    for (int i = 0; i < nSlices; ++i) {
        if (hChi2pidSlices[i]->GetEntries() == 0) continue; // Skip empty histograms

        cChi2->Clear(); // Clear the canvas for the new histogram
        hChi2pidSlices[i]->Scale(1.0 / hChi2pidSlices[i]->Integral());
        hChi2pidSlices[i]->SetStats(0);
        hChi2pidSlices[i]->SetLineColor(kBlue);
        hChi2pidSlices[i]->SetMarkerStyle(20);  // Circular marker
        hChi2pidSlices[i]->SetMarkerSize(1.0);  // Normal size
        hChi2pidSlices[i]->Draw("p");           // Draw as points
        //gPad->SetLogy(); // Optional: log scale for better visibility
        cChi2->Print(chi2PdfName.c_str()); // Save the current canvas as a new page
    }

    cChi2->Print((chi2PdfName + "]").c_str()); // Close the PDF

    timer.Stop();
    cout << "Total events processed: " << event_count << endl;
    cout << "Total electrons: " << electron_count << endl;

    cout << "Events with trigger electrons: " << trigger_electron_count << endl;
    cout << "Real time: " << timer.RealTime() << " s" << endl;

    // Cleanup
    delete hDtVsP;
    delete hChi2pidVsP;
    delete c1;
    delete c2;
    delete cDt;
    delete cChi2;
    for (auto h : hDtSlices) delete h;
    for (auto h : hChi2pidSlices) delete h;

    return 0;
}