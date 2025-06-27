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
#include "TH1F.h"
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

// Constants
const float MIN_PCAL_ENERGY = 0.06f; // GeV
const float NSIGMA_CUT = 5.0f;
const int HTCC_DETECTOR = 15; // HTCC detector ID
const int ECAL_DETECTOR = 7;  // ECAL detector ID 
const int ELECTRON_PID = 11;  // PID for electrons in REC::Particle

// Map type for bank indices
typedef std::map<int, std::vector<int>> IndexMap;

// Particle data structure
struct ParticleData {
    int charge;
    float p;          // Momentum
    int sector;       // ECAL sector (1-6)
    float nphe;       // HTCC photoelectrons
    float energy_pcal; // PCAL energy
    float energy_total; // PCAL + ECIN + ECOUT
    int pid;          // PID from REC::Particle
};

// CCDB parameters for sampling fraction (per sector)
struct SamplingFractionParams {
    float sf[4];  // sf1, sf2, sf3, sf4 for mean
    float sfs[4]; // sfs1, sfs2, sfs3, sfs4 for sigma
};
std::vector<SamplingFractionParams> sfParams(7); // Index 0 unused, 1-6 for sectors
float htccNpheCut = 0.0f; // Global HTCC nphe cut from CCDB

// Utility function to build index map
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

// Load CCDB parameters for run 18536
void loadCCDBParams() {
    ccdb::Calibration *calib = ccdb::CalibrationGenerator::CreateCalibration("mysql://clas12reader@clasdb.jlab.org/clas12", 18536, "default");
    if (!calib) {
        std::cerr << "Failed to create CCDB Calibration object!" << std::endl;
        exit(1);
    }

    // Fetch electron_sf parameters
    std::vector<std::vector<double>> sfValues;
    if (!calib->GetCalib(sfValues, "/calibration/eb/electron_sf")) {
        std::cerr << "Failed to get electron_sf data from /calibration/eb/electron_sf for run 18536! Using hardcoded defaults." << std::endl;
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
            if (sector < 1 || sector > 6) continue;
            if (row.size() < 11) continue;
            sfParams[sector].sf[0] = row[3];  // sf1
            sfParams[sector].sf[1] = row[4];  // sf2
            sfParams[sector].sf[2] = row[5];  // sf3
            sfParams[sector].sf[3] = row[6];  // sf4
            sfParams[sector].sfs[0] = row[7]; // sfs1
            sfParams[sector].sfs[1] = row[8]; // sfs2
            sfParams[sector].sfs[2] = row[9]; // sfs3
            sfParams[sector].sfs[3] = row[10]; // sfs4
        }
    }

    // Fetch HTCC nphe cut from /calibration/eb/htcc_matching
    std::vector<std::vector<double>> htccValues;
    if (!calib->GetCalib(htccValues, "/calibration/eb/htcc_matching")) {
        std::cerr << "Failed to get HTCC nphe cut from /calibration/eb/htcc_matching! Using default 2.0" << std::endl;
        htccNpheCut = 2.0f;
    } else if (htccValues.empty() || htccValues[0].size() < 7) {
        std::cerr << "Invalid HTCC data format! Using default 2.0" << std::endl;
        htccNpheCut = 2.0f;
    } else {
        htccNpheCut = htccValues[0][6]; // nphe is 7th column (index 6)
    }

    delete calib;
}

// Calculate sampling fraction mean
float getSamplingFractionMean(int sector, float measuredEnergy) {
    if (sector < 1 || sector > 6) return 0.0f;
    const auto& p = sfParams[sector].sf;
    return p[0] * (p[1] + p[2] / measuredEnergy + p[3] * std::pow(measuredEnergy, -2));
}

// Calculate sampling fraction sigma
float getSamplingFractionSigma(int sector, float measuredEnergy) {
    if (sector < 1 || sector > 6) return 0.0f;
    const auto& p = sfParams[sector].sfs;
    return p[0] * (p[1] + p[2] / measuredEnergy + p[3] * std::pow(measuredEnergy, -2));
}

// Calculate nSigma
float getSamplingFractionNSigma(float samplingFraction, float mean, float sigma) {
    if (sigma == 0) return 0.0f;
    return (samplingFraction - mean) / sigma;
}

// Extract particle data
ParticleData getParticleData(int partIdx, hipo::bank& PART, hipo::bank& CHER, hipo::bank& CALO,
                             IndexMap& cherMap, IndexMap& caloMap) {
    ParticleData pd = {0, 0.0f, -1, 0.0f, 0.0f, 0.0f, 0};
    pd.charge = PART.getByte("charge", partIdx);
    pd.pid = PART.getInt("pid", partIdx);
    float px = PART.getFloat("px", partIdx);
    float py = PART.getFloat("py", partIdx);
    float pz = PART.getFloat("pz", partIdx);
    pd.p = std::sqrt(px * px + py * py + pz * pz);

    if (caloMap.find(partIdx) != caloMap.end()) {
        for (int iCalRow : caloMap[partIdx]) {
            int det_cal = CALO.getByte("detector", iCalRow);
            if (det_cal != ECAL_DETECTOR) continue;
            pd.sector = CALO.getByte("sector", iCalRow);
            int layer_cal = CALO.getByte("layer", iCalRow);
            if (layer_cal == 1) pd.energy_pcal = CALO.getFloat("energy", iCalRow);
            else if (layer_cal == 4) pd.energy_total += CALO.getFloat("energy", iCalRow);
            else if (layer_cal == 7) pd.energy_total += CALO.getFloat("energy", iCalRow);
        }
        pd.energy_total += pd.energy_pcal;
    }

    if (cherMap.find(partIdx) != cherMap.end()) {
        for (int iCherRow : cherMap[partIdx]) {
            if (CHER.getByte("detector", iCherRow) == HTCC_DETECTOR) {
                pd.nphe = CHER.getFloat("nphe", iCherRow);
                break;
            }
        }
    }

    return pd;
}

// Check if particle is a simple electron
bool isSimpleElectron(const ParticleData& pd) {
    if (pd.charge != -1) return false;
    if (pd.sector < 1 || pd.sector > 6) return false;
    if (pd.nphe < htccNpheCut) return false;
    if (pd.energy_pcal < MIN_PCAL_ENERGY) return false;

    float sampling_fraction = pd.energy_total / pd.p;
    float mean = getSamplingFractionMean(pd.sector, pd.energy_total);
    float sigma = getSamplingFractionSigma(pd.sector, pd.energy_total);
    float sfNSigma = getSamplingFractionNSigma(sampling_fraction, mean, sigma);
    if (std::abs(sfNSigma) > NSIGMA_CUT) return false;

    return true;
}

int main() {
    TStopwatch timer;
    timer.Start();

    // Read directories from file
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

    // Histograms for simple electrons (red)
    TH1F* h1_nphe_tele = new TH1F("h1_nphe_tele", "NPHE; nphe; counts", 200, 0, 60);
    TH1F* h1_nphe_tele_1 = new TH1F("h1_nphe_tele_1", "NPHE (0-10); nphe; counts", 200, 0, 10);
    TH1F* h1_charge_tele = new TH1F("h1_charge_tele", "Charge; charge; counts", 200, -2, 2);
    TH1F* h1_status_tele = new TH1F("h1_status_tele", "Status; status; counts", 200, -4000, -2000);
    TH1F* h1_vz_tele = new TH1F("h1_vz_tele", "Vz; vz (cm); counts", 200, -25, 10);
    TH1F* h1_chi2pid_tele = new TH1F("h1_chi2pid_tele", "Chi2PID; chi2pid; counts", 200, -10, 10);
    TH1F* h1_energy_pcal_tele = new TH1F("h1_energy_pcal_tele", "PCAL Energy; energy_pcal; counts", 200, 0, 2);
    TH1F* h1_energy_pcal_tele_1 = new TH1F("h1_energy_pcal_tele_1", "PCAL Energy (0-0.3 GeV); energy_pcal; counts", 200, 0, 0.3);
    TH1F* h1_sampling_fraction_tele = new TH1F("h1_sampling_fraction_tele", "Sampling Fraction; sampling_fraction; counts", 200, 0, 0.6);

    // Histograms for PID 11 electrons (blue)
    TH1F* h1_nphe_pid11 = new TH1F("h1_nphe_pid11", "NPHE; nphe; counts", 200, 0, 60);
    TH1F* h1_nphe_pid11_1 = new TH1F("h1_nphe_pid11_1", "NPHE (0-10); nphe; counts", 200, 0, 10);
    TH1F* h1_charge_pid11 = new TH1F("h1_charge_pid11", "Charge; charge; counts", 200, -2, 2);
    TH1F* h1_status_pid11 = new TH1F("h1_status_pid11", "Status; status; counts", 200, -4000, -2000);
    TH1F* h1_vz_pid11 = new TH1F("h1_vz_pid11", "Vz; vz (cm); counts", 200, -25, 10);
    TH1F* h1_chi2pid_pid11 = new TH1F("h1_chi2pid_pid11", "Chi2PID; chi2pid; counts", 200, -10, 10);
    TH1F* h1_energy_pcal_pid11 = new TH1F("h1_energy_pcal_pid11", "PCAL Energy; energy_pcal; counts", 200, 0, 2);
    TH1F* h1_energy_pcal_pid11_1 = new TH1F("h1_energy_pcal_pid11_1", "PCAL Energy (0-0.3 GeV); energy_pcal; counts", 200, 0, 0.3);
    TH1F* h1_sampling_fraction_pid11 = new TH1F("h1_sampling_fraction_pid11", "Sampling Fraction; sampling_fraction; counts", 200, 0, 0.6);

    // Set colors and transparency
    h1_nphe_tele->SetLineColor(kRed); h1_nphe_pid11->SetLineColor(kBlue);
    h1_nphe_tele_1->SetLineColor(kRed); h1_nphe_pid11_1->SetLineColor(kBlue);
    h1_charge_tele->SetLineColor(kRed); h1_charge_pid11->SetLineColor(kBlue);
    h1_status_tele->SetLineColor(kRed); h1_status_pid11->SetLineColor(kBlue);
    h1_vz_tele->SetLineColor(kRed); h1_vz_pid11->SetLineColor(kBlue);
    h1_chi2pid_tele->SetLineColor(kRed); h1_chi2pid_pid11->SetLineColor(kBlue);
    h1_energy_pcal_tele->SetLineColor(kRed); h1_energy_pcal_pid11->SetLineColor(kBlue);
    h1_energy_pcal_tele_1->SetLineColor(kRed); h1_energy_pcal_pid11_1->SetLineColor(kBlue);
    h1_sampling_fraction_tele->SetLineColor(kRed); h1_sampling_fraction_pid11->SetLineColor(kBlue);
    // Add transparency (0.5 alpha)
    h1_nphe_tele->SetLineWidth(2); h1_nphe_pid11->SetLineWidth(2);
    h1_nphe_tele_1->SetLineWidth(2); h1_nphe_pid11_1->SetLineWidth(2);
    h1_charge_tele->SetLineWidth(2); h1_charge_pid11->SetLineWidth(2);
    h1_status_tele->SetLineWidth(2); h1_status_pid11->SetLineWidth(2);
    h1_vz_tele->SetLineWidth(2); h1_vz_pid11->SetLineWidth(2);
    h1_chi2pid_tele->SetLineWidth(2); h1_chi2pid_pid11->SetLineWidth(2);
    h1_energy_pcal_tele->SetLineWidth(2); h1_energy_pcal_pid11->SetLineWidth(2);
    h1_energy_pcal_tele_1->SetLineWidth(2); h1_energy_pcal_pid11_1->SetLineWidth(2);
    h1_sampling_fraction_tele->SetLineWidth(2); h1_sampling_fraction_pid11->SetLineWidth(2);

    int total_events_processed = 0;
    int total_electron_number_simple = 0;
    int total_electron_number_pid11 = 0;
    const int maxEvents = 100000;
    int event_count = 0;

    // Debug variables
    int mismatchCount = 0; // Count particles that are PID 11 but not simple, or vice versa

    // Load CCDB parameters once
    loadCCDBParams();

    for (const auto& dir : directories) {
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
            cout << "  Opening file: " << file << endl;
            hipo::reader reader;
            reader.open(file.c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);

            while (reader.next()) {
                if (event_count >= maxEvents) break;
                event_count++;
                total_events_processed++;

                hipo::event event;
                reader.read(event);

                hipo::bank PART(factory.getSchema("REC::Particle"));
                hipo::bank CHER(factory.getSchema("REC::Cherenkov"));
                hipo::bank CALO(factory.getSchema("REC::Calorimeter"));
                event.getStructure(PART);
                event.getStructure(CHER);
                event.getStructure(CALO);

                IndexMap cherMap = loadMapByIndex(CHER, "pindex");
                IndexMap caloMap = loadMapByIndex(CALO, "pindex");

                for (int i = 0; i < PART.getRows(); i++) {
                    ParticleData pd = getParticleData(i, PART, CHER, CALO, cherMap, caloMap);

                    bool isPid11 = (pd.pid == ELECTRON_PID);
                    bool isSimple = isSimpleElectron(pd);

                    // PID 11 electrons
                    if (isPid11) {
                        total_electron_number_pid11++;
                        h1_nphe_pid11->Fill(pd.nphe);
                        h1_nphe_pid11_1->Fill(pd.nphe);
                        h1_charge_pid11->Fill(pd.charge);
                        h1_status_pid11->Fill(PART.getShort("status", i));
                        h1_vz_pid11->Fill(PART.getFloat("vz", i));
                        h1_chi2pid_pid11->Fill(PART.getFloat("chi2pid", i));
                        h1_energy_pcal_pid11->Fill(pd.energy_pcal);
                        h1_energy_pcal_pid11_1->Fill(pd.energy_pcal);
                        h1_sampling_fraction_pid11->Fill(pd.energy_total / pd.p);
                    }

                    // Simple electrons
                    if (isSimple) {
                        total_electron_number_simple++;
                        h1_nphe_tele->Fill(pd.nphe);
                        h1_nphe_tele_1->Fill(pd.nphe);
                        h1_charge_tele->Fill(pd.charge);
                        h1_status_tele->Fill(PART.getShort("status", i));
                        h1_vz_tele->Fill(PART.getFloat("vz", i));
                        h1_chi2pid_tele->Fill(PART.getFloat("chi2pid", i));
                        h1_energy_pcal_tele->Fill(pd.energy_pcal);
                        h1_energy_pcal_tele_1->Fill(pd.energy_pcal);
                        h1_sampling_fraction_tele->Fill(pd.energy_total / pd.p);
                    }

                    // Check for mismatches
                    if (isPid11 != isSimple) {
                        mismatchCount++;
                        if (mismatchCount <= 10) { // Limit output to first 10 for brevity
                            cout << "Mismatch at event " << event_count << ", particle " << i 
                                 << ": PID 11 = " << isPid11 << ", Simple = " << isSimple << endl;
                        }
                    }
                }
            }
        }
    }

    // Save histograms to ROOT file
    TFile* outFile = new TFile("output.root", "RECREATE");
    h1_nphe_tele->Write(); h1_nphe_pid11->Write();
    h1_nphe_tele_1->Write(); h1_nphe_pid11_1->Write();
    h1_charge_tele->Write(); h1_charge_pid11->Write();
    h1_status_tele->Write(); h1_status_pid11->Write();
    h1_vz_tele->Write(); h1_vz_pid11->Write();
    h1_chi2pid_tele->Write(); h1_chi2pid_pid11->Write();
    h1_energy_pcal_tele->Write(); h1_energy_pcal_pid11->Write();
    h1_energy_pcal_tele_1->Write(); h1_energy_pcal_pid11_1->Write();
    h1_sampling_fraction_tele->Write(); h1_sampling_fraction_pid11->Write();
    outFile->Close();

    // Save histograms as high-resolution PDFs with both sets overlaid
    gStyle->SetOptStat(0); // Disable stats box
    TCanvas* c1 = new TCanvas("c1", "Histograms", 800, 600);
    c1->SetHighLightColor(2);

    // Arrays for iteration
    TH1F* histSimple[] = {
        h1_nphe_tele, h1_nphe_tele_1, h1_charge_tele, h1_status_tele,
        h1_vz_tele, h1_chi2pid_tele, h1_energy_pcal_tele, h1_energy_pcal_tele_1,
        h1_sampling_fraction_tele
    };
    TH1F* histPid11[] = {
        h1_nphe_pid11, h1_nphe_pid11_1, h1_charge_pid11, h1_status_pid11,
        h1_vz_pid11, h1_chi2pid_pid11, h1_energy_pcal_pid11, h1_energy_pcal_pid11_1,
        h1_sampling_fraction_pid11
    };
    const char* pdfNames[] = {
        "nphe_comparison.pdf", "nphe_1_comparison.pdf", "charge_comparison.pdf", "status_comparison.pdf",
        "vz_comparison.pdf", "chi2pid_comparison.pdf", "energy_pcal_comparison.pdf", "energy_pcal_1_comparison.pdf",
        "sampling_fraction_comparison.pdf"
    };

    for (int i = 0; i < 9; i++) {
        c1->Clear();
        c1->SetGrid();

        // Draw PID 11 first with transparency, then simple electrons
        histPid11[i]->SetLineStyle(1); // Solid line
        histPid11[i]->SetFillStyle(0); // No fill
        histPid11[i]->Draw("HIST"); // Draw blue first
        histSimple[i]->SetLineStyle(2); // Dashed line for red
        histSimple[i]->SetFillStyle(0); // No fill
        histSimple[i]->Draw("HIST SAME"); // Overlay red

        // Add legend
        TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg->AddEntry(histPid11[i], "PID 11 Electrons", "l");
        leg->AddEntry(histSimple[i], "Simple Electrons", "l");
        leg->Draw();

        // Optional: Log scale for better visibility of differences
        // c1->SetLogy(); // Uncomment if you want log scale

        c1->Print(pdfNames[i], "pdf");
        delete leg;
    }

    // Print summary
    cout << "Total events processed: " << total_events_processed << endl;
    cout << "Total simple electrons found: " << total_electron_number_simple << endl;
    cout << "Total PID 11 electrons found: " << total_electron_number_pid11 << endl;
    cout << "Number of mismatches (PID 11 != Simple): " << mismatchCount << endl;
    cout << "Time elapsed: " << timer.RealTime() << " seconds" << endl;

    // Clean up
    delete h1_nphe_tele; delete h1_nphe_pid11;
    delete h1_nphe_tele_1; delete h1_nphe_pid11_1;
    delete h1_charge_tele; delete h1_charge_pid11;
    delete h1_status_tele; delete h1_status_pid11;
    delete h1_vz_tele; delete h1_vz_pid11;
    delete h1_chi2pid_tele; delete h1_chi2pid_pid11;
    delete h1_energy_pcal_tele; delete h1_energy_pcal_pid11;
    delete h1_energy_pcal_tele_1; delete h1_energy_pcal_pid11_1;
    delete h1_sampling_fraction_tele; delete h1_sampling_fraction_pid11;
    delete outFile;
    delete c1;

    return 0;
}