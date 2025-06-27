/* #include <cstdlib>
#include <iostream>
#include <filesystem>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <TFile.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "reader.h"

namespace fs = std::filesystem;

// Structure to hold particle data
struct ParticleData {
    float px, py, pz, p, vz;
    short status;
    int pid;
    float chi2pid;
};

// Structure to hold binning parameters
struct BinConfig {
    int nBins;
    double min, max;
    std::string axisLabel;
};

// Function to create bins
std::vector<double> createBins(double min, double max, int nBins) {
    std::vector<double> bins;
    double width = (max - min) / nBins;
    for (int i = 0; i <= nBins; ++i) {
        bins.push_back(min + i * width);
    }
    return bins;
}

// Function to calculate p_t^2
double calculatePt2(TLorentzVector P_h, TLorentzVector q) {
    TVector3 q3 = q.Vect();
    TVector3 p_h3 = P_h.Vect();
    TVector3 cross_qph = q3.Cross(p_h3);
    double magnitudeCrossProductSquared = cross_qph.Mag2();
    double magnitudeQSquared = q3.Mag2();
    return magnitudeCrossProductSquared / magnitudeQSquared;
}

// Function to extract particle data
ParticleData getParticleData(int partIdx, hipo::bank& PART) {
    ParticleData pd;
    pd.px = PART.getFloat("px", partIdx);
    pd.py = PART.getFloat("py", partIdx);
    pd.pz = PART.getFloat("pz", partIdx);
    pd.p = sqrt(pd.px * pd.px + pd.py * pd.py + pd.pz * pd.pz);
    pd.vz = PART.getFloat("vz", partIdx);
    pd.status = PART.getShort("status", partIdx);
    pd.pid = PART.getInt("pid", partIdx);
    pd.chi2pid = PART.getFloat("chi2pid", partIdx);
    return pd;
}

// Function to load HIPO files
std::vector<std::string> loadHipoFiles(const std::string& dir) {
    std::vector<std::string> files;
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (entry.path().extension() == ".hipo") {
            files.push_back(entry.path().string());
        }
    }
    return files;
}

int main(int argc, char* argv[]) {
    // Check command-line arguments
    if (argc != 7) {
        std::cerr << "Usage: ./main <variable> <target1> <target2> <target3> <config> <pion_type>\n";
        std::cerr << "Example: ./main Q2 CxC LD2 CuSn OB PP\n";
        return 1;
    }

    std::string variable = argv[1]; // nu, Q2, z, pt2
    std::string target1 = argv[2];  // CxC
    std::string target2 = argv[3];  // LD2
    std::string target3 = argv[4];  // CuSn
    std::string config = argv[5];   // OB or IB
    std::string pion_type = argv[6]; // PP or NP

    // Validate inputs
    std::vector<std::string> validVars = {"nu", "Q2", "z", "pt2"};
    if (std::find(validVars.begin(), validVars.end(), variable) == validVars.end()) {
        std::cerr << "Invalid variable. Use: nu, Q2, z, pt2\n";
        return 1;
    }
    if (config != "OB" && config != "IB") {
        std::cerr << "Invalid config. Use: OB or IB\n";
        return 1;
    }
    if (pion_type != "PP" && pion_type != "NP") {
        std::cerr << "Invalid pion type. Use: PP or NP\n";
        return 1;
    }

    // Beam and target setup
    double beamEnergy = 10.532; // GeV
    auto db = TDatabasePDG::Instance();
    double protonMass = db->GetParticle(2212)->Mass();
    TLorentzVector beam(0, 0, beamEnergy, beamEnergy);
    TLorentzVector target(0, 0, 0, protonMass);

    // Binning configuration
    std::map<std::string, BinConfig> binConfigs = {
        {"nu", {10, 0.1, 1.0, "#nu (GeV)"}},
        {"Q2", {10, 1.0, 9.0, "Q^{2} (GeV^{2})"}},
        {"z", {10, 0.3, 0.7, "z"}},
        {"pt2", {10, 0.0, 1.5, "p_{t}^{2} (GeV^{2})"}}
    };
    BinConfig bc = binConfigs[variable];

    // Histograms for electrons and pions
    auto* h_LD2_ele = new TH1F("h_LD2_ele", ("LD2 Electrons;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_LD2_pi = new TH1F("h_LD2_pi", ("LD2 Pions;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_Cu_ele = new TH1F("h_Cu_ele", ("Cu Electrons;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_Cu_pi = new TH1F("h_Cu_pi", ("Cu Pions;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_Sn_ele = new TH1F("h_Sn_ele", ("Sn Electrons;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_Sn_pi = new TH1F("h_Sn_pi", ("Sn Pions;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_CxC_ele = new TH1F("h_CxC_ele", ("CxC Electrons;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_CxC_pi = new TH1F("h_CxC_pi", ("CxC Pions;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);

    // Total electron counts for z and pt2
    double total_LD2_ele = 0, total_Cu_ele = 0, total_Sn_ele = 0, total_CxC_ele = 0;

    // Base directory
    std::string baseDir = "/lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v11";
    std::map<std::string, std::string> targetDirs = {
        {"CxC", baseDir + "/CxC"},
        {"LD2", baseDir + "/LD2"},
        {"CuSn", baseDir + "/CuSn"}
    };

    // Vz cuts for outbending (OB)
    std::map<std::string, std::pair<double, double>> vzCutsEleOB = {
        {"Cu", {-10.5023, -6.4529}},
        {"Sn", {-6.0086, 5.0000}},
        {"LD2", {-20.0000, 5.0000}},
        {"CxC", {-10.5319, 5.0000}}
    };
    std::map<std::string, std::pair<double, double>> vzCutsPiPOB = {
        {"Cu", {-9.76024, -5.1676}},
        {"Sn", {-4.6413, 5.0000}},
        {"LD2", {-20.0000, 5.0000}},
        {"CxC", {-10.5319, 5.0000}}
    };
    std::map<std::string, std::pair<double, double>> vzCutsPiMOB = {
        {"Cu", {-10.9553, -5.8338}},
        {"Sn", {-5.4629, 5.0000}},
        {"LD2", {-20.0000, 5.0000}},
        {"CxC", {-10.3501, 5.0000}}
    };

    // Placeholder vz cuts for inbending (IB)
    std::map<std::string, std::pair<double, double>> vzCutsEleIB = {
        {"Cu", {0, 0}}, {"Sn", {0, 0}}, {"LD2", {0, 0}}, {"CxC", {0, 0}}
    };
    std::map<std::string, std::pair<double, double>> vzCutsPiPIB = {
        {"Cu", {0, 0}}, {"Sn", {0, 0}}, {"LD2", {0, 0}}, {"CxC", {0, 0}}
    };
    std::map<std::string, std::pair<double, double>> vzCutsPiMIB = {
        {"Cu", {0, 0}}, {"Sn", {0, 0}}, {"LD2", {0, 0}}, {"CxC", {0, 0}}
    };

    // Select vz cuts
    auto vzCutsEle = (config == "OB") ? vzCutsEleOB : vzCutsEleIB;
    auto vzCutsPi = (pion_type == "PP") ? ((config == "OB") ? vzCutsPiPOB : vzCutsPiPIB)
                                        : ((config == "OB") ? vzCutsPiMOB : vzCutsPiMIB);
    int pion_pid = (pion_type == "PP") ? 211 : -211;

    // Process each target
    std::vector<std::string> targets = {target1, target2, target3};
    for (const auto& tgt : targets) {
        if (targetDirs.find(tgt) == targetDirs.end()) {
            std::cerr << "Invalid target: " << tgt << "\n";
            continue;
        }

        // Load HIPO files
        std::vector<std::string> hipoFiles = loadHipoFiles(targetDirs[tgt]);
        if (hipoFiles.empty()) {
            std::cerr << "No HIPO files found in " << targetDirs[tgt] << "\n";
            continue;
        }

        // Counters for tracking
        long totalElectrons = 0, triggerElectrons = 0, kinematicElectrons = 0;
        long totalPions = 0, triggerPions = 0, kinematicPions = 0;
        long totalElectronsCu = 0, triggerElectronsCu = 0, kinematicElectronsCu = 0;
        long totalPionsCu = 0, triggerPionsCu = 0, kinematicPionsCu = 0;
        long totalElectronsSn = 0, triggerElectronsSn = 0, kinematicElectronsSn = 0;
        long totalPionsSn = 0, triggerPionsSn = 0, kinematicPionsSn = 0;

        std::cout << "Processing target: " << tgt << "\n";

        for (const auto& file : hipoFiles) {
            hipo::reader reader;
            reader.open(file.c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);
            hipo::bank RUN(factory.getSchema("RUN::config"));
            hipo::bank PART(factory.getSchema("REC::Particle"));
            hipo::event event;

            // Check torus from the first event
            if (reader.next()) {
                reader.read(event);
                event.getStructure(RUN);
                float torus = RUN.getFloat("torus", 0);
                bool isInbending = (torus < 0);
                if ((config == "OB" && isInbending) || (config == "IB" && !isInbending)) {
                    std::cout << "Skipping file (wrong config): " << file << "\n";
                    continue; // Skip to next file
                }
                // Reset reader to process all events, including the first
                reader.open(file.c_str());
            } else {
                std::cout << "No events in file: " << file << "\n";
                continue;
            }
            // Initialize event count          
                
            int event_count = 0;

            // Read events
            while (reader.next()) {

                event_count++;
                reader.read(event);
                event.getStructure(PART);
            
                // Extract particles
                std::vector<ParticleData> particles(PART.getRows());
                for (int i = 0; i < PART.getRows(); ++i) {
                    particles[i] = getParticleData(i, PART);
                }

                bool hasTrigger = false;
                TLorentzVector el_temp;
                double nu, Q2, y, W;
                std::string targetType = tgt; // Default for non-CuSn targets

                // Electron selection
                for (int i = 0; i < particles.size(); ++i) {
                    const auto& p = particles[i];
                    if (p.pid != 11) continue;
                    totalElectrons++;
                    if (tgt == "CuSn") {
                        // Check vz range to classify as Cu or Sn before any common cuts
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);
                
                        if (isCu) totalElectronsCu++;
                        if (isSn)  totalElectronsSn++;
                    }

                    // Common electron cuts
                    if (p.status >= 0 || i != 0 || p.chi2pid < -5 || p.chi2pid > 5) continue;
                    if (abs(p.status) / 1000 != 2) continue;

                    // For CuSn, check vz against both Cu and Sn ranges
                    if (tgt == "CuSn") {
                        triggerElectrons++;
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);

                        if (!isCu && !isSn) continue; // Skip if not in either range

                        // Process Cu
                        if (isCu) {
                            
                            
                            triggerElectronsCu++;
                            el_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(11)->Mass());

                            // Kinematic variables
                            nu = beamEnergy - el_temp.E();
                            TLorentzVector q = beam - el_temp;
                            Q2 = -q.Mag2();
                            y = nu / beamEnergy;
                            W = sqrt(protonMass * protonMass + 2 * protonMass * nu - Q2);

                            // Kinematic cuts
                            if (Q2 <= 1 || y <= 0.25 || y >= 0.85 || W <= 2) continue;
                            kinematicElectronsCu++;
                            hasTrigger = true;

                            // Fill electron histogram or counter for Cu
                            if (variable == "nu") h_Cu_ele->Fill(nu);
                            else if (variable == "Q2") h_Cu_ele->Fill(Q2);
                            else total_Cu_ele++;
                        }

                        // Process Sn
                        if (isSn) {
                           
                            
                            triggerElectronsSn++;
                            el_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(11)->Mass());

                            // Kinematic variables
                            nu = beamEnergy - el_temp.E();
                            TLorentzVector q = beam - el_temp;
                            Q2 = -q.Mag2();
                            y = nu / beamEnergy;
                            W = sqrt(protonMass * protonMass + 2 * protonMass * nu - Q2);

                            // Kinematic cuts
                            if (Q2 <= 1 || y <= 0.25 || y >= 0.85 || W <= 2) continue;
                            kinematicElectronsSn++;
                            hasTrigger = true;

                            // Fill electron histogram or counter for Sn
                            if (variable == "nu") h_Sn_ele->Fill(nu);
                            else if (variable == "Q2") h_Sn_ele->Fill(Q2);
                            else total_Sn_ele++;
                        }
                    } else {
                        // Non-CuSn targets (LD2, CxC)
                        if (p.vz < vzCutsEle[tgt].first || p.vz > vzCutsEle[tgt].second) continue;

                        triggerElectrons++;
                        el_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(11)->Mass());

                        // Kinematic variables
                        nu = beamEnergy - el_temp.E();
                        TLorentzVector q = beam - el_temp;
                        Q2 = -q.Mag2();
                        y = nu / beamEnergy;
                        W = sqrt(protonMass * protonMass + 2 * protonMass * nu - Q2);

                        // Kinematic cuts
                        if (Q2 <= 1 || y <= 0.25 || y >= 0.85 || W <= 2) continue;

                        kinematicElectrons++;
                        hasTrigger = true;

                        // Fill electron histogram or counter
                        if (variable == "nu") {
                            if (tgt == "LD2") h_LD2_ele->Fill(nu);
                            else if (tgt == "CxC") h_CxC_ele->Fill(nu);
                        } else if (variable == "Q2") {
                            if (tgt == "LD2") h_LD2_ele->Fill(Q2);
                            else if (tgt == "CxC") h_CxC_ele->Fill(Q2);
                        } else {
                            if (tgt == "LD2") total_LD2_ele++;
                            else if (tgt == "CxC") total_CxC_ele++;
                        }
                    }
                    break; // One trigger electron per event
                }

                if (!hasTrigger) continue;

                // Pion selection
                for (const auto& p : particles) {
                    if (p.pid != pion_pid) continue;
                    totalPions++;
                    if (tgt == "CuSn") {
                        // Check vz range to classify as Cu or Sn before any common cuts
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);
                
                        if (isCu) totalPionsCu++;
                        if (isSn)  totalPionsSn++;
                    }

                    if (abs(p.status) / 1000 != 2) continue;
                    if (p.chi2pid < -5 || p.chi2pid > 5) continue;

                    // For CuSn, check pion vz against both Cu and Sn ranges
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsPi["Cu"].first && p.vz <= vzCutsPi["Cu"].second);
                        bool isSn = (p.vz >= vzCutsPi["Sn"].first && p.vz <= vzCutsPi["Sn"].second);

                        if (!isCu && !isSn) continue;

                        TLorentzVector pi_temp;
                        pi_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(pion_pid)->Mass());
                        double z = pi_temp.E() / nu;
                        double pt2 = calculatePt2(pi_temp, beam - el_temp);
                        if (z <= 0.3 || z >= 0.7 || pt2 >= 1.2) continue;

                        // Process Cu pions
                        if (isCu) {
                            kinematicPionsCu++;
                            if (variable == "nu") h_Cu_pi->Fill(nu);
                            else if (variable == "Q2") h_Cu_pi->Fill(Q2);
                            else if (variable == "z") h_Cu_pi->Fill(z);
                            else if (variable == "pt2") h_Cu_pi->Fill(pt2);
                        }

                        // Process Sn pions
                        if (isSn) {
                            kinematicPionsSn++;
                            if (variable == "nu") h_Sn_pi->Fill(nu);
                            else if (variable == "Q2") h_Sn_pi->Fill(Q2);
                            else if (variable == "z") h_Sn_pi->Fill(z);
                            else if (variable == "pt2") h_Sn_pi->Fill(pt2);
                        }
                    } else {
                        // Non-CuSn targets
                        if (p.vz < vzCutsPi[tgt].first || p.vz > vzCutsPi[tgt].second) continue;
                        TLorentzVector pi_temp;
                        pi_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(pion_pid)->Mass());
                        double z = pi_temp.E() / nu;
                        double pt2 = calculatePt2(pi_temp, beam - el_temp);
                        if (z <= 0.3 || z >= 0.7 || pt2 >= 1.2) continue;

                        kinematicPions++;
                        if (variable == "nu") {
                            if (tgt == "LD2") h_LD2_pi->Fill(nu);
                            else if (tgt == "CxC") h_CxC_pi->Fill(nu);
                        } else if (variable == "Q2") {
                            if (tgt == "LD2") h_LD2_pi->Fill(Q2);
                            else if (tgt == "CxC") h_CxC_pi->Fill(Q2);
                        } else if (variable == "z") {
                            if (tgt == "LD2") h_LD2_pi->Fill(z);
                            else if (tgt == "CxC") h_CxC_pi->Fill(z);
                        } else if (variable == "pt2") {
                            if (tgt == "LD2") h_LD2_pi->Fill(pt2);
                            else if (tgt == "CxC") h_CxC_pi->Fill(pt2);
                        }
                    }
                }
            }
            if (event_count > 10000) {
                std::cout << "Processed " << event_count << " events in file: " << file << "\n";
                break; // Limit to first 10000 events for testing
            }
        }

        // Print counts
        std::cout << tgt << " Electron counts: Total=" << totalElectrons
                  << ", Trigger=" << triggerElectrons
                  << ", Kinematic=" << kinematicElectrons << "\n";
        std::cout << tgt << " Pion counts: Total=" << totalPions
                  << ", Kinematic=" << kinematicPions << "\n";
        if (tgt == "CuSn") {
            std::cout << "Cu Electron counts: Total=" << totalElectronsCu
                      << ", Kinematic=" << kinematicElectronsCu << "\n";
            std::cout << "Cu Pion counts: Total=" << totalPionsCu
                      << ", Kinematic=" << kinematicPionsCu << "\n";
            std::cout << "Sn Electron counts: Total=" << totalElectronsSn
                      << ", Kinematic=" << kinematicElectronsSn << "\n";
            std::cout << "Sn Pion counts: Total=" << totalPionsSn
                      << ", Kinematic=" << kinematicPionsSn << "\n";
        }
    }

    // Calculate multiplicity ratio R
    TGraphErrors* g_Cu = new TGraphErrors();
    TGraphErrors* g_Sn = new TGraphErrors();
    TGraphErrors* g_CxC = new TGraphErrors();
    int pointCu = 0, pointSn = 0, pointCxC = 0;
    double maxY = -1e30, minY = 1e30;

    bool useTotalElectrons = (variable == "z" || variable == "pt2");
    for (int bin = 1; bin <= bc.nBins; ++bin) {
        double binCenter = h_LD2_pi->GetBinCenter(bin);
        double ld2_pi = h_LD2_pi->GetBinContent(bin);
        double cu_pi = h_Cu_pi->GetBinContent(bin);
        double sn_pi = h_Sn_pi->GetBinContent(bin);
        double cxc_pi = h_CxC_pi->GetBinContent(bin);

        double ld2_ele = useTotalElectrons ? total_LD2_ele : h_LD2_ele->GetBinContent(bin);
        double cu_ele = useTotalElectrons ? total_Cu_ele : h_Cu_ele->GetBinContent(bin);
        double sn_ele = useTotalElectrons ? total_Sn_ele : h_Sn_ele->GetBinContent(bin);
        double cxc_ele = useTotalElectrons ? total_CxC_ele : h_CxC_ele->GetBinContent(bin);

        // Cu ratio
        if (ld2_pi > 0 && cu_pi > 0 && ld2_ele > 0 && cu_ele > 0) {
            double R = (cu_pi / cu_ele) / (ld2_pi / ld2_ele);
            double err = R * sqrt(1.0 / cu_pi + 1.0 / cu_ele + 1.0 / ld2_pi + 1.0 / ld2_ele);
            g_Cu->SetPoint(pointCu, binCenter, R);
            g_Cu->SetPointError(pointCu, 0, err);
            pointCu++;
            maxY = std::max(maxY, R + err);
            minY = std::min(minY, R - err);
        }

        // Sn ratio
        if (ld2_pi > 0 && sn_pi > 0 && ld2_ele > 0 && sn_ele > 0) {
            double R = (sn_pi / sn_ele) / (ld2_pi / ld2_ele);
            double err = R * sqrt(1.0 / sn_pi + 1.0 / sn_ele + 1.0 / ld2_pi + 1.0 / ld2_ele);
            g_Sn->SetPoint(pointSn, binCenter, R);
            g_Sn->SetPointError(pointSn, 0, err);
            pointSn++;
            maxY = std::max(maxY, R + err);
            minY = std::min(minY, R - err);
        }

        // CxC ratio
        if (ld2_pi > 0 && cxc_pi > 0 && ld2_ele > 0 && cxc_ele > 0) {
            double R = (cxc_pi / cxc_ele) / (ld2_pi / ld2_ele);
            double err = R * sqrt(1.0 / cxc_pi + 1.0 / cxc_ele + 1.0 / ld2_pi + 1.0 / ld2_ele);
            g_CxC->SetPoint(pointCxC, binCenter, R);
            g_CxC->SetPointError(pointCxC, 0, err);
            pointCxC++;
            maxY = std::max(maxY, R + err);
            minY = std::min(minY, R - err);
        }
    }

    // Create canvas and plot
    TCanvas* canvas = new TCanvas("canvas", ("Multiplicity Ratio vs " + variable).c_str(), 800, 600);
    g_Cu->SetTitle(("Multiplicity Ratio vs " + bc.axisLabel).c_str());
    g_Cu->GetXaxis()->SetTitle(bc.axisLabel.c_str());
    g_Cu->GetYaxis()->SetTitle(pion_type == "PP" ? "R_{M}^{#pi^{+}}" : "R_{M}^{#pi^{-}}");
    double margin = 0.1 * (maxY - minY);
    g_Cu->SetMinimum(minY - margin);
    g_Cu->SetMaximum(maxY + margin);
    g_Cu->SetMarkerStyle(20);
    g_Cu->SetMarkerColor(kRed);
    g_Cu->SetLineColor(kRed);
    g_Cu->Draw("AP");

    g_Sn->SetMarkerStyle(21);
    g_Sn->SetMarkerColor(kBlue);
    g_Sn->SetLineColor(kBlue);
    g_Sn->Draw("P SAME");

    g_CxC->SetMarkerStyle(22);
    g_CxC->SetMarkerColor(kGreen);
    g_CxC->SetLineColor(kGreen);
    g_CxC->Draw("P SAME");

    TLegend* legend = new TLegend(0.11, 0.9, 0.3, 0.8);
    legend->AddEntry(g_CxC, "CxC", "p");
    legend->AddEntry(g_Cu, "Cu", "p");
    legend->AddEntry(g_Sn, "Sn", "p");
    legend->SetTextSize(0.03);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();

    std::string output = "MultiplicityRatio_" + config + "_" + pion_type + "_" + variable + ".png";
    canvas->SaveAs(output.c_str());

    return 0;
} */

/* #include <cstdlib>
#include <iostream>
#include <filesystem>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "reader.h"

namespace fs = std::filesystem;

// Structure to hold particle data
struct ParticleData {
    float px, py, pz, p, vz;
    short status;
    int pid;
    float chi2pid;
};

// Structure to hold binning parameters
struct BinConfig {
    int nBins;
    double min, max;
    std::string axisLabel;
};

// Function to create bins
std::vector<double> createBins(double min, double max, int nBins) {
    std::vector<double> bins;
    double width = (max - min) / nBins;
    for (int i = 0; i <= nBins; ++i) {
        bins.push_back(min + i * width);
    }
    return bins;
}

// Function to calculate p_t^2
double calculatePt2(TLorentzVector P_h, TLorentzVector q) {
    TVector3 q3 = q.Vect();
    TVector3 p_h3 = P_h.Vect();
    TVector3 cross_qph = q3.Cross(p_h3);
    double magnitudeCrossProductSquared = cross_qph.Mag2();
    double magnitudeQSquared = q3.Mag2();
    return magnitudeCrossProductSquared / magnitudeQSquared;
}

// Function to extract particle data
ParticleData getParticleData(int partIdx, hipo::bank& PART) {
    ParticleData pd;
    pd.px = PART.getFloat("px", partIdx);
    pd.py = PART.getFloat("py", partIdx);
    pd.pz = PART.getFloat("pz", partIdx);
    pd.p = sqrt(pd.px * pd.px + pd.py * pd.py + pd.pz * pd.pz);
    pd.vz = PART.getFloat("vz", partIdx);
    pd.status = PART.getShort("status", partIdx);
    pd.pid = PART.getInt("pid", partIdx);
    pd.chi2pid = PART.getFloat("chi2pid", partIdx);
    return pd;
}

// Function to load HIPO files
std::vector<std::string> loadHipoFiles(const std::string& dir) {
    std::vector<std::string> files;
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (entry.path().extension() == ".hipo") {
            files.push_back(entry.path().string());
        }
    }
    return files;
}

int main(int argc, char* argv[]) {
    // Check command-line arguments
    if (argc != 7) {
        std::cerr << "Usage: ./main <variable> <target1> <target2> <target3> <config> <pion_type>\n";
        std::cerr << "Example: ./main Q2 CxC LD2 CuSn OB PP\n";
        return 1;
    }

    std::string variable = argv[1]; // nu, Q2, z, pt2
    std::string target1 = argv[2];  // CxC
    std::string target2 = argv[3];  // LD2
    std::string target3 = argv[4];  // CuSn
    std::string config = argv[5];   // OB or IB
    std::string pion_type = argv[6]; // PP or NP

    // Validate inputs
    std::vector<std::string> validVars = {"nu", "Q2", "z", "pt2"};
    if (std::find(validVars.begin(), validVars.end(), variable) == validVars.end()) {
        std::cerr << "Invalid variable. Use: nu, Q2, z, pt2\n";
        return 1;
    }
    if (config != "OB" && config != "IB") {
        std::cerr << "Invalid config. Use: OB or IB\n";
        return 1;
    }
    if (pion_type != "PP" && pion_type != "NP") {
        std::cerr << "Invalid pion type. Use: PP or NP\n";
        return 1;
    }

    // Open log file
    std::string txtOutput = "ProcessedFiles_" + config + "_" + pion_type + "_" + variable + ".txt";
    std::ofstream logFile(txtOutput);
    if (!logFile.is_open()) {
        std::cerr << "Failed to open " << txtOutput << " for writing\n";
    }

    // Beam and target setup
    double beamEnergy = 10.532; // GeV
    auto db = TDatabasePDG::Instance();
    double protonMass = db->GetParticle(2212)->Mass();
    TLorentzVector beam(0, 0, beamEnergy, beamEnergy);
    TLorentzVector target(0, 0, 0, protonMass);

    // Binning configuration
    std::map<std::string, BinConfig> binConfigs = {
        {"nu", {10, 0.1, 1.0, "#nu (GeV)"}},
        {"Q2", {10, 1.0, 9.0, "Q^{2} (GeV^{2})"}},
        {"z", {10, 0.2, 1.0, "z"}},
        {"pt2", {10, 0.0, 1.5, "p_{t}^{2} (GeV^{2})"}}
    };
    BinConfig bc = binConfigs[variable];

    // Histograms for electrons and pions
    auto* h_LD2_ele = new TH1F("h_LD2_ele", ("LD2 Electrons;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_LD2_pi = new TH1F("h_LD2_pi", ("LD2 Pions;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_Cu_ele = new TH1F("h_Cu_ele", ("Cu Electrons;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_Cu_pi = new TH1F("h_Cu_pi", ("Cu Pions;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_Sn_ele = new TH1F("h_Sn_ele", ("Sn Electrons;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_Sn_pi = new TH1F("h_Sn_pi", ("Sn Pions;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_CxC_ele = new TH1F("h_CxC_ele", ("CxC Electrons;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_CxC_pi = new TH1F("h_CxC_pi", ("CxC Pions;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);

    // Total electron counts for z and pt2
    double total_LD2_ele = 0, total_Cu_ele = 0, total_Sn_ele = 0, total_CxC_ele = 0;

    // Base directory
    std::string baseDir = "/lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v11";
    std::map<std::string, std::string> targetDirs = {
        {"CxC", baseDir + "/CxC"},
        {"LD2", baseDir + "/LD2"},
        {"CuSn", baseDir + "/CuSn"}
    };

    // Vz cuts for outbending (OB)
    std::map<std::string, std::pair<double, double>> vzCutsEleOB = {
        {"Cu", {-10.5023, -6.4529}},
        {"Sn", {-6.0086, 5.0000}},
        {"LD2", {-20.0000, 5.0000}},
        {"CxC", {-10.5319, 5.0000}}
    };
    std::map<std::string, std::pair<double, double>> vzCutsPiPOB = {
        {"Cu", {-9.76024, -5.1676}},
        {"Sn", {-4.6413, 5.0000}},
        {"LD2", {-20.0000, 5.0000}},
        {"CxC", {-10.5319, 5.0000}}
    };
    std::map<std::string, std::pair<double, double>> vzCutsPiMOB = {
        {"Cu", {-10.9553, -5.8338}},
        {"Sn", {-5.4629, 5.0000}},
        {"LD2", {-20.0000, 5.0000}},
        {"CxC", {-10.3501, 5.0000}}
    };

    // Placeholder vz cuts for inbending (IB)
    std::map<std::string, std::pair<double, double>> vzCutsEleIB = {
        {"Cu", {-9.22, -5.46}}, {"Sn", {-4.101, 5}}, {"LD2", {-15, 5}}, {"CxC", {-9.128, 5}}
    };
    std::map<std::string, std::pair<double, double>> vzCutsPiPIB = {
        {"Cu", {-10.69, -6.50}}, {"Sn", {-6.13, 5}}, {"LD2", {-15, 5}}, {"CxC", {-10.28, 5}}
    };
    std::map<std::string, std::pair<double, double>> vzCutsPiMIB = {
        {"Cu", {-9.74, -5.52}}, {"Sn", {-4.742, 5}}, {"LD2", {-15, 5}}, {"CxC", {-10.35, 5}}
    };

    // Select vz cuts
    auto vzCutsEle = (config == "OB") ? vzCutsEleOB : vzCutsEleIB;
    auto vzCutsPi = (pion_type == "PP") ? ((config == "OB") ? vzCutsPiPOB : vzCutsPiPIB)
                                        : ((config == "OB") ? vzCutsPiMOB : vzCutsPiMIB);
    int pion_pid = (pion_type == "PP") ? 211 : -211;

    const long maxEvents = 100000; // Total event limit per target

    // Process each target
    std::vector<std::string> inputTargets = {target1, target2, target3};
    for (const auto& tgt : inputTargets) {
        if (targetDirs.find(tgt) == targetDirs.end()) {
            std::cerr << "Invalid target: " << tgt << "\n";
            if (logFile.is_open()) {
                logFile << "Target: " << tgt << "\nInvalid target\n\n";
            }
            continue;
        }

        // Load HIPO files
        std::vector<std::string> hipoFiles = loadHipoFiles(targetDirs[tgt]);
        if (hipoFiles.empty()) {
            std::cout << "No .hipo files found for target " << tgt << "\n";
            if (logFile.is_open()) {
                logFile << "Target: " << tgt << "\nNo .hipo files found\n\n";
            }
            continue;
        }

        // Log the start of processing for this target
        std::cout << "Processing .hipo files for target " << tgt << ":\n";
        if (logFile.is_open()) {
            logFile << "Target: " << tgt << "\n";
        }

        // Counters for tracking
        long totalEvents = 0, triggerEvents = 0, kinematicEvents = 0;
        long totalPions = 0, kinematicPions = 0;
        long totalEventsCu = 0, triggerEventsCu = 0, kinematicEventsCu = 0;
        long totalPionsCu = 0, kinematicPionsCu = 0;
        long totalEventsSn = 0, triggerEventsSn = 0, kinematicEventsSn = 0;
        long totalPionsSn = 0, kinematicPionsSn = 0;
        long event_count = 0;
        long event_countCu = 0, event_countSn = 0;

        for (const auto& file : hipoFiles) {
            // Check if we've reached the event limit for this target
            if (tgt == "CuSn") {
                if (event_countCu >= maxEvents && event_countSn >= maxEvents) {
                    std::cout << "Reached " << maxEvents << " kinematic events for Cu and Sn in " << tgt << "\n";
                    break;
                }
            } else {
                if (event_count >= maxEvents) {
                    std::cout << "Reached " << maxEvents << " kinematic events for " << tgt << "\n";
                    break;
                }
            }

            // Log the file being processed
            std::cout << "  - " << file << "\n";
            if (logFile.is_open()) {
                logFile << "  - " << file << "\n";
            }

            hipo::reader reader;
            reader.open(file.c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);
            hipo::bank RUN(factory.getSchema("RUN::config"));
            hipo::bank PART(factory.getSchema("REC::Particle"));
            hipo::event event;
            bool isFirstEvent = true;

            while (reader.next()) {
                if (tgt == "CuSn") {
                    if (event_countCu >= maxEvents && event_countSn >= maxEvents) break;
                } else {
                    if (event_count >= maxEvents) break; // Stop if total event limit reached
                }

                reader.read(event);
                event.getStructure(PART);

                // Check torus from the first event
                if (isFirstEvent) {
                    event.getStructure(RUN);
                    float torus = RUN.getFloat("torus", 0);
                    bool isInbending = (torus < 0);
                    if ((config == "OB" && isInbending) || (config == "IB" && !isInbending)) {
                        std::cout << "Skipping file (wrong config): " << file << "\n";
                        if (logFile.is_open()) {
                            logFile << "  Skipped (wrong config)\n";
                        }
                        break; // Skip to next file
                    }
                    isFirstEvent = false;
                }

                // Extract particles
                std::vector<ParticleData> particles(PART.getRows());
                for (int i = 0; i < PART.getRows(); ++i) {
                    particles[i] = getParticleData(i, PART);
                }

                bool hasTrigger = false;
                TLorentzVector el_temp;
                double nu, Q2, y, W;
                std::string targetType = tgt; // Default for non-CuSn targets

                // Electron selection
                for (int i = 0; i < particles.size(); ++i) {
                    const auto& p = particles[i];
                    if (p.pid != 11) continue;
                    totalEvents++;
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);
                        if (isCu) totalEventsCu++;
                        if (isSn) totalEventsSn++;
                    }

                    // Common electron cuts
                    if (p.status >= 0 || i != 0 || p.chi2pid < -5 || p.chi2pid > 5) continue;
                    if (abs(p.status) / 1000 != 2) continue;

                    // For CuSn, check vz against both Cu and Sn ranges
                    if (tgt == "CuSn") {
                        triggerEvents++;
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);

                        if (!isCu && !isSn) continue;

                        // Process Cu
                        if (isCu && event_countCu < maxEvents) {
                            triggerEventsCu++;
                            el_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(11)->Mass());

                            // Kinematic variables
                            nu = beamEnergy - el_temp.E();
                            TLorentzVector q = beam - el_temp;
                            Q2 = -q.Mag2();
                            y = nu / beamEnergy;
                            W = sqrt(protonMass * protonMass + 2 * protonMass * nu - Q2);

                            // Kinematic cuts
                            if (Q2 <= 1 || y <= 0.25 || y >= 0.85 || W <= 2) continue;
                            event_countCu++;
                            kinematicEventsCu++;
                            hasTrigger = true;

                            // Fill electron histogram or counter for Cu
                            if (variable == "nu") h_Cu_ele->Fill(nu);
                            else if (variable == "Q2") h_Cu_ele->Fill(Q2);
                            else total_Cu_ele++;
                        }

                        // Process Sn
                        if (isSn && event_countSn < maxEvents) {
                            triggerEventsSn++;
                            el_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(11)->Mass());

                            // Kinematic variables
                            nu = beamEnergy - el_temp.E();
                            TLorentzVector q = beam - el_temp;
                            Q2 = -q.Mag2();
                            y = nu / beamEnergy;
                            W = sqrt(protonMass * protonMass + 2 * protonMass * nu - Q2);

                            // Kinematic cuts
                            if (Q2 <= 1 || y <= 0.25 || y >= 0.85 || W <= 2) continue;
                            event_countSn++;
                            kinematicEventsSn++;
                            hasTrigger = true;

                            // Fill electron histogram or counter for Sn
                            if (variable == "nu") h_Sn_ele->Fill(nu);
                            else if (variable == "Q2") h_Sn_ele->Fill(Q2);
                            else total_Sn_ele++;
                        }
                    } else {
                        // Non-CuSn targets (LD2, CxC)
                        if (p.vz < vzCutsEle[tgt].first || p.vz > vzCutsEle[tgt].second) continue;

                        triggerEvents++;
                        el_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(11)->Mass());

                        // Kinematic variables
                        nu = beamEnergy - el_temp.E();
                        TLorentzVector q = beam - el_temp;
                        Q2 = -q.Mag2();
                        y = nu / beamEnergy;
                        W = sqrt(protonMass * protonMass + 2 * protonMass * nu - Q2);

                        // Kinematic cuts
                        if (Q2 <= 1 || y <= 0.25 || y >= 0.85 || W <= 2) continue;

                        event_count++;
                        kinematicEvents++;
                        hasTrigger = true;

                        // Fill electron histogram or counter
                        if (variable == "nu") {
                            if (tgt == "LD2") h_LD2_ele->Fill(nu);
                            else if (tgt == "CxC") h_CxC_ele->Fill(nu);
                        } else if (variable == "Q2") {
                            if (tgt == "LD2") h_LD2_ele->Fill(Q2);
                            else if (tgt == "CxC") h_CxC_ele->Fill(Q2);
                        } else {
                            if (tgt == "LD2") total_LD2_ele++;
                            else if (tgt == "CxC") total_CxC_ele++;
                        }
                    }
                    break; // One trigger electron per event
                }

                if (!hasTrigger) continue;

                // Pion selection
                for (const auto& p : particles) {
                    if (p.pid != pion_pid) continue;
                    totalPions++;
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);
                        if (isCu) totalPionsCu++;
                        if (isSn) totalPionsSn++;
                    }

                    if (abs(p.status) / 1000 != 2) continue;
                    if (p.chi2pid < -5 || p.chi2pid > 5) continue;

                    // For CuSn, check pion vz against both Cu and Sn ranges
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsPi["Cu"].first && p.vz <= vzCutsPi["Cu"].second);
                        bool isSn = (p.vz >= vzCutsPi["Sn"].first && p.vz <= vzCutsPi["Sn"].second);

                        if (!isCu && !isSn) continue;

                        TLorentzVector pi_temp;
                        pi_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(pion_pid)->Mass());
                        double z = pi_temp.E() / nu;
                        double pt2 = calculatePt2(pi_temp, beam - el_temp);
                        if (z <= 0.2 || z >= 1.0 || pt2 >= 1.2) continue;

                        // Process Cu pions
                        if (isCu) {
                            kinematicPionsCu++;
                            if (variable == "nu") h_Cu_pi->Fill(nu);
                            else if (variable == "Q2") h_Cu_pi->Fill(Q2);
                            else if (variable == "z") h_Cu_pi->Fill(z);
                            else if (variable == "pt2") h_Cu_pi->Fill(pt2);
                        }

                        // Process Sn pions
                        if (isSn) {
                            kinematicPionsSn++;
                            if (variable == "nu") h_Sn_pi->Fill(nu);
                            else if (variable == "Q2") h_Sn_pi->Fill(Q2);
                            else if (variable == "z") h_Sn_pi->Fill(z);
                            else if (variable == "pt2") h_Sn_pi->Fill(pt2);
                        }
                    } else {
                        // Non-CuSn targets
                        if (p.vz < vzCutsPi[tgt].first || p.vz > vzCutsPi[tgt].second) continue;
                        TLorentzVector pi_temp;
                        pi_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(pion_pid)->Mass());
                        double z = pi_temp.E() / nu;
                        double pt2 = calculatePt2(pi_temp, beam - el_temp);
                        if (z <= 0.2 || z >= 1.0 || pt2 >= 1.2) continue;

                        kinematicPions++;
                        if (variable == "nu") {
                            if (tgt == "LD2") h_LD2_pi->Fill(nu);
                            else if (tgt == "CxC") h_CxC_pi->Fill(nu);
                        } else if (variable == "Q2") {
                            if (tgt == "LD2") h_LD2_pi->Fill(Q2);
                            else if (tgt == "CxC") h_CxC_pi->Fill(Q2);
                        } else if (variable == "z") {
                            if (tgt == "LD2") h_LD2_pi->Fill(z);
                            else if (tgt == "CxC") h_CxC_pi->Fill(z);
                        } else if (variable == "pt2") {
                            if (tgt == "LD2") h_LD2_pi->Fill(pt2);
                            else if (tgt == "CxC") h_CxC_pi->Fill(pt2);
                        }
                    }
                }
            }
        }

        // Print counts
        std::cout << tgt << " Event counts: Total=" << totalEvents
                  << ", Trigger=" << triggerEvents
                  << ", Kinematic=" << kinematicEvents << "\n";
        std::cout << tgt << " Pion counts: Total=" << totalPions
                  << ", Kinematic=" << kinematicPions << "\n";
        if (tgt == "CuSn") {
            std::cout << "Cu Event counts: Total=" << totalEventsCu
                      << ", Trigger=" << triggerEventsCu
                      << ", Kinematic=" << kinematicEventsCu << "\n";
            std::cout << "Cu Pion counts: Total=" << totalPionsCu
                      << ", Kinematic=" << kinematicPionsCu << "\n";
            std::cout << "Sn Event counts: Total=" << totalEventsSn
                      << ", Trigger=" << triggerEventsSn
                      << ", Kinematic=" << kinematicEventsSn << "\n";
            std::cout << "Sn Pion counts: Total=" << totalPionsSn
                      << ", Kinematic=" << kinematicPionsSn << "\n";
        }
    }

    // Close log file
    if (logFile.is_open()) {
        logFile.close();
        std::cout << "Saved processed files to " << txtOutput << "\n";
    }

    // Calculate multiplicity ratio R
    TGraphErrors* g_Cu = new TGraphErrors();
    TGraphErrors* g_Sn = new TGraphErrors();
    TGraphErrors* g_CxC = new TGraphErrors();
    int pointCu = 0, pointSn = 0, pointCxC = 0;
    double maxY = -1e30, minY = 1e30;

    bool useTotalElectrons = (variable == "z" || variable == "pt2");
    for (int bin = 1; bin <= bc.nBins; ++bin) {
        double binCenter = h_LD2_pi->GetBinCenter(bin);
        double ld2_pi = h_LD2_pi->GetBinContent(bin);
        double cu_pi = h_Cu_pi->GetBinContent(bin);
        double sn_pi = h_Sn_pi->GetBinContent(bin);
        double cxc_pi = h_CxC_pi->GetBinContent(bin);

        double ld2_ele = useTotalElectrons ? total_LD2_ele : h_LD2_ele->GetBinContent(bin);
        double cu_ele = useTotalElectrons ? total_Cu_ele : h_Cu_ele->GetBinContent(bin);
        double sn_ele = useTotalElectrons ? total_Sn_ele : h_Sn_ele->GetBinContent(bin);
        double cxc_ele = useTotalElectrons ? total_CxC_ele : h_CxC_ele->GetBinContent(bin);

        // Cu ratio
        if (ld2_pi > 0 && cu_pi > 0 && ld2_ele > 0 && cu_ele > 0) {
            double R = (cu_pi / cu_ele) / (ld2_pi / ld2_ele);
            double err = R * sqrt(1.0 / cu_pi + 1.0 / cu_ele + 1.0 / ld2_pi + 1.0 / ld2_ele);
            g_Cu->SetPoint(pointCu, binCenter, R);
            g_Cu->SetPointError(pointCu, 0, err);
            pointCu++;
            maxY = std::max(maxY, R + err);
            minY = std::min(minY, R - err);
        }

        // Sn ratio
        if (ld2_pi > 0 && sn_pi > 0 && ld2_ele > 0 && sn_ele > 0) {
            double R = (sn_pi / sn_ele) / (ld2_pi / ld2_ele);
            double err = R * sqrt(1.0 / sn_pi + 1.0 / sn_ele + 1.0 / ld2_pi + 1.0 / ld2_ele);
            g_Sn->SetPoint(pointSn, binCenter, R);
            g_Sn->SetPointError(pointSn, 0, err);
            pointSn++;
            maxY = std::max(maxY, R + err);
            minY = std::min(minY, R - err);
        }

        // CxC ratio
        if (ld2_pi > 0 && cxc_pi > 0 && ld2_ele > 0 && cxc_ele > 0) {
            double R = (cxc_pi / cxc_ele) / (ld2_pi / ld2_ele);
            double err = R * sqrt(1.0 / cxc_pi + 1.0 / cxc_ele + 1.0 / ld2_pi + 1.0 / ld2_ele);
            g_CxC->SetPoint(pointCxC, binCenter, R);
            g_CxC->SetPointError(pointCxC, 0, err);
            pointCxC++;
            maxY = std::max(maxY, R + err);
            minY = std::min(minY, R - err);
        }
    }

    // Create canvas and plot
    TCanvas* canvas = new TCanvas("canvas", ("Multiplicity Ratio vs " + variable).c_str(),1000, 800);
    canvas->SetLeftMargin(0.15); // Increased left margin
    canvas->SetBottomMargin(0.15); // Increased bottom margin
    g_Cu->SetTitle(("Multiplicity Ratio vs " + bc.axisLabel).c_str());
    g_Cu->GetXaxis()->SetTitle(bc.axisLabel.c_str());
    g_Cu->GetYaxis()->SetTitle(pion_type == "PP" ? "R_{M}^{#pi^{+}}" : "R_{M}^{#pi^{-}}");
    double margin = 0.1 * (maxY - minY);
    g_Cu->SetMinimum(minY - margin);
    g_Cu->SetMaximum(1.3);
    g_Cu->SetMarkerStyle(20);
    g_Cu->SetMarkerColor(kRed);
    g_Cu->SetLineColor(kRed);
    g_Cu->Draw("AP");

    g_Sn->SetMarkerStyle(21);
    g_Sn->SetMarkerColor(kBlue);
    g_Sn->SetLineColor(kBlue);
    g_Sn->Draw("P SAME");

    g_CxC->SetMarkerStyle(22);
    g_CxC->SetMarkerColor(kGreen);
    g_CxC->SetLineColor(kGreen);
    g_CxC->Draw("P SAME");

 
TLegend* legend = new TLegend(0.8, 0.78, 0.95, 0.88); 
    legend->AddEntry(g_CxC, "CxC", "p");
    legend->AddEntry(g_Cu, "Cu", "p");
    legend->AddEntry(g_Sn, "Sn", "p");
    legend->SetTextSize(0.03);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();

    std::string output = "MultiplicityRatio_" + config + "_" + pion_type + "_" + variable + ".pdf";
    canvas->SaveAs(output.c_str());

    return 0;
} */
#include <cstdlib>
#include <iostream>
#include <filesystem>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "reader.h"

namespace fs = std::filesystem;

// Structure to hold particle data
struct ParticleData {
    float px, py, pz, p, vz;
    short status;
    int pid;
    float chi2pid;
};

// Structure to hold binning parameters
struct BinConfig {
    int nBins;
    double min, max;
    std::string axisLabel;
};

// Function to create bins
std::vector<double> createBins(double min, double max, int nBins) {
    std::vector<double> bins;
    double width = (max - min) / nBins;
    for (int i = 0; i <= nBins; ++i) {
        bins.push_back(min + i * width);
    }
    return bins;
}

// Function to calculate p_t^2
double calculatePt2(TLorentzVector P_h, TLorentzVector q) {
    TVector3 q3 = q.Vect();
    TVector3 p_h3 = P_h.Vect();
    TVector3 cross_qph = q3.Cross(p_h3);
    double magnitudeCrossProductSquared = cross_qph.Mag2();
    double magnitudeQSquared = q3.Mag2();
    return magnitudeCrossProductSquared / magnitudeQSquared;
}

// Function to extract particle data
ParticleData getParticleData(int partIdx, hipo::bank& PART) {
    ParticleData pd;
    pd.px = PART.getFloat("px", partIdx);
    pd.py = PART.getFloat("py", partIdx);
    pd.pz = PART.getFloat("pz", partIdx);
    pd.p = sqrt(pd.px * pd.px + pd.py * pd.py + pd.pz * pd.pz);
    pd.vz = PART.getFloat("vz", partIdx);
    pd.status = PART.getShort("status", partIdx);
    pd.pid = PART.getInt("pid", partIdx);
    pd.chi2pid = PART.getFloat("chi2pid", partIdx);
    return pd;
}

// Function to load HIPO files
std::vector<std::string> loadHipoFiles(const std::string& dir) {
    std::vector<std::string> files;
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (entry.path().extension() == ".hipo") {
            files.push_back(entry.path().string());
        }
    }
    return files;
}

int main(int argc, char* argv[]) {
    // Check command-line arguments
    if (argc != 7) {
        std::cerr << "Usage: ./main <variable> <target1> <target2> <target3> <config> <pion_type>\n";
        std::cerr << "Example: ./main Q2 CxC LD2 CuSn OB PP\n";
        return 1;
    }

    std::string variable = argv[1]; // nu, Q2, z, pt2
    std::string target1 = argv[2];  // CxC
    std::string target2 = argv[3];  // LD2
    std::string target3 = argv[4];  // CuSn
    std::string config = argv[5];   // OB or IB
    std::string pion_type = argv[6]; // PP or NP

    // Validate inputs
    std::vector<std::string> validVars = {"nu", "Q2", "z", "pt2"};
    if (std::find(validVars.begin(), validVars.end(), variable) == validVars.end()) {
        std::cerr << "Invalid variable. Use: nu, Q2, z, pt2\n";
        return 1;
    }
    if (config != "OB" && config != "IB") {
        std::cerr << "Invalid config. Use: OB or IB\n";
        return 1;
    }
    if (pion_type != "PP" && pion_type != "NP") {
        std::cerr << "Invalid pion type. Use: PP or NP\n";
        return 1;
    }

    // Open log file
    std::string txtOutput = "ProcessedFiles_" + config + "_" + pion_type + "_" + variable + ".txt";
    std::ofstream logFile(txtOutput);
    if (!logFile.is_open()) {
        std::cerr << "Failed to open " << txtOutput << " for writing\n";
    }

    // Beam and target setup
    double beamEnergy = 10.532; // GeV
    auto db = TDatabasePDG::Instance();
    double protonMass = db->GetParticle(2212)->Mass();
    TLorentzVector beam(0, 0, beamEnergy, beamEnergy);
    TLorentzVector target(0, 0, 0, protonMass);

    // Binning configuration
    std::map<std::string, BinConfig> binConfigs = {
        {"nu", {10, 0.1, 1.0, "#nu (GeV)"}},
        {"Q2", {10, 1.0, 9.0, "Q^{2} (GeV^{2})"}},
        {"z", {10, 0.2, 1.0, "z"}},
        {"pt2", {10, 0.0, 1.5, "p_{t}^{2} (GeV^{2})"}}
    };
    BinConfig bc = binConfigs[variable];

    // Histograms for electrons and pions
    auto* h_LD2_ele = new TH1F("h_LD2_ele", ("LD2 Electrons;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_LD2_pi = new TH1F("h_LD2_pi", ("LD2 Pions;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_Cu_ele = new TH1F("h_Cu_ele", ("Cu Electrons;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_Cu_pi = new TH1F("h_Cu_pi", ("Cu Pions;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_Sn_ele = new TH1F("h_Sn_ele", ("Sn Electrons;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_Sn_pi = new TH1F("h_Sn_pi", ("Sn Pions;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_CxC_ele = new TH1F("h_CxC_ele", ("CxC Electrons;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);
    auto* h_CxC_pi = new TH1F("h_CxC_pi", ("CxC Pions;" + bc.axisLabel + ";Counts").c_str(), bc.nBins, bc.min, bc.max);

    // Total electron counts for z and pt2
    double total_LD2_ele = 0, total_Cu_ele = 0, total_Sn_ele = 0, total_CxC_ele = 0;

    // Base directory
    std::string baseDir = "/lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v11";
    std::map<std::string, std::string> targetDirs = {
        {"CxC", baseDir + "/CxC"},
        {"LD2", baseDir + "/LD2"},
        {"CuSn", baseDir + "/CuSn"}
    };

    // Vz cuts for outbending (OB)
    std::map<std::string, std::pair<double, double>> vzCutsEleOB = {
        {"Cu", {-10.5023, -6.4529}},
        {"Sn", {-6.0086, 5.0000}},
        {"LD2", {-20.0000, 5.0000}},
        {"CxC", {-10.5319, 5.0000}}
    };
    std::map<std::string, std::pair<double, double>> vzCutsPiPOB = {
        {"Cu", {-9.76024, -5.1676}},
        {"Sn", {-4.6413, 5.0000}},
        {"LD2", {-20.0000, 5.0000}},
        {"CxC", {-10.5319, 5.0000}}
    };
    std::map<std::string, std::pair<double, double>> vzCutsPiMOB = {
        {"Cu", {-10.9553, -5.8338}},
        {"Sn", {-5.4629, 5.0000}},
        {"LD2", {-20.0000, 5.0000}},
        {"CxC", {-10.3501, 5.0000}}
    };

    // Vz cuts for inbending (IB)
    std::map<std::string, std::pair<double, double>> vzCutsEleIB = {
        {"Cu", {-9.22, -5.46}}, {"Sn", {-4.101, 5}}, {"LD2", {-15, 5}}, {"CxC", {-9.128, 5}}
    };
    std::map<std::string, std::pair<double, double>> vzCutsPiPIB = {
        {"Cu", {-10.69, -6.50}}, {"Sn", {-6.13, 5}}, {"LD2", {-15, 5}}, {"CxC", {-10.28, 5}}
    };
    std::map<std::string, std::pair<double, double>> vzCutsPiMIB = {
        {"Cu", {-9.74, -5.52}}, {"Sn", {-4.742, 5}}, {"LD2", {-15, 5}}, {"CxC", {-10.35, 5}}
    };

    // Select vz cuts
    auto vzCutsEle = (config == "OB") ? vzCutsEleOB : vzCutsEleIB;
    auto vzCutsPi = (pion_type == "PP") ? ((config == "OB") ? vzCutsPiPOB : vzCutsPiPIB)
                                        : ((config == "OB") ? vzCutsPiMOB : vzCutsPiMIB);
    int pion_pid = (pion_type == "PP") ? 211 : -211;

    // Process each target
    std::vector<std::string> inputTargets = {target1, target2, target3};
    for (const auto& tgt : inputTargets) {
        if (targetDirs.find(tgt) == targetDirs.end()) {
            std::cerr << "Invalid target: " << tgt << "\n";
            if (logFile.is_open()) {
                logFile << "Target: " << tgt << "\nInvalid target\n\n";
            }
            continue;
        }

        // Load HIPO files
        std::vector<std::string> hipoFiles = loadHipoFiles(targetDirs[tgt]);
        if (hipoFiles.empty()) {
            std::cout << "No .hipo files found for target " << tgt << "\n";
            if (logFile.is_open()) {
                logFile << "Target: " << tgt << "\nNo .hipo files found\n\n";
            }
            continue;
        }

        // Log the start of processing for this target
        std::cout << "Processing .hipo files for target " << tgt << ":\n";
        if (logFile.is_open()) {
            logFile << "Target: " << tgt << "\n";
        }

        // Counters for tracking
        long totalElectrons = 0, triggerElectrons = 0, kinematicElectrons = 0;
        long totalPions = 0, kinematicPions = 0;
        long totalElectronsCu = 0, triggerElectronsCu = 0, kinematicElectronsCu = 0;
        long totalPionsCu = 0, kinematicPionsCu = 0;
        long totalElectronsSn = 0, triggerElectronsSn = 0, kinematicElectronsSn = 0;
        long totalPionsSn = 0, kinematicPionsSn = 0;
        long event_count = 0;

        // Set max events based on target
        long maxEvents = (tgt == "CuSn") ? 200000000 : 100000000;

        for (const auto& file : hipoFiles) {
            if (event_count >= maxEvents) {
                std::cout << "Reached " << maxEvents << " events for " << tgt << "\n";
                if (logFile.is_open()) {
                    logFile << "Reached " << maxEvents << " events\n\n";
                }
                break;
            }

            // Log the file being processed
            std::cout << "  - " << file << "\n";
            if (logFile.is_open()) {
                logFile << "  - " << file << "\n";
            }

            hipo::reader reader;
            reader.open(file.c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);
            hipo::bank RUN(factory.getSchema("RUN::config"));
            hipo::bank PART(factory.getSchema("REC::Particle"));
            hipo::event event;
            bool isFirstEvent = true;

            while (reader.next()) {
                if (event_count >= maxEvents) break;

                event_count++;
                reader.read(event);
                event.getStructure(PART);

                // Check torus from the first event
                if (isFirstEvent) {
                    event.getStructure(RUN);
                    float torus = RUN.getFloat("torus", 0);
                    bool isInbending = (torus < 0);
                    if ((config == "OB" && isInbending) || (config == "IB" && !isInbending)) {
                        std::cout << "Skipping file (wrong config): " << file << "\n";
                        if (logFile.is_open()) {
                            logFile << "  Skipped (wrong config)\n";
                        }
                        break; // Skip to next file
                    }
                    isFirstEvent = false;
                }

                // Extract particles
                std::vector<ParticleData> particles(PART.getRows());
                for (int i = 0; i < PART.getRows(); ++i) {
                    particles[i] = getParticleData(i, PART);
                }

                bool hasTrigger = false;
                TLorentzVector el_temp;
                double nu, Q2, y, W;

                // Electron selection
                for (int i = 0; i < particles.size(); ++i) {
                    const auto& p = particles[i];
                    if (p.pid != 11) continue;
                    totalElectrons++;
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);
                        if (isCu) totalElectronsCu++;
                        if (isSn) totalElectronsSn++;
                    }

                    // Common electron cuts
                    if (p.status >= 0 || i != 0 || p.chi2pid < -5 || p.chi2pid > 5) continue;
                    if (abs(p.status) / 1000 != 2) continue;

                    // For CuSn, check vz against both Cu and Sn ranges
                    if (tgt == "CuSn") {
                        triggerElectrons++;
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);

                        if (!isCu && !isSn) continue;

                        // Process Cu
                        if (isCu) {
                            triggerElectronsCu++;
                            el_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(11)->Mass());

                            // Kinematic variables
                            nu = beamEnergy - el_temp.E();
                            TLorentzVector q = beam - el_temp;
                            Q2 = -q.Mag2();
                            y = nu / beamEnergy;
                            W = sqrt(protonMass * protonMass + 2 * protonMass * nu - Q2);

                            // Kinematic cuts
                            if (Q2 <= 1 || y <= 0.25 || y >= 0.85 || W <= 2) continue;
                            kinematicElectronsCu++;
                            hasTrigger = true;

                            // Fill electron histogram or counter for Cu
                            if (variable == "nu") h_Cu_ele->Fill(nu);
                            else if (variable == "Q2") h_Cu_ele->Fill(Q2);
                            else total_Cu_ele++;
                        }

                        // Process Sn
                        if (isSn) {
                            triggerElectronsSn++;
                            el_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(11)->Mass());

                            // Kinematic variables
                            nu = beamEnergy - el_temp.E();
                            TLorentzVector q = beam - el_temp;
                            Q2 = -q.Mag2();
                            y = nu / beamEnergy;
                            W = sqrt(protonMass * protonMass + 2 * protonMass * nu - Q2);

                            // Kinematic cuts
                            if (Q2 <= 1 || y <= 0.25 || y >= 0.85 || W <= 2) continue;
                            kinematicElectronsSn++;
                            hasTrigger = true;

                            // Fill electron histogram or counter for Sn
                            if (variable == "nu") h_Sn_ele->Fill(nu);
                            else if (variable == "Q2") h_Sn_ele->Fill(Q2);
                            else total_Sn_ele++;
                        }
                    } else {
                        // Non-CuSn targets (LD2, CxC)
                        if (p.vz < vzCutsEle[tgt].first || p.vz > vzCutsEle[tgt].second) continue;

                        triggerElectrons++;
                        el_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(11)->Mass());

                        // Kinematic variables
                        nu = beamEnergy - el_temp.E();
                        TLorentzVector q = beam - el_temp;
                        Q2 = -q.Mag2();
                        y = nu / beamEnergy;
                        W = sqrt(protonMass * protonMass + 2 * protonMass * nu - Q2);

                        // Kinematic cuts
                        if (Q2 <= 1 || y <= 0.25 || y >= 0.85 || W <= 2) continue;

                        kinematicElectrons++;
                        hasTrigger = true;

                        // Fill electron histogram or counter
                        if (variable == "nu") {
                            if (tgt == "LD2") h_LD2_ele->Fill(nu);
                            else if (tgt == "CxC") h_CxC_ele->Fill(nu);
                        } else if (variable == "Q2") {
                            if (tgt == "LD2") h_LD2_ele->Fill(Q2);
                            else if (tgt == "CxC") h_CxC_ele->Fill(Q2);
                        } else {
                            if (tgt == "LD2") total_LD2_ele++;
                            else if (tgt == "CxC") total_CxC_ele++;
                        }
                    }
                    break; // One trigger electron per event
                }

                if (!hasTrigger) continue;

                // Pion selection
                for (const auto& p : particles) {
                    if (p.pid != pion_pid) continue;
                    totalPions++;
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);
                        if (isCu) totalPionsCu++;
                        if (isSn) totalPionsSn++;
                    }

                    if (abs(p.status) / 1000 != 2) continue;
                    if (p.chi2pid < -5 || p.chi2pid > 5) continue;

                    // For CuSn, check pion vz against both Cu and Sn ranges
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsPi["Cu"].first && p.vz <= vzCutsPi["Cu"].second);
                        bool isSn = (p.vz >= vzCutsPi["Sn"].first && p.vz <= vzCutsPi["Sn"].second);

                        if (!isCu && !isSn) continue;

                        TLorentzVector pi_temp;
                        pi_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(pion_pid)->Mass());
                        double z = pi_temp.E() / nu;
                        double pt2 = calculatePt2(pi_temp, beam - el_temp);
                        if (z <= 0.2 || z >= 1.0 || pt2 >= 1.2) continue;

                        // Process Cu pions
                        if (isCu) {
                            kinematicPionsCu++;
                            if (variable == "nu") h_Cu_pi->Fill(nu);
                            else if (variable == "Q2") h_Cu_pi->Fill(Q2);
                            else if (variable == "z") h_Cu_pi->Fill(z);
                            else if (variable == "pt2") h_Cu_pi->Fill(pt2);
                        }

                        // Process Sn pions
                        if (isSn) {
                            kinematicPionsSn++;
                            if (variable == "nu") h_Sn_pi->Fill(nu);
                            else if (variable == "Q2") h_Sn_pi->Fill(Q2);
                            else if (variable == "z") h_Sn_pi->Fill(z);
                            else if (variable == "pt2") h_Sn_pi->Fill(pt2);
                        }
                    } else {
                        // Non-CuSn targets
                        if (p.vz < vzCutsPi[tgt].first || p.vz > vzCutsPi[tgt].second) continue;
                        TLorentzVector pi_temp;
                        pi_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(pion_pid)->Mass());
                        double z = pi_temp.E() / nu;
                        double pt2 = calculatePt2(pi_temp, beam - el_temp);
                        if (z <= 0.2 || z >= 1.0 || pt2 >= 1.2) continue;

                        kinematicPions++;
                        if (variable == "nu") {
                            if (tgt == "LD2") h_LD2_pi->Fill(nu);
                            else if (tgt == "CxC") h_CxC_pi->Fill(nu);
                        } else if (variable == "Q2") {
                            if (tgt == "LD2") h_LD2_pi->Fill(Q2);
                            else if (tgt == "CxC") h_CxC_pi->Fill(Q2);
                        } else if (variable == "z") {
                            if (tgt == "LD2") h_LD2_pi->Fill(z);
                            else if (tgt == "CxC") h_CxC_pi->Fill(z);
                        } else if (variable == "pt2") {
                            if (tgt == "LD2") h_LD2_pi->Fill(pt2);
                            else if (tgt == "CxC") h_CxC_pi->Fill(pt2);
                        }
                    }
                }
            }
        }

        // Print counts
        std::cout << tgt << " Electron counts: Total=" << totalElectrons
                  << ", Trigger=" << triggerElectrons
                  << ", Kinematic=" << kinematicElectrons << "\n";
        std::cout << tgt << " Pion counts: Total=" << totalPions
                  << ", Kinematic=" << kinematicPions << "\n";
        if (tgt == "CuSn") {
            std::cout << "Cu Electron counts: Total=" << totalElectronsCu
                      << ", Trigger=" << triggerElectronsCu
                      << ", Kinematic=" << kinematicElectronsCu << "\n";
            std::cout << "Cu Pion counts: Total=" << totalPionsCu
                      << ", Kinematic=" << kinematicPionsCu << "\n";
            std::cout << "Sn Electron counts: Total=" << totalElectronsSn
                      << ", Trigger=" << triggerElectronsSn
                      << ", Kinematic=" << kinematicElectronsSn << "\n";
            std::cout << "Sn Pion counts: Total=" << totalPionsSn
                      << ", Kinematic=" << kinematicPionsSn << "\n";
        }
        if (logFile.is_open()) {
            logFile << "Electron counts: Total=" << totalElectrons
                    << ", Trigger=" << triggerElectrons
                    << ", Kinematic=" << kinematicElectrons << "\n";
            logFile << "Pion counts: Total=" << totalPions
                    << ", Kinematic=" << kinematicPions << "\n";
            if (tgt == "CuSn") {
                logFile << "Cu Electron counts: Total=" << totalElectronsCu
                        << ", Trigger=" << triggerElectronsCu
                        << ", Kinematic=" << kinematicElectronsCu << "\n";
                logFile << "Cu Pion counts: Total=" << totalPionsCu
                        << ", Kinematic=" << kinematicPionsCu << "\n";
                logFile << "Sn Electron counts: Total=" << totalElectronsSn
                        << ", Trigger=" << triggerElectronsSn
                        << ", Kinematic=" << kinematicElectronsSn << "\n";
                logFile << "Sn Pion counts: Total=" << totalPionsSn
                        << ", Kinematic=" << kinematicPionsSn << "\n";
            }
            logFile << "\n";
        }
    }

    // Close log file
    if (logFile.is_open()) {
        logFile.close();
        std::cout << "Saved processed files to " << txtOutput << "\n";
    }

    // Calculate multiplicity ratio R
    TGraphErrors* g_Cu = new TGraphErrors();
    TGraphErrors* g_Sn = new TGraphErrors();
    TGraphErrors* g_CxC = new TGraphErrors();
    int pointCu = 0, pointSn = 0, pointCxC = 0;
    double maxY = -1e30, minY = 1e30;

    bool useTotalElectrons = (variable == "z" || variable == "pt2");
    for (int bin = 1; bin <= bc.nBins; ++bin) {
        double binCenter = h_LD2_pi->GetBinCenter(bin);
        double ld2_pi = h_LD2_pi->GetBinContent(bin);
        double cu_pi = h_Cu_pi->GetBinContent(bin);
        double sn_pi = h_Sn_pi->GetBinContent(bin);
        double cxc_pi = h_CxC_pi->GetBinContent(bin);

        double ld2_ele = useTotalElectrons ? total_LD2_ele : h_LD2_ele->GetBinContent(bin);
        double cu_ele = useTotalElectrons ? total_Cu_ele : h_Cu_ele->GetBinContent(bin);
        double sn_ele = useTotalElectrons ? total_Sn_ele : h_Sn_ele->GetBinContent(bin);
        double cxc_ele = useTotalElectrons ? total_CxC_ele : h_CxC_ele->GetBinContent(bin);

        // Cu ratio
        if (ld2_pi > 0 && cu_pi > 0 && ld2_ele > 0 && cu_ele > 0) {
            double R = (cu_pi / cu_ele) / (ld2_pi / ld2_ele);
            double err = R * sqrt(1.0 / cu_pi + 1.0 / cu_ele + 1.0 / ld2_pi + 1.0 / ld2_ele);
            g_Cu->SetPoint(pointCu, binCenter, R);
            g_Cu->SetPointError(pointCu, 0, err);
            pointCu++;
            maxY = std::max(maxY, R + err);
            minY = std::min(minY, R - err);
        }

        // Sn ratio
        if (ld2_pi > 0 && sn_pi > 0 && ld2_ele > 0 && sn_ele > 0) {
            double R = (sn_pi / sn_ele) / (ld2_pi / ld2_ele);
            double err = R * sqrt(1.0 / sn_pi + 1.0 / sn_ele + 1.0 / ld2_pi + 1.0 / ld2_ele);
            g_Sn->SetPoint(pointSn, binCenter, R);
            g_Sn->SetPointError(pointSn, 0, err);
            pointSn++;
            maxY = std::max(maxY, R + err);
            minY = std::min(minY, R - err);
        }

        // CxC ratio
        if (ld2_pi > 0 && cxc_pi > 0 && ld2_ele > 0 && cxc_ele > 0) {
            double R = (cxc_pi / cxc_ele) / (ld2_pi / ld2_ele);
            double err = R * sqrt(1.0 / cxc_pi + 1.0 / cxc_ele + 1.0 / ld2_pi + 1.0 / ld2_ele);
            g_CxC->SetPoint(pointCxC, binCenter, R);
            g_CxC->SetPointError(pointCxC, 0, err);
            pointCxC++;
            maxY = std::max(maxY, R + err);
            minY = std::min(minY, R - err);
        }
    }

    // Create canvas and plot
    TCanvas* canvas = new TCanvas("canvas", ("Multiplicity Ratio vs " + variable).c_str(), 1000, 800);
    canvas->SetLeftMargin(0.15); // Increased left margin
    canvas->SetBottomMargin(0.15); // Increased bottom margin
    g_Cu->SetTitle(("Multiplicity Ratio vs " + bc.axisLabel).c_str());
    g_Cu->GetXaxis()->SetTitle(bc.axisLabel.c_str());
    g_Cu->GetYaxis()->SetTitle(pion_type == "PP" ? "R_{M}^{#pi^{+}}" : "R_{M}^{#pi^{-}}");
    double margin = 0.1 * (maxY - minY);
    g_Cu->SetMinimum(0.3);
    g_Cu->SetMaximum(1.1);
    g_Cu->SetMarkerStyle(20);
    g_Cu->SetMarkerColor(kRed);
    g_Cu->SetLineColor(kRed);
    g_Cu->Draw("AP");

    g_Sn->SetMarkerStyle(21);
    g_Sn->SetMarkerColor(kBlue);
    g_Sn->SetLineColor(kBlue);
    g_Sn->Draw("P SAME");

    g_CxC->SetMarkerStyle(22);
    g_CxC->SetMarkerColor(kGreen);
    g_CxC->SetLineColor(kGreen);
    g_CxC->Draw("P SAME");

    TLegend* legend = new TLegend(0.8, 0.78, 0.95, 0.88);
    legend->AddEntry(g_CxC, "CxC", "p");
    legend->AddEntry(g_Cu, "Cu", "p");
    legend->AddEntry(g_Sn, "Sn", "p");
    legend->SetTextSize(0.03);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();

    std::string base_output = "MultiplicityRatio_" + config + "_" + pion_type + "_" + variable;
canvas->SaveAs((base_output + ".pdf").c_str());
TFile* outFile = new TFile((base_output + ".root").c_str(), "RECREATE");
canvas->Write("MultiplicityCanvas");
outFile->Close();

    return 0;
}