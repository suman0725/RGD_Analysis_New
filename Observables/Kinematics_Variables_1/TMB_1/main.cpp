/* #include <cstdlib>
#include <iostream>
#include <filesystem>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <TAxis.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include "reader.h"

namespace fs = std::filesystem;

// Structure to hold particle data
struct ParticleData {
    float px, py, pz, p, vz;
    short status;
    int pid;
    float chi2pid;
};

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
    if (argc != 7) {
        std::cerr << "Usage: ./main <target1> <target2> <target3> <config> <pion_type> <kin_var>\n";
        std::cerr << "Example: ./main CxC LD2 CuSn OB PP z\n";
        std::cerr << "kin_var options: z, Q2, nu\n";
        return 1;
    }

    std::string target1 = argv[1];  // CxC
    std::string target2 = argv[2];  // LD2
    std::string target3 = argv[3];  // CuSn
    std::string config = argv[4];   // OB or IB
    std::string pion_type = argv[5]; // PP or NP
    std::string kin_var = argv[6];  // z, Q2, or nu

    if (config != "OB" && config != "IB") {
        std::cerr << "Invalid config. Use: OB or IB\n";
        return 1;
    }
    if (pion_type != "PP" && pion_type != "NP") {
        std::cerr << "Invalid pion type. Use: PP or NP\n";
        return 1;
    }
    if (kin_var != "z" && kin_var != "Q2" && kin_var != "nu") {
        std::cerr << "Invalid kinematic variable. Use: z, Q2, or nu\n";
        return 1;
    }

    double beamEnergy = 10.532; // GeV
    auto db = TDatabasePDG::Instance();
    double protonMass = db->GetParticle(2212)->Mass();
    TLorentzVector beam(0, 0, beamEnergy, beamEnergy);
    TLorentzVector target(0, 0, 0, protonMass);

    // Binning configuration for kinematic variable
    const int kinBins = 10;
    double kinMin, kinMax;
    std::string kinLabel;
    if (kin_var == "z") {
        kinMin = 0.3;
        kinMax = 0.7;
        kinLabel = "z";
    } else if (kin_var == "Q2") {
        kinMin = 1.0;
        kinMax = 9.0;
        kinLabel = "Q^{2} (GeV^{2})";
    } else { // nu
        kinMin = 0.0;
        kinMax = 10.0;
        kinLabel = "#nu (GeV)";
    }

    // Binning configuration for p_T^2
    const int pt2Bins = 10;
    const double pt2Min = 0.0;
    const double pt2Max = 1.2;

    // Store p_T^2 sums and counts per kinematic and p_T^2 bin
    std::map<std::string, std::vector<std::vector<double>>> pt2_sums; // target -> [kinBin][pt2Bin]
    std::map<std::string, std::vector<std::vector<double>>> pt2_counts; // target -> [kinBin][pt2Bin]
    std::vector<std::string> targets = {"LD2", "Cu", "Sn", "CxC"};
    for (const auto& tgt : targets) {
        pt2_sums[tgt] = std::vector<std::vector<double>>(kinBins, std::vector<double>(pt2Bins, 0.0));
        pt2_counts[tgt] = std::vector<std::vector<double>>(kinBins, std::vector<double>(pt2Bins, 0.0));
    }

    // Total electron counts
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
        {"Cu", {-10.5023, -6.4529}}, {"Sn", {-6.0086, 5.0000}}, {"LD2", {-20.0000, 5.0000}}, {"CxC", {-10.5319, 5.0000}}
    };
    std::map<std::string, std::pair<double, double>> vzCutsPiPOB = {
        {"Cu", {-9.76024, -5.1676}}, {"Sn", {-4.6413, 5.0000}}, {"LD2", {-20.0000, 5.0000}}, {"CxC", {-10.5319, 5.0000}}
    };
    std::map<std::string, std::pair<double, double>> vzCutsPiMOB = {
        {"Cu", {-10.9553, -5.8338}}, {"Sn", {-5.4629, 5.0000}}, {"LD2", {-20.0000, 5.0000}}, {"CxC", {-10.3501, 5.0000}}
    };

    // Placeholder vz cuts for inbending (IB)
    std::map<std::string, std::pair<double, double>> vzCutsEleIB = {{"Cu", {0, 0}}, {"Sn", {0, 0}}, {"LD2", {0, 0}}, {"CxC", {0, 0}}};
    std::map<std::string, std::pair<double, double>> vzCutsPiPIB = {{"Cu", {0, 0}}, {"Sn", {0, 0}}, {"LD2", {0, 0}}, {"CxC", {0, 0}}};
    std::map<std::string, std::pair<double, double>> vzCutsPiMIB = {{"Cu", {0, 0}}, {"Sn", {0, 0}}, {"LD2", {0, 0}}, {"CxC", {0, 0}}};

    auto vzCutsEle = (config == "OB") ? vzCutsEleOB : vzCutsEleIB;
    auto vzCutsPi = (pion_type == "PP") ? ((config == "OB") ? vzCutsPiPOB : vzCutsPiPIB)
                                       : ((config == "OB") ? vzCutsPiMOB : vzCutsPiMIB);
    int pion_pid = (pion_type == "PP") ? 211 : -211;

    std::vector<std::string> inputTargets = {target1, target2, target3};
    for (const auto& tgt : inputTargets) {
        if (targetDirs.find(tgt) == targetDirs.end()) continue;

        std::vector<std::string> hipoFiles = loadHipoFiles(targetDirs[tgt]);
        if (hipoFiles.empty()) continue;

        long totalElectrons = 0, triggerElectrons = 0, kinematicElectrons = 0;
        long totalPions = 0, triggerPions = 0, kinematicPions = 0;
        int event_count = 0;
        for (const auto& file : hipoFiles) {
            hipo::reader reader;
            reader.open(file.c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);
            hipo::bank RUN(factory.getSchema("RUN::config"));
            hipo::bank PART(factory.getSchema("REC::Particle"));
            hipo::event event;
            bool isFirstEvent = true;
            

            while (reader.next()) {
                event_count++;
                if (event_count > 1000) break;
                reader.read(event);
                event.getStructure(PART);
                event.getStructure(RUN);

                if (isFirstEvent) {
                    float torus = RUN.getFloat("torus", 0);
                    bool isInbending = (torus < 0);
                    if ((config == "OB" && isInbending) || (config == "IB" && !isInbending)) break;
                    isFirstEvent = false;
                }

                std::vector<ParticleData> particles(PART.getRows());
                for (int i = 0; i < PART.getRows(); ++i) {
                    particles[i] = getParticleData(i, PART);
                }

                bool hasTrigger = false;
                TLorentzVector el_temp;
                double nu = 0.0;
                TLorentzVector q;

                for (const auto& p : particles) {
                    if (p.pid != 11) continue;
                    totalElectrons++;
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);
                        if (isCu) totalElectrons++;
                        if (isSn) totalElectrons++;
                    }

                    if (p.status >= 0 || abs(p.status) / 1000 != 2 || p.chi2pid < -5 || p.chi2pid > 5) continue;
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);
                        if (!isCu && !isSn) continue;
                    } else if (p.vz < vzCutsEle[tgt].first || p.vz > vzCutsEle[tgt].second) continue;

                    triggerElectrons++;
                    el_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(11)->Mass());
                    nu = beamEnergy - el_temp.E();
                    q = beam - el_temp;
                    double Q2 = -q.Mag2();
                    double y = nu / beamEnergy;
                    double W = sqrt(protonMass * protonMass + 2 * protonMass * nu - Q2);
                    if (Q2 <= 1 || y <= 0.25 || y >= 0.85 || W <= 2) continue;

                    kinematicElectrons++;
                    hasTrigger = true;
                    if (tgt == "LD2") total_LD2_ele++;
                    else if (tgt == "CxC") total_CxC_ele++;
                    else if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);
                        if (isCu) total_Cu_ele++;
                        if (isSn) total_Sn_ele++;
                    }
                    break;
                }

                if (!hasTrigger) continue;

                for (const auto& p : particles) {
                    if (p.pid != pion_pid) continue;
                    totalPions++;
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsPi["Cu"].first && p.vz <= vzCutsPi["Cu"].second);
                        bool isSn = (p.vz >= vzCutsPi["Sn"].first && p.vz <= vzCutsPi["Sn"].second);
                        if (isCu || isSn) {}
                    }

                    if (abs(p.status) / 1000 != 2 || p.chi2pid < -10 || p.chi2pid > 10) continue;
                    TLorentzVector pi_temp;
                    pi_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(pion_pid)->Mass());
                    double z = pi_temp.E() / nu;
                    double pt2 = calculatePt2(pi_temp, q);
                    double Q2_val = -q.Mag2();

                    // Apply cuts
                    if (z <= 0.3 || z >= 0.7 || pt2 >= 1.2) continue;

                    triggerPions++;
                    double kin_value;
                    if (kin_var == "z") kin_value = z;
                    else if (kin_var == "Q2") kin_value = Q2_val;
                    else kin_value = nu;

                    // Apply kinematic variable cuts
                    if (kin_value < kinMin || kin_value > kinMax) continue;

                    int kinBin = static_cast<int>((kin_value - kinMin) / ((kinMax - kinMin) / kinBins));
                    if (kinBin < 0 || kinBin >= kinBins) continue;

                    int pt2Bin = static_cast<int>((pt2 - pt2Min) / ((pt2Max - pt2Min) / pt2Bins));
                    if (pt2Bin < 0 || pt2Bin >= pt2Bins) continue;

                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsPi["Cu"].first && p.vz <= vzCutsPi["Cu"].second);
                        bool isSn = (p.vz >= vzCutsPi["Sn"].first && p.vz <= vzCutsPi["Sn"].second);
                        if (!isCu && !isSn) continue;
                        if (isCu) {
                            pt2_sums["Cu"][kinBin][pt2Bin] += pt2;
                            pt2_counts["Cu"][kinBin][pt2Bin] += 1.0;
                        }
                        if (isSn) {
                            pt2_sums["Sn"][kinBin][pt2Bin] += pt2;
                            pt2_counts["Sn"][kinBin][pt2Bin] += 1.0;
                        }
                    } else {
                        if (p.vz < vzCutsPi[tgt].first || p.vz > vzCutsPi[tgt].second) continue;
                        pt2_sums[tgt][kinBin][pt2Bin] += pt2;
                        pt2_counts[tgt][kinBin][pt2Bin] += 1.0;
                    }
                    kinematicPions++;
                }
            }
            

        }


        std::cout << tgt << " Electron counts: Total=" << totalElectrons
                  << ", Trigger=" << triggerElectrons
                  << ", Kinematic=" << kinematicElectrons << "\n";
        std::cout << tgt << " Pion counts: Total=" << totalPions
                  << ", Trigger=" << triggerPions
                  << ", Kinematic=" << kinematicPions << "\n";
    }

    // Calculate <p_T^2> using weighted average over p_T^2 bins
    std::map<std::string, std::vector<double>> avg_pT2, err_pT2;
    for (const auto& tgt : targets) {
        avg_pT2[tgt].resize(kinBins, 0.0);
        err_pT2[tgt].resize(kinBins, 0.0);
        for (int kinBin = 0; kinBin < kinBins; ++kinBin) {
            double sum_weighted_pT2 = 0.0, sum_weights = 0.0, sum_pT2_sq = 0.0;
            for (int pt2Bin = 0; pt2Bin < pt2Bins; ++pt2Bin) {
                double N_j = pt2_counts[tgt][kinBin][pt2Bin];
                if (N_j <= 0) continue;
                double pt2_j = pt2_sums[tgt][kinBin][pt2Bin] / N_j; // Average p_T^2 in this sub-bin
                sum_weighted_pT2 += pt2_j * N_j; // Weighted contribution
                sum_weights += N_j; // Total weight (number of entries)
                sum_pT2_sq += N_j * pt2_j * pt2_j; // Weighted sum of (p_T^2)^2
            }
            if (sum_weights > 0) {
                avg_pT2[tgt][kinBin] = sum_weighted_pT2 / sum_weights;
                double avg_pT2_sq = sum_pT2_sq / sum_weights;
                double variance = avg_pT2_sq - avg_pT2[tgt][kinBin] * avg_pT2[tgt][kinBin];
                err_pT2[tgt][kinBin] = (variance >= 0) ? sqrt(variance / sum_weights) : 0.0;
            }

            double kinLow = kinMin + kinBin * (kinMax - kinMin) / kinBins;
            double kinHigh = kinMin + (kinBin + 1) * (kinMax - kinMin) / kinBins;
            std::cout << "Target " << tgt << ", " << kin_var << " Bin " << (kinBin + 1)
                      << " (" << kinLow << " to " << kinHigh << ")"
                      << ", <p_T^2> = " << avg_pT2[tgt][kinBin] << " GeV^2"
                      << ", Error = " << err_pT2[tgt][kinBin] << " GeV^2"
                      << ", Entries = " << sum_weights << "\n";
        }
    }

    // Compute Delta p_T^2
    TGraphErrors* g_Cu = new TGraphErrors();
    TGraphErrors* g_Sn = new TGraphErrors();
    TGraphErrors* g_CxC = new TGraphErrors();
    int pointCu = 0, pointSn = 0, pointCxC = 0;
    double maxY = -1e30, minY = 1e30;

    for (int kinBin = 0; kinBin < kinBins; ++kinBin) {
        double kinCenter = kinMin + (kinBin + 0.5) * (kinMax - kinMin) / kinBins;

        if (avg_pT2["Cu"][kinBin] > 0 && avg_pT2["LD2"][kinBin] > 0) {
            double delta_pT2 = avg_pT2["Cu"][kinBin] - avg_pT2["LD2"][kinBin];
            double err = sqrt(err_pT2["Cu"][kinBin] * err_pT2["Cu"][kinBin] + err_pT2["LD2"][kinBin] * err_pT2["LD2"][kinBin]);
            g_Cu->SetPoint(pointCu, kinCenter, delta_pT2);
            g_Cu->SetPointError(pointCu, 0, err);
            pointCu++;
            maxY = std::max(maxY, delta_pT2 + err);
            minY = std::min(minY, delta_pT2 - err);
        }

        if (avg_pT2["Sn"][kinBin] > 0 && avg_pT2["LD2"][kinBin] > 0) {
            double delta_pT2 = avg_pT2["Sn"][kinBin] - avg_pT2["LD2"][kinBin];
            double err = sqrt(err_pT2["Sn"][kinBin] * err_pT2["Sn"][kinBin] + err_pT2["LD2"][kinBin] * err_pT2["LD2"][kinBin]);
            g_Sn->SetPoint(pointSn, kinCenter, delta_pT2);
            g_Sn->SetPointError(pointSn, 0, err);
            pointSn++;
            maxY = std::max(maxY, delta_pT2 + err);
            minY = std::min(minY, delta_pT2 - err);
        }

        if (avg_pT2["CxC"][kinBin] > 0 && avg_pT2["LD2"][kinBin] > 0) {
            double delta_pT2 = avg_pT2["CxC"][kinBin] - avg_pT2["LD2"][kinBin];
            double err = sqrt(err_pT2["CxC"][kinBin] * err_pT2["CxC"][kinBin] + err_pT2["LD2"][kinBin] * err_pT2["LD2"][kinBin]);
            g_CxC->SetPoint(pointCxC, kinCenter, delta_pT2);
            g_CxC->SetPointError(pointCxC, 0, err);
            pointCxC++;
            maxY = std::max(maxY, delta_pT2 + err);
            minY = std::min(minY, delta_pT2 - err);
        }
    }

    // Plot
    TCanvas* canvas = new TCanvas("canvas", ("Transverse Momentum Broadening vs " + kin_var).c_str(), 800, 600);
    g_Cu->SetTitle(("Transverse Momentum Broadening vs " + kin_var).c_str());
    g_Cu->GetXaxis()->SetTitle(kinLabel.c_str());
    g_Cu->GetYaxis()->SetTitle(pion_type == "PP" ? "#Delta p_{T}^{2} (#pi^{+}) (GeV^{2})" : "#Delta p_{T}^{2} (#pi^{-}) (GeV^{2})");
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

    TLegend* legend = new TLegend(0.83, 0.75, 0.9, 0.9);
    legend->AddEntry(g_CxC, "CxC", "p");
    legend->AddEntry(g_Cu, "Cu", "p");
    legend->AddEntry(g_Sn, "Sn", "p");
    legend->SetTextSize(0.03);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();

    std::string output = "TransverseBroadening_" + config + "_" + pion_type + "_" + kin_var + ".png";
    canvas->SaveAs(output.c_str());

    return 0;
} 



  */


 /*  #include <cstdlib>
#include <iostream>
#include <filesystem>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <ctime>
#include <TAxis.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include <TH1F.h>
#include <TLine.h>
#include <TROOT.h>
#include "reader.h"

namespace fs = std::filesystem;

// Structure to hold particle data
struct ParticleData {
    float px, py, pz, p, vz;
    short status;
    int pid;
    float chi2pid;
};

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
    if (argc != 7) {
        std::cerr << "Usage: ./main <target1> <target2> <target3> <config> <pion_type> <kin_var>\n";
        std::cerr << "Example: ./main CxC LD2 CuSn OB PP z\n";
        std::cerr << "kin_var options: z, Q2, nu\n";
        return 1;
    }

    std::string target1 = argv[1];  // CxC
    std::string target2 = argv[2];  // LD2
    std::string target3 = argv[3];  // CuSn
    std::string config = argv[4];   // OB or IB
    std::string pion_type = argv[5]; // PP or NP
    std::string kin_var = argv[6];  // z, Q2, or nu

    if (config != "OB" && config != "IB") {
        std::cerr << "Invalid config. Use: OB or IB\n";
        return 1;
    }
    if (pion_type != "PP" && pion_type != "NP") {
        std::cerr << "Invalid pion type. Use: PP or NP\n";
        return 1;
    }
    if (kin_var != "z" && kin_var != "Q2" && kin_var != "nu") {
        std::cerr << "Invalid kinematic variable. Use: z, Q2, or nu\n";
        return 1;
    }

    // Generate unique suffix for histogram names to avoid redefinition
    std::time_t now = std::time(nullptr);
    std::string unique_suffix = std::to_string(now) + "_" + config + "_" + pion_type;

    double beamEnergy = 10.532; // GeV
    auto db = TDatabasePDG::Instance();
    double protonMass = db->GetParticle(2212)->Mass();
    TLorentzVector beam(0, 0, beamEnergy, beamEnergy);
    TLorentzVector target(0, 0, 0, protonMass);

    // Binning configuration for kinematic variable
    const int kinBins = 10;
    double kinMin, kinMax;
    std::string kinLabel;
    if (kin_var == "z") {
        kinMin = 0.3;
        kinMax = 0.7;
        kinLabel = "z";
    } else if (kin_var == "Q2") {
        kinMin = 1.0;
        kinMax = 9.0;
        kinLabel = "Q^{2} (GeV^{2})";
    } else { // nu
        kinMin = 0.0;
        kinMax = 10.0;
        kinLabel = "#nu (GeV)";
    }

    // Binning configuration for p_T^2
    const int pt2Bins = 10;
    const double pt2Min = 0.0;
    const double pt2Max = 1.2;

    // Store p_T^2 sums and counts per kinematic and p_T^2 bin
    std::map<std::string, std::vector<std::vector<double>>> pt2_sums; // target -> [kinBin][pt2Bin]
    std::map<std::string, std::vector<std::vector<double>>> pt2_counts; // target -> [kinBin][pt2Bin]
    std::vector<std::string> targets = {"LD2", "Cu", "Sn", "CxC"};

    // Histograms for kinematic variables
    std::map<std::string, TH1F*> h_electrons, h_pions;
    // Initialize histograms before the target loop to avoid redefinition
    for (const auto& tgt : targets) {
        pt2_sums[tgt] = std::vector<std::vector<double>>(kinBins, std::vector<double>(pt2Bins, 0.0));
        pt2_counts[tgt] = std::vector<std::vector<double>>(kinBins, std::vector<double>(pt2Bins, 0.0));
        h_electrons[tgt + "_Q2"] = new TH1F((tgt + "_Q2_" + unique_suffix).c_str(), (tgt + ": Q^{2};Q^{2} (GeV^{2});Counts").c_str(), 100, 0, 10);
        h_electrons[tgt + "_y"] = new TH1F((tgt + "_y_" + unique_suffix).c_str(), (tgt + ": y;y;Counts").c_str(), 100, 0, 1);
        h_electrons[tgt + "_W"] = new TH1F((tgt + "_W_" + unique_suffix).c_str(), (tgt + ": W;W (GeV);Counts").c_str(), 100, 0, 5);
        h_electrons[tgt + "_vz"] = new TH1F((tgt + "_vz_" + unique_suffix).c_str(), (tgt + ": v_{z};v_{z} (cm);Counts").c_str(), 100, -20, 10);
        h_pions[tgt + "_z"] = new TH1F((tgt + "_z_" + unique_suffix).c_str(), (tgt + ": z;z;Counts").c_str(), 100, 0, 1);
        h_pions[tgt + "_pt2"] = new TH1F((tgt + "_pt2_" + unique_suffix).c_str(), (tgt + ": p_{T}^{2};p_{T}^{2} (GeV^{2});Counts").c_str(), 100, 0, 2);
        h_pions[tgt + "_vz"] = new TH1F((tgt + "_vz_" + unique_suffix).c_str(), (tgt + ": v_{z};v_{z} (cm);Counts").c_str(), 100, -20, 10);
        // Register histograms with ROOT to manage them
        gROOT->Add(h_electrons[tgt + "_Q2"]);
        gROOT->Add(h_electrons[tgt + "_y"]);
        gROOT->Add(h_electrons[tgt + "_W"]);
        gROOT->Add(h_electrons[tgt + "_vz"]);
        gROOT->Add(h_pions[tgt + "_z"]);
        gROOT->Add(h_pions[tgt + "_pt2"]);
        gROOT->Add(h_pions[tgt + "_vz"]);
    }

    // Total electron counts
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
        {"Cu", {-10.5023, -6.4529}}, {"Sn", {-6.0086, 5.0000}}, {"LD2", {-20.0000, 5.0000}}, {"CxC", {-10.5319, 5.0000}}
    };
    std::map<std::string, std::pair<double, double>> vzCutsPiPOB = {
        {"Cu", {-9.76024, -5.1676}}, {"Sn", {-4.6413, 5.0000}}, {"LD2", {-20.0000, 5.0000}}, {"CxC", {-10.5319, 5.0000}}
    };
    std::map<std::string, std::pair<double, double>> vzCutsPiMOB = {
        {"Cu", {-10.9553, -5.8338}}, {"Sn", {-5.4629, 5.0000}}, {"LD2", {-20.0000, 5.0000}}, {"CxC", {-10.3501, 5.0000}}
    };

    // Placeholder vz cuts for inbending (IB)
    std::map<std::string, std::pair<double, double>> vzCutsEleIB = {{"Cu", {0, 0}}, {"Sn", {0, 0}}, {"LD2", {0, 0}}, {"CxC", {0, 0}}};
    std::map<std::string, std::pair<double, double>> vzCutsPiPIB = {{"Cu", {0, 0}}, {"Sn", {0, 0}}, {"LD2", {0, 0}}, {"CxC", {0, 0}}};
    std::map<std::string, std::pair<double, double>> vzCutsPiMIB = {{"Cu", {0, 0}}, {"Sn", {0, 0}}, {"LD2", {0, 0}}, {"CxC", {0, 0}}};

    auto vzCutsEle = (config == "OB") ? vzCutsEleOB : vzCutsEleIB;
    auto vzCutsPi = (pion_type == "PP") ? ((config == "OB") ? vzCutsPiPOB : vzCutsPiPIB)
                                       : ((config == "OB") ? vzCutsPiMOB : vzCutsPiMIB);
    int pion_pid = (pion_type == "PP") ? 211 : -211;

    std::vector<std::string> inputTargets = {target1, target2, target3};
    for (const auto& tgt : inputTargets) {
        if (targetDirs.find(tgt) == targetDirs.end()) continue;

        std::vector<std::string> hipoFiles = loadHipoFiles(targetDirs[tgt]);
        if (hipoFiles.empty()) continue;

        long totalElectrons = 0, triggerElectrons = 0, kinematicElectrons = 0;
        long totalPions = 0, triggerPions = 0, kinematicPions = 0;
        int event_count = 0;
        for (const auto& file : hipoFiles) {
            hipo::reader reader;
            reader.open(file.c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);
            hipo::bank RUN(factory.getSchema("RUN::config"));
            hipo::bank PART(factory.getSchema("REC::Particle"));
            hipo::event event;
            bool isFirstEvent = true;

            while (reader.next()) {
                event_count++;
                if (event_count > 1000000) break;
                reader.read(event);
                event.getStructure(PART);
                event.getStructure(RUN);

                if (isFirstEvent) {
                    float torus = RUN.getFloat("torus", 0);
                    bool isInbending = (torus < 0);
                    if ((config == "OB" && isInbending) || (config == "IB" && !isInbending)) break;
                    isFirstEvent = false;
                }

                std::vector<ParticleData> particles(PART.getRows());
                for (int i = 0; i < PART.getRows(); ++i) {
                    particles[i] = getParticleData(i, PART);
                }

                bool hasTrigger = false;
                TLorentzVector el_temp;
                double nu = 0.0;
                TLorentzVector q;

                for (const auto& p : particles) {
                    if (p.pid != 11) continue;
                    totalElectrons++;
                    if (h_electrons[tgt + "_vz"]) h_electrons[tgt + "_vz"]->Fill(p.vz);
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);
                        if (isCu) {
                            totalElectrons++;
                            if (h_electrons["Cu_vz"]) h_electrons["Cu_vz"]->Fill(p.vz);
                        }
                        if (isSn) {
                            totalElectrons++;
                            if (h_electrons["Sn_vz"]) h_electrons["Sn_vz"]->Fill(p.vz);
                        }
                    }

                    if (p.status >= 0 || abs(p.status) / 1000 != 2 || p.chi2pid < -5 || p.chi2pid > 5) continue;
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);
                        if (!isCu && !isSn) continue;
                    } else if (p.vz < vzCutsEle[tgt].first || p.vz > vzCutsEle[tgt].second) continue;

                    triggerElectrons++;
                    el_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(11)->Mass());
                    nu = beamEnergy - el_temp.E();
                    q = beam - el_temp;
                    double Q2 = -q.Mag2();
                    double y = nu / beamEnergy;
                    double W = sqrt(protonMass * protonMass + 2 * protonMass * nu - Q2);
                    if (h_electrons[tgt + "_Q2"]) h_electrons[tgt + "_Q2"]->Fill(Q2);
                    if (h_electrons[tgt + "_y"]) h_electrons[tgt + "_y"]->Fill(y);
                    if (h_electrons[tgt + "_W"]) h_electrons[tgt + "_W"]->Fill(W);
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);
                        if (isCu) {
                            if (h_electrons["Cu_Q2"]) h_electrons["Cu_Q2"]->Fill(Q2);
                            if (h_electrons["Cu_y"]) h_electrons["Cu_y"]->Fill(y);
                            if (h_electrons["Cu_W"]) h_electrons["Cu_W"]->Fill(W);
                        }
                        if (isSn) {
                            if (h_electrons["Sn_Q2"]) h_electrons["Sn_Q2"]->Fill(Q2);
                            if (h_electrons["Sn_y"]) h_electrons["Sn_y"]->Fill(y);
                            if (h_electrons["Sn_W"]) h_electrons["Sn_W"]->Fill(W);
                        }
                    }
                    if (Q2 <= 1 || y <= 0.25 || y >= 0.85 || W <= 2) continue;

                    kinematicElectrons++;
                    hasTrigger = true;
                    if (tgt == "LD2") total_LD2_ele++;
                    else if (tgt == "CxC") total_CxC_ele++;
                    else if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);
                        if (isCu) total_Cu_ele++;
                        if (isSn) total_Sn_ele++;
                    }
                    break;
                }

                if (!hasTrigger) continue;

                for (const auto& p : particles) {
                    if (p.pid != pion_pid) continue;
                    totalPions++;
                    if (h_pions[tgt + "_vz"]) h_pions[tgt + "_vz"]->Fill(p.vz);
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsPi["Cu"].first && p.vz <= vzCutsPi["Cu"].second);
                        bool isSn = (p.vz >= vzCutsPi["Sn"].first && p.vz <= vzCutsPi["Sn"].second);
                        if (isCu && h_pions["Cu_vz"]) h_pions["Cu_vz"]->Fill(p.vz);
                        if (isSn && h_pions["Sn_vz"]) h_pions["Sn_vz"]->Fill(p.vz);
                    }

                    if (abs(p.status) / 1000 != 2 || p.chi2pid < -10 || p.chi2pid > 10) continue;
                    TLorentzVector pi_temp;
                    pi_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(pion_pid)->Mass());
                    double z = pi_temp.E() / nu;
                    double pt2 = calculatePt2(pi_temp, q);
                    double Q2_val = -q.Mag2();
                    if (h_pions[tgt + "_z"]) h_pions[tgt + "_z"]->Fill(z);
                    if (h_pions[tgt + "_pt2"]) h_pions[tgt + "_pt2"]->Fill(pt2);
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsPi["Cu"].first && p.vz <= vzCutsPi["Cu"].second);
                        bool isSn = (p.vz >= vzCutsPi["Sn"].first && p.vz <= vzCutsPi["Sn"].second);
                        if (isCu) {
                            if (h_pions["Cu_z"]) h_pions["Cu_z"]->Fill(z);
                            if (h_pions["Cu_pt2"]) h_pions["Cu_pt2"]->Fill(pt2);
                        }
                        if (isSn) {
                            if (h_pions["Sn_z"]) h_pions["Sn_z"]->Fill(z);
                            if (h_pions["Sn_pt2"]) h_pions["Sn_pt2"]->Fill(pt2);
                        }
                    }

                    // Apply cuts
                    if (z <= 0.3 || z >= 0.7 || pt2 >= 1.2) continue;

                    triggerPions++;
                    double kin_value;
                    if (kin_var == "z") kin_value = z;
                    else if (kin_var == "Q2") kin_value = Q2_val;
                    else kin_value = nu;

                    // Apply kinematic variable cuts
                    if (kin_value < kinMin || kin_value > kinMax) continue;

                    int kinBin = static_cast<int>((kin_value - kinMin) / ((kinMax - kinMin) / kinBins));
                    if (kinBin < 0 || kinBin >= kinBins) continue;

                    int pt2Bin = static_cast<int>((pt2 - pt2Min) / ((pt2Max - pt2Min) / pt2Bins));
                    if (pt2Bin < 0 || pt2Bin >= pt2Bins) continue;

                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsPi["Cu"].first && p.vz <= vzCutsPi["Cu"].second);
                        bool isSn = (p.vz >= vzCutsPi["Sn"].first && p.vz <= vzCutsPi["Sn"].second);
                        if (!isCu && !isSn) continue;
                        if (isCu) {
                            pt2_sums["Cu"][kinBin][pt2Bin] += pt2;
                            pt2_counts["Cu"][kinBin][pt2Bin] += 1.0;
                        }
                        if (isSn) {
                            pt2_sums["Sn"][kinBin][pt2Bin] += pt2;
                            pt2_counts["Sn"][kinBin][pt2Bin] += 1.0;
                        }
                    } else {
                        if (p.vz < vzCutsPi[tgt].first || p.vz > vzCutsPi[tgt].second) continue;
                        pt2_sums[tgt][kinBin][pt2Bin] += pt2;
                        pt2_counts[tgt][kinBin][pt2Bin] += 1.0;
                    }
                    kinematicPions++;
                }
            }
        }

        std::cout << tgt << " Electron counts: Total=" << totalElectrons
                  << ", Trigger=" << triggerElectrons
                  << ", Kinematic=" << kinematicElectrons << "\n";
        std::cout << tgt << " Pion counts: Total=" << totalPions
                  << ", Trigger=" << triggerPions
                  << ", Kinematic=" << kinematicPions << "\n";
    }

    // Calculate <p_T^2> using weighted average over p_T^2 bins
    std::map<std::string, std::vector<double>> avg_pT2, err_pT2;
    for (const auto& tgt : targets) {
        avg_pT2[tgt].resize(kinBins, 0.0);
        err_pT2[tgt].resize(kinBins, 0.0);
        for (int kinBin = 0; kinBin < kinBins; ++kinBin) {
            double sum_weighted_pT2 = 0.0, sum_weights = 0.0, sum_pT2_sq = 0.0;
            for (int pt2Bin = 0; pt2Bin < pt2Bins; ++pt2Bin) {
                double N_j = pt2_counts[tgt][kinBin][pt2Bin];
                if (N_j <= 0) continue;
                double pt2_j = pt2_sums[tgt][kinBin][pt2Bin] / N_j; // Average p_T^2 in this sub-bin
                sum_weighted_pT2 += pt2_j * N_j; // Weighted contribution
                sum_weights += N_j; // Total weight (number of entries)
                sum_pT2_sq += N_j * pt2_j * pt2_j; // Weighted sum of (p_T^2)^2
            }
            if (sum_weights > 0) {
                avg_pT2[tgt][kinBin] = sum_weighted_pT2 / sum_weights;
                double avg_pT2_sq = sum_pT2_sq / sum_weights;
                double variance = avg_pT2_sq - avg_pT2[tgt][kinBin] * avg_pT2[tgt][kinBin];
                err_pT2[tgt][kinBin] = (variance >= 0) ? sqrt(variance / sum_weights) : 0.0;
            }

            double kinLow = kinMin + kinBin * (kinMax - kinMin) / kinBins;
            double kinHigh = kinMin + (kinBin + 1) * (kinMax - kinMin) / kinBins;
            std::cout << "Target " << tgt << ", " << kin_var << " Bin " << (kinBin + 1)
                      << " (" << kinLow << " to " << kinHigh << ")"
                      << ", <p_T^2> = " << avg_pT2[tgt][kinBin] << " GeV^2"
                      << ", Error = " << err_pT2[tgt][kinBin] << " GeV^2"
                      << ", Entries = " << sum_weights << "\n";
        }
    }

    // Compute Delta p_T^2
    TGraphErrors* g_Cu = new TGraphErrors();
    TGraphErrors* g_Sn = new TGraphErrors();
    TGraphErrors* g_CxC = new TGraphErrors();
    int pointCu = 0, pointSn = 0, pointCxC = 0;
    double maxY = -1e30, minY = 1e30;

    for (int kinBin = 0; kinBin < kinBins; ++kinBin) {
        double kinCenter = kinMin + (kinBin + 0.5) * (kinMax - kinMin) / kinBins;

        if (avg_pT2["Cu"][kinBin] > 0 && avg_pT2["LD2"][kinBin] > 0) {
            double delta_pT2 = avg_pT2["Cu"][kinBin] - avg_pT2["LD2"][kinBin];
            double err = sqrt(err_pT2["Cu"][kinBin] * err_pT2["Cu"][kinBin] + err_pT2["LD2"][kinBin] * err_pT2["LD2"][kinBin]);
            g_Cu->SetPoint(pointCu, kinCenter, delta_pT2);
            g_Cu->SetPointError(pointCu, 0, err);
            pointCu++;
            maxY = std::max(maxY, delta_pT2 + err);
            minY = std::min(minY, delta_pT2 - err);
        }

        if (avg_pT2["Sn"][kinBin] > 0 && avg_pT2["LD2"][kinBin] > 0) {
            double delta_pT2 = avg_pT2["Sn"][kinBin] - avg_pT2["LD2"][kinBin];
            double err = sqrt(err_pT2["Sn"][kinBin] * err_pT2["Sn"][kinBin] + err_pT2["LD2"][kinBin] * err_pT2["LD2"][kinBin]);
            g_Sn->SetPoint(pointSn, kinCenter, delta_pT2);
            g_Sn->SetPointError(pointSn, 0, err);
            pointSn++;
            maxY = std::max(maxY, delta_pT2 + err);
            minY = std::min(minY, delta_pT2 - err);
        }

        if (avg_pT2["CxC"][kinBin] > 0 && avg_pT2["LD2"][kinBin] > 0) {
            double delta_pT2 = avg_pT2["CxC"][kinBin] - avg_pT2["LD2"][kinBin];
            double err = sqrt(err_pT2["CxC"][kinBin] * err_pT2["CxC"][kinBin] + err_pT2["LD2"][kinBin] * err_pT2["LD2"][kinBin]);
            g_CxC->SetPoint(pointCxC, kinCenter, delta_pT2);
            g_CxC->SetPointError(pointCxC, 0, err);
            pointCxC++;
            maxY = std::max(maxY, delta_pT2 + err);
            minY = std::min(minY, delta_pT2 - err);
        }
    }

    // Plot transverse momentum broadening
    TCanvas* canvas = new TCanvas("canvas", ("Transverse Momentum Broadening vs " + kin_var).c_str(), 800, 600);
    canvas->SetLeftMargin(0.15);  // Try increasing this from the default (~0.1)
    canvas->SetBottomMargin(0.15); // Try increasing this from the default (~0.1)
    g_Cu->SetTitle(("Transverse Momentum Broadening vs " + kin_var).c_str());
    g_Cu->GetXaxis()->SetTitle(kinLabel.c_str());
    g_Cu->GetYaxis()->SetTitle(pion_type == "PP" ? "#Delta p_{T}^{2} (#pi^{+}) (GeV^{2})" : "#Delta p_{T}^{2} (#pi^{-}) (GeV^{2})");
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

    TLegend* legend = new TLegend(0.83, 0.75, 0.9, 0.9);
    legend->AddEntry(g_CxC, "CxC", "p");
    legend->AddEntry(g_Cu, "Cu", "p");
    legend->AddEntry(g_Sn, "Sn", "p");
    legend->SetTextSize(0.03);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();

    std::string output_pdf = "TransverseBroadening_" + config + "_" + pion_type + "_" + kin_var + ".pdf";
canvas->SaveAs(output_pdf.c_str());

    // Plot electron kinematics for each target separately
    for (const auto& tgt : targets) {
        TCanvas* c_electrons = new TCanvas(("c_electrons_" + tgt).c_str(), ("Electron Kinematics: " + tgt).c_str(), 1200, 800);
        c_electrons->Divide(2, 2);

        c_electrons->cd(1);
        if (h_electrons[tgt + "_Q2"]) {
            h_electrons[tgt + "_Q2"]->SetLineColor(kBlue);
            h_electrons[tgt + "_Q2"]->Draw("HIST");
            TLine* line_Q2 = new TLine(1.0, 0, 1.0, h_electrons[tgt + "_Q2"]->GetMaximum());
            line_Q2->SetLineColor(kBlack);
            line_Q2->SetLineStyle(2);
            line_Q2->Draw();
        }

        c_electrons->cd(2);
        if (h_electrons[tgt + "_y"]) {
            h_electrons[tgt + "_y"]->SetLineColor(kBlue);
            h_electrons[tgt + "_y"]->Draw("HIST");
            TLine* line_y1 = new TLine(0.25, 0, 0.25, h_electrons[tgt + "_y"]->GetMaximum());
            TLine* line_y2 = new TLine(0.85, 0, 0.85, h_electrons[tgt + "_y"]->GetMaximum());
            line_y1->SetLineColor(kBlack);
            line_y2->SetLineColor(kBlack);
            line_y1->SetLineStyle(2);
            line_y2->SetLineStyle(2);
            line_y1->Draw();
            line_y2->Draw();
        }

        c_electrons->cd(3);
        if (h_electrons[tgt + "_W"]) {
            h_electrons[tgt + "_W"]->SetLineColor(kBlue);
            h_electrons[tgt + "_W"]->Draw("HIST");
            TLine* line_W = new TLine(2.0, 0, 2.0, h_electrons[tgt + "_W"]->GetMaximum());
            line_W->SetLineColor(kBlack);
            line_W->SetLineStyle(2);
            line_W->Draw();
        }

        c_electrons->cd(4);
        if (h_electrons[tgt + "_vz"]) {
            h_electrons[tgt + "_vz"]->SetLineColor(kBlue);
            h_electrons[tgt + "_vz"]->Draw("HIST");
            TLine* line_vz1 = new TLine(vzCutsEle[tgt].first, 0, vzCutsEle[tgt].first, h_electrons[tgt + "_vz"]->GetMaximum());
            TLine* line_vz2 = new TLine(vzCutsEle[tgt].second, 0, vzCutsEle[tgt].second, h_electrons[tgt + "_vz"]->GetMaximum());
            line_vz1->SetLineColor(kBlack);
            line_vz2->SetLineColor(kBlack);
            line_vz1->SetLineStyle(2);
            line_vz2->SetLineStyle(2);
            if (config == "OB") {
                line_vz1->Draw();
                line_vz2->Draw();
            }
        }

        c_electrons->SaveAs(("Electron_Kinematics_" + config + "_" + pion_type + "_" + tgt + ".pdf").c_str());
        delete c_electrons;
    }

    // Plot pion kinematics for each target separately
    for (const auto& tgt : targets) {
        TCanvas* c_pions = new TCanvas(("c_pions_" + tgt).c_str(), ("Pion Kinematics: " + tgt).c_str(), 1200, 800);
        c_pions->Divide(2, 2);

        c_pions->cd(1);
        if (h_pions[tgt + "_z"]) {
            h_pions[tgt + "_z"]->SetLineColor(kBlue);
            h_pions[tgt + "_z"]->Draw("HIST");
            TLine* line_z1 = new TLine(0.3, 0, 0.3, h_pions[tgt + "_z"]->GetMaximum());
            TLine* line_z2 = new TLine(0.7, 0, 0.7, h_pions[tgt + "_z"]->GetMaximum());
            line_z1->SetLineColor(kBlack);
            line_z2->SetLineColor(kBlack);
            line_z1->SetLineStyle(2);
            line_z2->SetLineStyle(2);
            line_z1->Draw();
            line_z2->Draw();
        }

        c_pions->cd(2);
        if (h_pions[tgt + "_pt2"]) {
            h_pions[tgt + "_pt2"]->SetLineColor(kBlue);
            h_pions[tgt + "_pt2"]->Draw("HIST");
            TLine* line_pt2 = new TLine(1.2, 0, 1.2, h_pions[tgt + "_pt2"]->GetMaximum());
            line_pt2->SetLineColor(kBlack);
            line_pt2->SetLineStyle(2);
            line_pt2->Draw();
        }

        c_pions->cd(3);
        if (h_pions[tgt + "_vz"]) {
            h_pions[tgt + "_vz"]->SetLineColor(kBlue);
            h_pions[tgt + "_vz"]->Draw("HIST");
            TLine* line_pivz1 = new TLine(vzCutsPi[tgt].first, 0, vzCutsPi[tgt].first, h_pions[tgt + "_vz"]->GetMaximum());
            TLine* line_pivz2 = new TLine(vzCutsPi[tgt].second, 0, vzCutsPi[tgt].second, h_pions[tgt + "_vz"]->GetMaximum());
            line_pivz1->SetLineColor(kBlack);
            line_pivz2->SetLineColor(kBlack);
            line_pivz1->SetLineStyle(2);
            line_pivz2->SetLineStyle(2);
            if (config == "OB") {
                line_pivz1->Draw();
                line_pivz2->Draw();
            }
        }

        c_pions->SaveAs(("Pion_Kinematics_" + config + "_" + pion_type + "_" + tgt + ".pdf").c_str());
        delete c_pions;
    }

    // Save histograms to ROOT file
    TFile* f_out = new TFile(("Kinematics_" + config + "_" + pion_type + ".root").c_str(), "RECREATE");
    if (f_out && f_out->IsOpen()) {
        for (const auto& tgt : targets) {
            if (h_electrons[tgt + "_Q2"]) h_electrons[tgt + "_Q2"]->Write();
            if (h_electrons[tgt + "_y"]) h_electrons[tgt + "_y"]->Write();
            if (h_electrons[tgt + "_W"]) h_electrons[tgt + "_W"]->Write();
            if (h_electrons[tgt + "_vz"]) h_electrons[tgt + "_vz"]->Write();
            if (h_pions[tgt + "_z"]) h_pions[tgt + "_z"]->Write();
            if (h_pions[tgt + "_pt2"]) h_pions[tgt + "_pt2"]->Write();
            if (h_pions[tgt + "_vz"]) h_pions[tgt + "_vz"]->Write();
        }
        if (canvas) canvas->Write();
        f_out->Close();
        delete f_out;
    }

    // Clean up histograms
    for (const auto& pair : h_electrons) {
        if (pair.second) {
            gROOT->Remove(pair.second);
            delete pair.second;
        }
    }
    for (const auto& pair : h_pions) {
        if (pair.second) {
            gROOT->Remove(pair.second);
            delete pair.second;
        }
    }
    delete g_Cu;
    delete g_Sn;
    delete g_CxC;
    delete canvas;
    delete legend;

    return 0;
} */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <ctime>
#include <TAxis.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include <TH1F.h>
#include <TLine.h>
#include <TROOT.h>
#include "reader.h"

namespace fs = std::filesystem;

// Structure to hold particle data
struct ParticleData {
    float px, py, pz, p, vz;
    short status;
    int pid;
    float chi2pid;
};

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
    if (argc != 7) {
        std::cerr << "Usage: ./main <target1> <target2> <target3> <config> <pion_type> <kin_var>\n";
        std::cerr << "Example: ./main CxC LD2 CuSn OB PP z\n";
        std::cerr << "kin_var options: z, Q2, nu\n";
        return 1;
    }

    std::string target1 = argv[1];  // CxC
    std::string target2 = argv[2];  // LD2
    std::string target3 = argv[3];  // CuSn
    std::string config = argv[4];   // OB or IB
    std::string pion_type = argv[5]; // PP or NP
    std::string kin_var = argv[6];  // z, Q2, or nu

    if (config != "OB" && config != "IB") {
        std::cerr << "Invalid config. Use: OB or IB\n";
        return 1;
    }
    if (pion_type != "PP" && pion_type != "NP") {
        std::cerr << "Invalid pion type. Use: PP or NP\n";
        return 1;
    }
    if (kin_var != "z" && kin_var != "Q2" && kin_var != "nu") {
        std::cerr << "Invalid kinematic variable. Use: z, Q2, or nu\n";
        return 1;
    }

    // Generate unique suffix for histogram names to avoid redefinition
    std::time_t now = std::time(nullptr);
    std::string unique_suffix = std::to_string(now) + "_" + config + "_" + pion_type;

    // Open a text file to log processed .hipo files
    std::ofstream logFile("processed_hipo_files.txt", std::ios::out);
    if (!logFile.is_open()) {
        std::cerr << "Warning: Could not open processed_hipo_files.txt for writing\n";
    }

    double beamEnergy = 10.532; // GeV
    auto db = TDatabasePDG::Instance();
    double protonMass = db->GetParticle(2212)->Mass();
    TLorentzVector beam(0, 0, beamEnergy, beamEnergy);
    TLorentzVector target(0, 0, 0, protonMass);

    // Binning configuration for kinematic variable
    const int kinBins = 10;
    double kinMin, kinMax;
    std::string kinLabel;
    if (kin_var == "z") {
        kinMin = 0.3;
        kinMax = 0.7;
        kinLabel = "z";
    } else if (kin_var == "Q2") {
        kinMin = 1.0;
        kinMax = 9.0;
        kinLabel = "Q^{2} (GeV^{2})";
    } else { // nu
        kinMin = 0.0;
        kinMax = 10.0;
        kinLabel = "#nu (GeV)";
    }

    // Binning configuration for p_T^2
    const int pt2Bins = 10;
    const double pt2Min = 0.0;
    const double pt2Max = 1.2;

    // Store p_T^2 sums and counts per kinematic and p_T^2 bin
    std::map<std::string, std::vector<std::vector<double>>> pt2_sums; // target -> [kinBin][pt2Bin]
    std::map<std::string, std::vector<std::vector<double>>> pt2_counts; // target -> [kinBin][pt2Bin]
    std::vector<std::string> targets = {"LD2", "Cu", "Sn", "CxC"};

    // Histograms for kinematic variables
    std::map<std::string, TH1F*> h_electrons, h_pions;
    // Initialize histograms before the target loop to avoid redefinition
    for (const auto& tgt : targets) {
        pt2_sums[tgt] = std::vector<std::vector<double>>(kinBins, std::vector<double>(pt2Bins, 0.0));
        pt2_counts[tgt] = std::vector<std::vector<double>>(kinBins, std::vector<double>(pt2Bins, 0.0));
        h_electrons[tgt + "_Q2"] = new TH1F((tgt + "_Q2_" + unique_suffix).c_str(), (tgt + ": Q^{2};Q^{2} (GeV^{2});Counts").c_str(), 100, 0, 10);
        h_electrons[tgt + "_y"] = new TH1F((tgt + "_y_" + unique_suffix).c_str(), (tgt + ": y;y;Counts").c_str(), 100, 0, 1);
        h_electrons[tgt + "_W"] = new TH1F((tgt + "_W_" + unique_suffix).c_str(), (tgt + ": W;W (GeV);Counts").c_str(), 100, 0, 5);
        h_electrons[tgt + "_vz"] = new TH1F((tgt + "_vz_" + unique_suffix).c_str(), (tgt + ": v_{z};v_{z} (cm);Counts").c_str(), 100, -20, 10);
        h_pions[tgt + "_z"] = new TH1F((tgt + "_z_" + unique_suffix).c_str(), (tgt + ": z;z;Counts").c_str(), 100, 0, 1);
        h_pions[tgt + "_pt2"] = new TH1F((tgt + "_pt2_" + unique_suffix).c_str(), (tgt + ": p_{T}^{2};p_{T}^{2} (GeV^{2});Counts").c_str(), 100, 0, 2);
        h_pions[tgt + "_vz"] = new TH1F((tgt + "_vz_" + unique_suffix).c_str(), (tgt + ": v_{z};v_{z} (cm);Counts").c_str(), 100, -20, 10);
        // Register histograms with ROOT to manage them
        gROOT->Add(h_electrons[tgt + "_Q2"]);
        gROOT->Add(h_electrons[tgt + "_y"]);
        gROOT->Add(h_electrons[tgt + "_W"]);
        gROOT->Add(h_electrons[tgt + "_vz"]);
        gROOT->Add(h_pions[tgt + "_z"]);
        gROOT->Add(h_pions[tgt + "_pt2"]);
        gROOT->Add(h_pions[tgt + "_vz"]);
    }

    // Total electron counts
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
        {"Cu", {-10.5023, -6.4529}}, {"Sn", {-6.0086, 5.0000}}, {"LD2", {-20.0000, 5.0000}}, {"CxC", {-10.5319, 5.0000}}
    };
    std::map<std::string, std::pair<double, double>> vzCutsPiPOB = {
        {"Cu", {-9.76024, -5.1676}}, {"Sn", {-4.6413, 5.0000}}, {"LD2", {-20.0000, 5.0000}}, {"CxC", {-10.5319, 5.0000}}
    };
    std::map<std::string, std::pair<double, double>> vzCutsPiMOB = {
        {"Cu", {-10.9553, -5.8338}}, {"Sn", {-5.4629, 5.0000}}, {"LD2", {-20.0000, 5.0000}}, {"CxC", {-10.3501, 5.0000}}
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

    auto vzCutsEle = (config == "OB") ? vzCutsEleOB : vzCutsEleIB;
    auto vzCutsPi = (pion_type == "PP") ? ((config == "OB") ? vzCutsPiPOB : vzCutsPiPIB)
                                       : ((config == "OB") ? vzCutsPiMOB : vzCutsPiMIB);
    int pion_pid = (pion_type == "PP") ? 211 : -211;

    const long maxEvents = 100000000; // Total event limit per target

    std::vector<std::string> inputTargets = {target1, target2, target3};
    for (const auto& tgt : inputTargets) {
        if (targetDirs.find(tgt) == targetDirs.end()) continue;

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

        long totalElectrons = 0, triggerElectrons = 0, kinematicElectrons = 0;
        long totalPions = 0, triggerPions = 0, kinematicPions = 0;
        long event_count = 0;

        for (const auto& file : hipoFiles) {
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
                if (event_count >= maxEvents) break; // Stop if total event limit reached
                event_count++;
                reader.read(event);
                event.getStructure(PART);
                event.getStructure(RUN);

                if (isFirstEvent) {
                    float torus = RUN.getFloat("torus", 0);
                    bool isInbending = (torus < 0);
                    if ((config == "OB" && isInbending) || (config == "IB" && !isInbending)) break;
                    isFirstEvent = false;
                }

                std::vector<ParticleData> particles(PART.getRows());
                for (int i = 0; i < PART.getRows(); ++i) {
                    particles[i] = getParticleData(i, PART);
                }

                bool hasTrigger = false;
                TLorentzVector el_temp;
                double nu = 0.0;
                TLorentzVector q;

                for (const auto& p : particles) {
                    if (p.pid != 11) continue;
                    totalElectrons++;
                    if (h_electrons[tgt + "_vz"]) h_electrons[tgt + "_vz"]->Fill(p.vz);
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);
                        if (isCu) {
                            totalElectrons++;
                            if (h_electrons["Cu_vz"]) h_electrons["Cu_vz"]->Fill(p.vz);
                        }
                        if (isSn) {
                            totalElectrons++;
                            if (h_electrons["Sn_vz"]) h_electrons["Sn_vz"]->Fill(p.vz);
                        }
                    }

                    if (p.status >= 0 || abs(p.status) / 1000 != 2 || p.chi2pid < -5 || p.chi2pid > 5) continue;
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);
                        if (!isCu && !isSn) continue;
                    } else if (p.vz < vzCutsEle[tgt].first || p.vz > vzCutsEle[tgt].second) continue;

                    triggerElectrons++;
                    el_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(11)->Mass());
                    nu = beamEnergy - el_temp.E();
                    q = beam - el_temp;
                    double Q2 = -q.Mag2();
                    double y = nu / beamEnergy;
                    double W = sqrt(protonMass * protonMass + 2 * protonMass * nu - Q2);
                    if (h_electrons[tgt + "_Q2"]) h_electrons[tgt + "_Q2"]->Fill(Q2);
                    if (h_electrons[tgt + "_y"]) h_electrons[tgt + "_y"]->Fill(y);
                    if (h_electrons[tgt + "_W"]) h_electrons[tgt + "_W"]->Fill(W);
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);
                        if (isCu) {
                            if (h_electrons["Cu_Q2"]) h_electrons["Cu_Q2"]->Fill(Q2);
                            if (h_electrons["Cu_y"]) h_electrons["Cu_y"]->Fill(y);
                            if (h_electrons["Cu_W"]) h_electrons["Cu_W"]->Fill(W);
                        }
                        if (isSn) {
                            if (h_electrons["Sn_Q2"]) h_electrons["Sn_Q2"]->Fill(Q2);
                            if (h_electrons["Sn_y"]) h_electrons["Sn_y"]->Fill(y);
                            if (h_electrons["Sn_W"]) h_electrons["Sn_W"]->Fill(W);
                        }
                    }
                    if (Q2 <= 1 || y <= 0.25 || y >= 0.85 || W <= 2) continue;

                    kinematicElectrons++;
                    hasTrigger = true;
                    if (tgt == "LD2") total_LD2_ele++;
                    else if (tgt == "CxC") total_CxC_ele++;
                    else if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        bool isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);
                        if (isCu) total_Cu_ele++;
                        if (isSn) total_Sn_ele++;
                    }
                    break;
                }

                if (!hasTrigger) continue;

                for (const auto& p : particles) {
                    if (p.pid != pion_pid) continue;
                    totalPions++;
                    if (h_pions[tgt + "_vz"]) h_pions[tgt + "_vz"]->Fill(p.vz);
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsPi["Cu"].first && p.vz <= vzCutsPi["Cu"].second);
                        bool isSn = (p.vz >= vzCutsPi["Sn"].first && p.vz <= vzCutsPi["Sn"].second);
                        if (isCu && h_pions["Cu_vz"]) h_pions["Cu_vz"]->Fill(p.vz);
                        if (isSn && h_pions["Sn_vz"]) h_pions["Sn_vz"]->Fill(p.vz);
                    }

                    if (abs(p.status) / 1000 != 2 || p.chi2pid < -10 || p.chi2pid > 10) continue;
                    TLorentzVector pi_temp;
                    pi_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(pion_pid)->Mass());
                    double z = pi_temp.E() / nu;
                    double pt2 = calculatePt2(pi_temp, q);
                    double Q2_val = -q.Mag2();
                    if (h_pions[tgt + "_z"]) h_pions[tgt + "_z"]->Fill(z);
                    if (h_pions[tgt + "_pt2"]) h_pions[tgt + "_pt2"]->Fill(pt2);
                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsPi["Cu"].first && p.vz <= vzCutsPi["Cu"].second);
                        bool isSn = (p.vz >= vzCutsPi["Sn"].first && p.vz <= vzCutsPi["Sn"].second);
                        if (isCu) {
                            if (h_pions["Cu_z"]) h_pions["Cu_z"]->Fill(z);
                            if (h_pions["Cu_pt2"]) h_pions["Cu_pt2"]->Fill(pt2);
                        }
                        if (isSn) {
                            if (h_pions["Sn_z"]) h_pions["Sn_z"]->Fill(z);
                            if (h_pions["Sn_pt2"]) h_pions["Sn_pt2"]->Fill(pt2);
                        }
                    }

                    // Apply cuts
                    if (z <= 0.3 || z >= 0.7 || pt2 >= 1.2) continue;

                    triggerPions++;
                    double kin_value;
                    if (kin_var == "z") kin_value = z;
                    else if (kin_var == "Q2") kin_value = Q2_val;
                    else kin_value = nu;

                    // Apply kinematic variable cuts
                    if (kin_value < kinMin || kin_value > kinMax) continue;

                    int kinBin = static_cast<int>((kin_value - kinMin) / ((kinMax - kinMin) / kinBins));
                    if (kinBin < 0 || kinBin >= kinBins) continue;

                    int pt2Bin = static_cast<int>((pt2 - pt2Min) / ((pt2Max - pt2Min) / pt2Bins));
                    if (pt2Bin < 0 || pt2Bin >= pt2Bins) continue;

                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsPi["Cu"].first && p.vz <= vzCutsPi["Cu"].second);
                        bool isSn = (p.vz >= vzCutsPi["Sn"].first && p.vz <= vzCutsPi["Sn"].second);
                        if (!isCu && !isSn) continue;
                        if (isCu) {
                            pt2_sums["Cu"][kinBin][pt2Bin] += pt2;
                            pt2_counts["Cu"][kinBin][pt2Bin] += 1.0;
                        }
                        if (isSn) {
                            pt2_sums["Sn"][kinBin][pt2Bin] += pt2;
                            pt2_counts["Sn"][kinBin][pt2Bin] += 1.0;
                        }
                    } else {
                        if (p.vz < vzCutsPi[tgt].first || p.vz > vzCutsPi[tgt].second) continue;
                        pt2_sums[tgt][kinBin][pt2Bin] += pt2;
                        pt2_counts[tgt][kinBin][pt2Bin] += 1.0;
                    }
                    kinematicPions++;
                }
            }
            if (event_count >= maxEvents) break; // Stop processing more files if limit reached
        }

        std::cout << tgt << " Electron counts: Total=" << totalElectrons
                  << ", Trigger=" << triggerElectrons
                  << ", Kinematic=" << kinematicElectrons << "\n";
        std::cout << tgt << " Pion counts: Total=" << totalPions
                  << ", Trigger=" << triggerPions
                  << ", Kinematic=" << kinematicPions << "\n";
        std::cout << tgt << " Total events processed: " << event_count << "\n";
        // Add a separator in the log file
        if (logFile.is_open()) {
            logFile << "Total events processed: " << event_count << "\n\n";
        }
    }

    // Close the log file
    if (logFile.is_open()) {
        logFile.close();
    }

    // Calculate <p_T^2> using weighted average over p_T^2 bins
    std::map<std::string, std::vector<double>> avg_pT2, err_pT2;
    for (const auto& tgt : targets) {
        avg_pT2[tgt].resize(kinBins, 0.0);
        err_pT2[tgt].resize(kinBins, 0.0);
        for (int kinBin = 0; kinBin < kinBins; ++kinBin) {
            double sum_weighted_pT2 = 0.0, sum_weights = 0.0, sum_pT2_sq = 0.0;
            for (int pt2Bin = 0; pt2Bin < pt2Bins; ++pt2Bin) {
                double N_j = pt2_counts[tgt][kinBin][pt2Bin];
                if (N_j <= 0) continue;
                double pt2_j = pt2_sums[tgt][kinBin][pt2Bin] / N_j; // Average p_T^2 in this sub-bin
                sum_weighted_pT2 += pt2_j * N_j; // Weighted contribution
                sum_weights += N_j; // Total weight (number of entries)
                sum_pT2_sq += N_j * pt2_j * pt2_j; // Weighted sum of (p_T^2)^2
            }
            if (sum_weights > 0) {
                avg_pT2[tgt][kinBin] = sum_weighted_pT2 / sum_weights;
                double avg_pT2_sq = sum_pT2_sq / sum_weights;
                double variance = avg_pT2_sq - avg_pT2[tgt][kinBin] * avg_pT2[tgt][kinBin];
                err_pT2[tgt][kinBin] = (variance >= 0) ? sqrt(variance / sum_weights) : 0.0;
            }

            double kinLow = kinMin + kinBin * (kinMax - kinMin) / kinBins;
            double kinHigh = kinMin + (kinBin + 1) * (kinMax - kinMin) / kinBins;
            std::cout << "Target " << tgt << ", " << kin_var << " Bin " << (kinBin + 1)
                      << " (" << kinLow << " to " << kinHigh << ")"
                      << ", <p_T^2> = " << avg_pT2[tgt][kinBin] << " GeV^2"
                      << ", Error = " << err_pT2[tgt][kinBin] << " GeV^2"
                      << ", Entries = " << sum_weights << "\n";
        }
    }

    // Compute Delta p_T^2
    TGraphErrors* g_Cu = new TGraphErrors();
    TGraphErrors* g_Sn = new TGraphErrors();
    TGraphErrors* g_CxC = new TGraphErrors();
    int pointCu = 0, pointSn = 0, pointCxC = 0;
    double maxY = -1e30, minY = 1e30;

    for (int kinBin = 0; kinBin < kinBins; ++kinBin) {
        double kinCenter = kinMin + (kinBin + 0.5) * (kinMax - kinMin) / kinBins;

        if (avg_pT2["Cu"][kinBin] > 0 && avg_pT2["LD2"][kinBin] > 0) {
            double delta_pT2 = avg_pT2["Cu"][kinBin] - avg_pT2["LD2"][kinBin];
            double err = sqrt(err_pT2["Cu"][kinBin] * err_pT2["Cu"][kinBin] + err_pT2["LD2"][kinBin] * err_pT2["LD2"][kinBin]);
            g_Cu->SetPoint(pointCu, kinCenter, delta_pT2);
            g_Cu->SetPointError(pointCu, 0, err);
            pointCu++;
            maxY = std::max(maxY, delta_pT2 + err);
            minY = std::min(minY, delta_pT2 - err);
        }

        if (avg_pT2["Sn"][kinBin] > 0 && avg_pT2["LD2"][kinBin] > 0) {
            double delta_pT2 = avg_pT2["Sn"][kinBin] - avg_pT2["LD2"][kinBin];
            double err = sqrt(err_pT2["Sn"][kinBin] * err_pT2["Sn"][kinBin] + err_pT2["LD2"][kinBin] * err_pT2["LD2"][kinBin]);
            g_Sn->SetPoint(pointSn, kinCenter, delta_pT2);
            g_Sn->SetPointError(pointSn, 0, err);
            pointSn++;
            maxY = std::max(maxY, delta_pT2 + err);
            minY = std::min(minY, delta_pT2 - err);
        }

        if (avg_pT2["CxC"][kinBin] > 0 && avg_pT2["LD2"][kinBin] > 0) {
            double delta_pT2 = avg_pT2["CxC"][kinBin] - avg_pT2["LD2"][kinBin];
            double err = sqrt(err_pT2["CxC"][kinBin] * err_pT2["CxC"][kinBin] + err_pT2["LD2"][kinBin] * err_pT2["LD2"][kinBin]);
            g_CxC->SetPoint(pointCxC, kinCenter, delta_pT2);
            g_CxC->SetPointError(pointCxC, 0, err);
            pointCxC++;
            maxY = std::max(maxY, delta_pT2 + err);
            minY = std::min(minY, delta_pT2 - err);
        }
    }

    // Plot transverse momentum broadening
    TCanvas* canvas = new TCanvas("canvas", ("Transverse Momentum Broadening vs " + kin_var).c_str(), 800, 600);
    canvas->SetLeftMargin(0.15);
    canvas->SetBottomMargin(0.15);
    g_Cu->SetTitle(("Transverse Momentum Broadening vs " + kin_var).c_str());
    g_Cu->GetXaxis()->SetTitle(kinLabel.c_str());
    g_Cu->GetYaxis()->SetTitle(pion_type == "PP" ? "#Delta p_{T}^{2} (#pi^{+}) (GeV^{2})" : "#Delta p_{T}^{2} (#pi^{-}) (GeV^{2})");
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

    TLegend* legend = new TLegend(0.83, 0.75, 0.9, 0.9);
    legend->AddEntry(g_CxC, "CxC", "p");
    legend->AddEntry(g_Cu, "Cu", "p");
    legend->AddEntry(g_Sn, "Sn", "p");
    legend->SetTextSize(0.03);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();

    std::string output = "TransverseBroadening_" + config + "_" + pion_type + "_" + kin_var + ".png";
    canvas->SaveAs(output.c_str());

    // Plot electron kinematics for each target separately
    for (const auto& tgt : targets) {
        TCanvas* c_electrons = new TCanvas(("c_electrons_" + tgt).c_str(), ("Electron Kinematics: " + tgt).c_str(), 1200, 800);
        c_electrons->Divide(2, 2);

        c_electrons->cd(1);
        if (h_electrons[tgt + "_Q2"]) {
            h_electrons[tgt + "_Q2"]->SetLineColor(kBlue);
            h_electrons[tgt + "_Q2"]->Draw("HIST");
            TLine* line_Q2 = new TLine(1.0, 0, 1.0, h_electrons[tgt + "_Q2"]->GetMaximum());
            line_Q2->SetLineColor(kBlack);
            line_Q2->SetLineStyle(2);
            line_Q2->Draw();
        }

        c_electrons->cd(2);
        if (h_electrons[tgt + "_y"]) {
            h_electrons[tgt + "_y"]->SetLineColor(kBlue);
            h_electrons[tgt + "_y"]->Draw("HIST");
            TLine* line_y1 = new TLine(0.25, 0, 0.25, h_electrons[tgt + "_y"]->GetMaximum());
            TLine* line_y2 = new TLine(0.85, 0, 0.85, h_electrons[tgt + "_y"]->GetMaximum());
            line_y1->SetLineColor(kBlack);
            line_y2->SetLineColor(kBlack);
            line_y1->SetLineStyle(2);
            line_y2->SetLineStyle(2);
            line_y1->Draw();
            line_y2->Draw();
        }

        c_electrons->cd(3);
        if (h_electrons[tgt + "_W"]) {
            h_electrons[tgt + "_W"]->SetLineColor(kBlue);
            h_electrons[tgt + "_W"]->Draw("HIST");
            TLine* line_W = new TLine(2.0, 0, 2.0, h_electrons[tgt + "_W"]->GetMaximum());
            line_W->SetLineColor(kBlack);
            line_W->SetLineStyle(2);
            line_W->Draw();
        }

        c_electrons->cd(4);
        if (h_electrons[tgt + "_vz"]) {
            h_electrons[tgt + "_vz"]->SetLineColor(kBlue);
            h_electrons[tgt + "_vz"]->Draw("HIST");
            TLine* line_vz1 = new TLine(vzCutsEle[tgt].first, 0, vzCutsEle[tgt].first, h_electrons[tgt + "_vz"]->GetMaximum());
            TLine* line_vz2 = new TLine(vzCutsEle[tgt].second, 0, vzCutsEle[tgt].second, h_electrons[tgt + "_vz"]->GetMaximum());
            line_vz1->SetLineColor(kBlack);
            line_vz2->SetLineColor(kBlack);
            line_vz1->SetLineStyle(2);
            line_vz2->SetLineStyle(2);
            if (config == "OB") {
                line_vz1->Draw();
                line_vz2->Draw();
            }
        }

        c_electrons->SaveAs(("Electron_Kinematics_" + config + "_" + pion_type + "_" + tgt + ".pdf").c_str());
        delete c_electrons;
    }

    // Plot pion kinematics for each target separately
    for (const auto& tgt : targets) {
        TCanvas* c_pions = new TCanvas(("c_pions_" + tgt).c_str(), ("Pion Kinematics: " + tgt).c_str(), 1200, 800);
        c_pions->Divide(2, 2);

        c_pions->cd(1);
        if (h_pions[tgt + "_z"]) {
            h_pions[tgt + "_z"]->SetLineColor(kBlue);
            h_pions[tgt + "_z"]->Draw("HIST");
            TLine* line_z1 = new TLine(0.3, 0, 0.3, h_pions[tgt + "_z"]->GetMaximum());
            TLine* line_z2 = new TLine(0.7, 0, 0.7, h_pions[tgt + "_z"]->GetMaximum());
            line_z1->SetLineColor(kBlack);
            line_z2->SetLineColor(kBlack);
            line_z1->SetLineStyle(2);
            line_z2->SetLineStyle(2);
            line_z1->Draw();
            line_z2->Draw();
        }

        c_pions->cd(2);
        if (h_pions[tgt + "_pt2"]) {
            h_pions[tgt + "_pt2"]->SetLineColor(kBlue);
            h_pions[tgt + "_pt2"]->Draw("HIST");
            TLine* line_pt2 = new TLine(1.2, 0, 1.2, h_pions[tgt + "_pt2"]->GetMaximum());
            line_pt2->SetLineColor(kBlack);
            line_pt2->SetLineStyle(2);
            line_pt2->Draw();
        }

        c_pions->cd(3);
        if (h_pions[tgt + "_vz"]) {
            h_pions[tgt + "_vz"]->SetLineColor(kBlue);
            h_pions[tgt + "_vz"]->Draw("HIST");
            TLine* line_pivz1 = new TLine(vzCutsPi[tgt].first, 0, vzCutsPi[tgt].first, h_pions[tgt + "_vz"]->GetMaximum());
            TLine* line_pivz2 = new TLine(vzCutsPi[tgt].second, 0, vzCutsPi[tgt].second, h_pions[tgt + "_vz"]->GetMaximum());
            line_pivz1->SetLineColor(kBlack);
            line_pivz2->SetLineColor(kBlack);
            line_pivz1->SetLineStyle(2);
            line_pivz2->SetLineStyle(2);
            if (config == "OB") {
                line_pivz1->Draw();
                line_pivz2->Draw();
            }
        }

        c_pions->SaveAs(("Pion_Kinematics_" + config + "_" + pion_type + "_" + tgt + ".pdf").c_str());
        delete c_pions;
    }

    // Save histograms to ROOT file
    TFile* f_out = new TFile(("Kinematics_" + config + "_" + pion_type + ".root").c_str(), "RECREATE");
    if (f_out && f_out->IsOpen()) {
        for (const auto& tgt : targets) {
            if (h_electrons[tgt + "_Q2"]) h_electrons[tgt + "_Q2"]->Write();
            if (h_electrons[tgt + "_y"]) h_electrons[tgt + "_y"]->Write();
            if (h_electrons[tgt + "_W"]) h_electrons[tgt + "_W"]->Write();
            if (h_electrons[tgt + "_vz"]) h_electrons[tgt + "_vz"]->Write();
            if (h_pions[tgt + "_z"]) h_pions[tgt + "_z"]->Write();
            if (h_pions[tgt + "_pt2"]) h_pions[tgt + "_pt2"]->Write();
            if (h_pions[tgt + "_vz"]) h_pions[tgt + "_vz"]->Write();
        }
        if (canvas) canvas->Write();
        f_out->Close();
        delete f_out;
    }

    // Clean up histograms
    for (const auto& pair : h_electrons) {
        if (pair.second) {
            gROOT->Remove(pair.second);
            delete pair.second;
        }
    }
    for (const auto& pair : h_pions) {
        if (pair.second) {
            gROOT->Remove(pair.second);
            delete pair.second;
        }
    }
    delete g_Cu;
    delete g_Sn;
    delete g_CxC;
    delete canvas;
    delete legend;

    return 0;
}