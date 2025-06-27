#include <cstdlib>
#include <iostream>
#include <filesystem>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <TFile.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include "reader.h"

namespace fs = std::filesystem;
using namespace std;

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
    // Check command-line arguments
    if (argc != 6) {
        std::cerr << "Usage: ./main <target1> <target2> <target3> <config> <pion_type>\n";
        std::cerr << "Example: ./main CxC LD2 CuSn OB PP\n";
        return 1;
    }

    std::string target1 = argv[1]; // CxC
    std::string target2 = argv[2]; // LD2
    std::string target3 = argv[3]; // CuSn
    std::string config = argv[4];  // OB or IB
    std::string pion_type = argv[5]; // PP or NP

    // Validate inputs
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

    // Hardcoded nu range
    const double nu_min = 4.45;
    const double nu_max = 6.05;

    // p_t^2 binning
    std::vector<double> pt2_bins = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15};
    int n_pt2_bins = pt2_bins.size() - 1;

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

    // 2D histograms for pions (z vs. p_t^2) and electron counters
    std::map<std::string, TH2F*> h_pi_2d;
    std::map<std::string, double> total_ele;
    for (const auto& tgt : {"LD2", "Cu", "Sn", "CxC"}) {
        h_pi_2d[tgt] = new TH2F((std::string("h_pi_") + tgt).c_str(), (tgt + std::string(" Pions; z; p_{t}^2 (GeV^2)")).c_str(),
                                10, 0.3, 0.7, n_pt2_bins, &pt2_bins[0]);
        total_ele[tgt] = 0; // Initialize electron counter
    }

    for (const auto& tgt : targets) {
        if (tgt != "CuSn" && tgt != "LD2" && tgt != "CxC") {
            std::cerr << "Invalid target: " << tgt << "\n";
            continue;
        }

        std::string targetDir = "/lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v11/" + tgt;
        std::vector<std::string> hipoFiles = loadHipoFiles(targetDir);
        if (hipoFiles.empty()) {
            std::cerr << "No HIPO files found in " << targetDir << "\n";
            continue;
        }

        std::cout << "Processing target: " << tgt << "\n";
        int file_count = 0;

        for (const auto& file : hipoFiles) {
            file_count++;
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
                    continue;
                }
                reader.open(file.c_str()); // Reset to process all events
            }

            int event_count = 0;
            while (reader.next() && event_count < 10000) {
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

                // Electron selection
                for (int i = 0; i < particles.size(); ++i) {
                    const auto& p = particles[i];
                    if (p.pid != 11) continue;

                    bool isCu = false, isSn = false;
                    if (tgt == "CuSn") {
                        isCu = (p.vz >= vzCutsEle["Cu"].first && p.vz <= vzCutsEle["Cu"].second);
                        isSn = (p.vz >= vzCutsEle["Sn"].first && p.vz <= vzCutsEle["Sn"].second);
                        if (!isCu && !isSn) continue;
                    } else {
                        if (p.vz < vzCutsEle[tgt].first || p.vz > vzCutsEle[tgt].second) continue;
                    }

                    if (p.status >= 0 || i != 0 || p.chi2pid < -5 || p.chi2pid > 5) continue;
                    if (abs(p.status) / 1000 != 2) continue;

                    el_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(11)->Mass());
                    nu = beamEnergy - el_temp.E();
                    if (nu <= nu_min || nu >= nu_max) continue;
                    TLorentzVector q = beam - el_temp;
                    Q2 = -q.Mag2();
                    y = nu / beamEnergy;
                    W = sqrt(protonMass * protonMass + 2 * protonMass * nu - Q2);

                    if (Q2 <= 1 || y <= 0.25 || y >= 0.85 || W <= 2) continue;
                    hasTrigger = true;

                    if (tgt == "CuSn") {
                        if (isCu) total_ele["Cu"]++;
                        if (isSn) total_ele["Sn"]++;
                    } else {
                        total_ele[tgt]++;
                    }
                    break;
                }

                if (!hasTrigger) continue;

                // Pion selection
                for (const auto& p : particles) {
                    if (p.pid != pion_pid) continue;
                    if (abs(p.status) / 1000 != 2 || p.chi2pid < -5 || p.chi2pid > 5) continue;

                    TLorentzVector pi_temp;
                    pi_temp.SetXYZM(p.px, p.py, p.pz, db->GetParticle(pion_pid)->Mass());
                    double z = pi_temp.E() / nu;
                    double pt2 = calculatePt2(pi_temp, beam - el_temp);
                    if (z <= 0.3 || z >= 0.7 || pt2 >= 1.2) continue;

                    if (tgt == "CuSn") {
                        bool isCu = (p.vz >= vzCutsPi["Cu"].first && p.vz <= vzCutsPi["Cu"].second);
                        bool isSn = (p.vz >= vzCutsPi["Sn"].first && p.vz <= vzCutsPi["Sn"].second);
                        if (isCu && pt2 >= 0.05 && pt2 < 1.15) h_pi_2d["Cu"]->Fill(z, pt2);
                        if (isSn && pt2 >= 0.05 && pt2 < 1.15) h_pi_2d["Sn"]->Fill(z, pt2);
                    } else {
                        if (p.vz >= vzCutsPi[tgt].first && p.vz <= vzCutsPi[tgt].second && pt2 >= 0.05 && pt2 < 1.15) {
                            h_pi_2d[tgt]->Fill(z, pt2);
                        }
                    }
                }
            }
            if (file_count == 1) break; // Only process one file per target for testing
            std::cout << "Processed " << event_count << " events in file: " << file << "\n";
        }
        
    }

    // Calculate and plot multiplicity ratio R
    TCanvas* canvas = new TCanvas("canvas", "Multiplicity Ratio", 800, 600);
    std::string pdf_output = "MultiplicityRatio_" + config + "_" + pion_type + "_z_nu4.45-6.05.pdf";

    for (int i = 0; i < n_pt2_bins; ++i) {
        double pt2_center = (pt2_bins[i] + pt2_bins[i + 1]) / 2;
        TGraphErrors* g_Cu = new TGraphErrors();
        TGraphErrors* g_Sn = new TGraphErrors();
        TGraphErrors* g_CxC = new TGraphErrors();
        int pointCu = 0, pointSn = 0, pointCxC = 0;
        double maxY = -1e30, minY = 1e30;

        // Project 2D histogram onto z for this pt2 bin
        TH1D* h_Cu_z = h_pi_2d["Cu"]->ProjectionX("h_Cu_z", i + 1, i + 1);
        TH1D* h_Sn_z = h_pi_2d["Sn"]->ProjectionX("h_Sn_z", i + 1, i + 1);
        TH1D* h_CxC_z = h_pi_2d["CxC"]->ProjectionX("h_CxC_z", i + 1, i + 1);
        TH1D* h_LD2_z = h_pi_2d["LD2"]->ProjectionX("h_LD2_z", i + 1, i + 1);

        for (int bin = 1; bin <= 10; ++bin) {
            double z = h_LD2_z->GetBinCenter(bin);
            double ld2_pi = h_LD2_z->GetBinContent(bin);
            double cu_pi = h_Cu_z->GetBinContent(bin);
            double sn_pi = h_Sn_z->GetBinContent(bin);
            double cxc_pi = h_CxC_z->GetBinContent(bin);
            double ld2_ele = total_ele["LD2"];
            double cu_ele = total_ele["Cu"];
            double sn_ele = total_ele["Sn"];
            double cxc_ele = total_ele["CxC"];

            // Cu ratio
            if (ld2_pi > 0 && cu_pi > 0 && ld2_ele > 0 && cu_ele > 0) {
                double R = (cu_pi / cu_ele) / (ld2_pi / ld2_ele);
                double err = R * sqrt(1.0 / cu_pi + 1.0 / cu_ele + 1.0 / ld2_pi + 1.0 / ld2_ele);
                g_Cu->SetPoint(pointCu, z, R);
                g_Cu->SetPointError(pointCu, 0, err);
                pointCu++;
                maxY = std::max(maxY, R + err);
                minY = std::min(minY, R - err);
            }

            // Sn ratio
            if (ld2_pi > 0 && sn_pi > 0 && ld2_ele > 0 && sn_ele > 0) {
                double R = (sn_pi / sn_ele) / (ld2_pi / ld2_ele);
                double err = R * sqrt(1.0 / sn_pi + 1.0 / sn_ele + 1.0 / ld2_pi + 1.0 / ld2_ele);
                g_Sn->SetPoint(pointSn, z, R);
                g_Sn->SetPointError(pointSn, 0, err);
                pointSn++;
                maxY = std::max(maxY, R + err);
                minY = std::min(minY, R - err);
            }

            // CxC ratio
            if (ld2_pi > 0 && cxc_pi > 0 && ld2_ele > 0 && cxc_ele > 0) {
                double R = (cxc_pi / cxc_ele) / (ld2_pi / ld2_ele);
                double err = R * sqrt(1.0 / cxc_pi + 1.0 / cxc_ele + 1.0 / ld2_pi + 1.0 / ld2_ele);
                g_CxC->SetPoint(pointCxC, z, R);
                g_CxC->SetPointError(pointCxC, 0, err);
                pointCxC++;
                maxY = std::max(maxY, R + err);
                minY = std::min(minY, R - err);
            }
        }

        // Plot setup
        canvas->Clear();
        g_Cu->SetTitle(("Multiplicity Ratio vs z (p_{t}^2 = " + std::to_string(pt2_center) + " GeV^2, 4.45 < #nu < 6.05 GeV)").c_str());
        g_Cu->GetXaxis()->SetTitle("z");
        g_Cu->GetYaxis()->SetTitle(pion_type == "PP" ? "R_{M}^{#pi^{+}}" : "R_{M}^{#pi^{-}}");
        double margin = 0.1 * (maxY - minY);
        g_Cu->SetMinimum(std::max(0.0, minY - margin));
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

        // Save to PDF
        if (i == 0) {
            canvas->Print((pdf_output + "(").c_str()); // Open PDF
        } else if (i == n_pt2_bins - 1) {
            canvas->Print((pdf_output + ")").c_str()); // Close PDF
        } else {
            canvas->Print(pdf_output.c_str()); // Add page
        }
    }

    return 0;
}