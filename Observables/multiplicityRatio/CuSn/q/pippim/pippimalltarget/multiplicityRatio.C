#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "HipoChain.h"
#include <cmath>
#include <filesystem>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <vector>
#include <TLegend.h> 
#include <TStyle.h>          



namespace fs = std::filesystem;
using namespace clas12;

void SetLorentzVector(TLorentzVector &p4, clas12::region_part_ptr rp) {
    double mass = TDatabasePDG::Instance()->GetParticle(rp->par()->getPid())->Mass();
    p4.SetXYZM(rp->par()->getPx(), rp->par()->getPy(), rp->par()->getPz(), mass);
}




void addHipoFiles(clas12root::HipoChain& chain, const fs::path& baseDir) {
    std::vector<fs::path> files;
    for (const auto& entry : fs::recursive_directory_iterator(baseDir)) {
        if (entry.path().extension() == ".hipo") {
            files.push_back(entry.path());
        }
    }
    std::sort(files.begin(), files.end());
    for (const auto& file : files) {
        chain.Add(file.string());
        std::cout << "Added file: " << file << std::endl;
    }
}



void multiplicityRatio() {
    auto db = TDatabasePDG::Instance();
    double beamEnergy = 10.532; // GeV
    TLorentzVector beam(0, 0, beamEnergy, beamEnergy); // Beam along z-axis
    TLorentzVector target(0, 0, 0, db->GetParticle(2212)->Mass()); // Proton target
    TLorentzVector el, pip, pim;

    // Bin for variable nu 
    int nBins = 6;
    double nuMin = 1.5;
    double nuMax = 7.5;

    // Histograms for LD, Cu, Sn targets (positive and negative pions)
    auto* h_LD2_ele = new TH1F("h_LD2_ele", "LD Electrons;nu;Counts", nBins, nuMin, nuMax);
    auto* h_LD2_piP = new TH1F("h_LD2_piP", "LD Positive Pions;nu;Counts", nBins, nuMin, nuMax);
    auto* h_LD2_piM = new TH1F("h_LD2_piM", "LD Negative Pions;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Cu_ele = new TH1F("h_Cu_ele", "Copper Electrons;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Cu_piP = new TH1F("h_Cu_piP", "Copper Positive Pions;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Cu_piM = new TH1F("h_Cu_piM", "Copper Negative Pions;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Sn_ele = new TH1F("h_Sn_ele", "Tin Electrons;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Sn_piP = new TH1F("h_Sn_piP", "Tin Positive Pions;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Sn_piM = new TH1F("h_Sn_piM", "Tin Negative Pions;nu;Counts", nBins, nuMin, nuMax);
    // Add histograms for Carbon target (positive and negative pions)
    auto* h_Carbon_ele = new TH1F("h_Carbon_ele", "Carbon Electrons;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Carbon_piP = new TH1F("h_Carbon_piP", "Carbon Positive Pions;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Carbon_piM = new TH1F("h_Carbon_piM", "Carbon Negative Pions;nu;Counts", nBins, nuMin, nuMax);

    // Load Hipo files
    clas12root::HipoChain chainLD2, chainC, chainCuSn;
    chainLD2.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v6/LD2/skim_run_018307.hipo");
    chainCuSn.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v6/CuSn/skim_run_018347.hipo");
    chainC.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v6/CxC/skim_run_018339.hipo");
    chainLD2.db()->turnOffQADB();
    chainCuSn.db()->turnOffQADB();
    chainC.db()->turnOffQADB();

    auto& c12 = chainLD2.C12ref();
    while (chainLD2.Next()) {
        auto particles = c12->getDetParticles();
        TLorentzVector el_temp, pip_temp, pim_temp;

        bool foundElectron = false;
        bool foundPionP = false;
        bool foundPionM = false;

        double nu = 0;
        double Q2 = 0; 
        TLorentzVector q;

        for (auto& p : particles) {
            int pid = p->par()->getPid();
            int status = p->par()->getStatus();
            int chi2Pid = p->par()->getChi2Pid();

            // Electron selection
            if (pid == 11 && status < 0 && chi2Pid < 5 && chi2Pid > -5) {
                SetLorentzVector(el_temp, p);
                foundElectron = true;
                nu = beam.Energy() - el_temp.Energy();
                q = beam - el_temp;
                Q2 = -q.Mag2(); 
                h_LD2_ele->Fill(Q2);
            }

            // Positive pion selection
            if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pip_temp, p);
                foundPionP = true;
                h_LD2_piP->Fill(Q2);
            }

            // Negative pion selection
            if (pid == -211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pim_temp, p);
                foundPionM = true;
                h_LD2_piM->Fill(Q2);
            }
        }
    }

    // Process Carbon target (similar to how Cu and Sn are processed)
auto& c12_Carbon = chainC.C12ref();  // Assuming you have loaded the Carbon chain
while (chainC.Next()) {
    auto particles = c12_Carbon->getDetParticles();
    TLorentzVector el_temp_Carbon, pip_temp_Carbon, pim_temp_Carbon;

    bool foundElectron_Carbon = false;
    bool foundPionP_Carbon = false;
    bool foundPionM_Carbon = false;

    double nu_Carbon = 0;
    double Q2_Carbon = 0; 
    TLorentzVector q_Carbon;

    for (auto& p : particles) {
        int pid = p->par()->getPid();
        int status = p->par()->getStatus();
        int chi2Pid = p->par()->getChi2Pid();

        double Vz = p->par()->getVz();

        // Process only for Carbon target (assume Carbon target Vz range if necessary)
        
            if (pid == 11 && status < 0 && chi2Pid < 5 && chi2Pid > -5) {
                SetLorentzVector(el_temp_Carbon, p);
                foundElectron_Carbon = true;
                nu_Carbon = beam.Energy() - el_temp_Carbon.Energy();
                q_Carbon = beam - el_temp_Carbon;
                Q2_Carbon = -q_Carbon.Mag2(); 
                h_Carbon_ele->Fill(Q2_Carbon);
            }

            if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pip_temp_Carbon, p);
                foundPionP_Carbon = true;
                h_Carbon_piP->Fill(Q2_Carbon);
            }

            if (pid == -211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pim_temp_Carbon, p);
                foundPionM_Carbon = true;
                h_Carbon_piM->Fill(Q2_Carbon);
            }
    
    }
}


    // Process Cu and Sn targets
    auto& c12_CuSn = chainCuSn.C12ref();
    while (chainCuSn.Next()) {
        auto particles = c12_CuSn->getDetParticles();
        TLorentzVector el_temp_CuSn, pip_temp_CuSn, pim_temp_CuSn;

        bool foundElectron_CuSn = false;
        bool foundPionP_CuSn = false;
        bool foundPionM_CuSn = false;

        double nu_CuSn = 0;
        double Q2_Cu= 0; 
        double Q2_Sn= 0; 
        TLorentzVector q_CuSn;

        for (auto& p : particles) {
            int pid = p->par()->getPid();
            int status = p->par()->getStatus();
            int chi2Pid = p->par()->getChi2Pid();

            // Get Vz value for separation
            double Vz = p->par()->getVz();

            // Determine target based on Vz range
            if (Vz >= -11.463 && Vz <= -6.576) { // Cu target
                if (pid == 11 && status < 0 && chi2Pid < 5 && chi2Pid > -5) {
                    SetLorentzVector(el_temp_CuSn, p);
                    foundElectron_CuSn = true;
                    nu_CuSn = beam.Energy() - el_temp_CuSn.Energy();
                    q_CuSn = beam - el_temp_CuSn;
                    Q2_Cu = -q_CuSn.Mag2(); 
                    h_Cu_ele->Fill(Q2_Cu);
                }

                if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                    SetLorentzVector(pip_temp_CuSn, p);
                    foundPionP_CuSn = true;
                    h_Cu_piP->Fill(Q2_Cu);
                }

                if (pid == -211 && chi2Pid > -10 && chi2Pid < 10) {
                    SetLorentzVector(pim_temp_CuSn, p);
                    foundPionM_CuSn = true;
                    h_Cu_piM->Fill(Q2_Cu);
                }
            } else if (Vz >= -6.137 && Vz <= 5) { // Sn target
                if (pid == 11 && status < 0 && chi2Pid < 5 && chi2Pid > -5) {
                    SetLorentzVector(el_temp_CuSn, p);
                    foundElectron_CuSn = true;
                    nu_CuSn = beam.Energy() - el_temp_CuSn.Energy();
                    q_CuSn = beam - el_temp_CuSn;
                    Q2_Sn = -q_CuSn.Mag2();
                    h_Sn_ele->Fill(Q2_Sn);
                }

                if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                    SetLorentzVector(pip_temp_CuSn, p);
                    foundPionP_CuSn = true;
                    h_Sn_piP->Fill(Q2_Sn);
                }

                if (pid == -211 && chi2Pid > -10 && chi2Pid < 10) {
                    SetLorentzVector(pim_temp_CuSn, p);
                    foundPionM_CuSn = true;
                    h_Sn_piM->Fill(Q2_Sn);
                }
            }
        }
    }









    // Process multiplicity ratio calculation and graphing for Cu and Sn in a similar way...
    // (Rest of your script follows without significant change)
    // Create TGraphErrors for Cu target
TGraphErrors* graph_multiplicity_piP_Cu = new TGraphErrors();
TGraphErrors* graph_multiplicity_piP_Sn = new TGraphErrors();
TGraphErrors* graph_multiplicity_piP_Carbon = new TGraphErrors();

int pointIndex_piP_Cu = 0, pointIndex_piP_Sn = 0, pointIndex_piP_Carbon = 0;
double maxY_piP = -1e30;  // Initialize with a small value
double minY_piP = 1e30;   // Initialize with a large value
double R_h_A_piP_Cu = 0.0;
double error_R_h_A_piP_Cu = 0.0;
double R_h_A_piP_Sn = 0.0;
double error_R_h_A_piP_Sn = 0.0;
double R_h_A_piP_Carbon = 0.0;
double error_R_h_A_piP_Carbon = 0.0;

for (int bin = 1; bin <= nBins; ++bin) {
    double binCenter = h_LD2_ele->GetBinCenter(bin);
    double count_LD2_ele = h_LD2_ele->GetBinContent(bin);
    double count_LD2_piP = h_LD2_piP->GetBinContent(bin);

    // Get counts for all targets
    double count_Cu_ele = h_Cu_ele->GetBinContent(bin);
    double count_Cu_piP = h_Cu_piP->GetBinContent(bin);
    double count_Sn_ele = h_Sn_ele->GetBinContent(bin);
    double count_Sn_piP = h_Sn_piP->GetBinContent(bin);
    double count_Carbon_ele = h_Carbon_ele->GetBinContent(bin);
    double count_Carbon_piP = h_Carbon_piP->GetBinContent(bin);

    // Calculate multiplicity ratio R_h_A for Cu
    if (count_LD2_ele > 0 && count_Cu_ele > 0 && count_LD2_piP > 0 && count_Cu_piP > 0) {
        R_h_A_piP_Cu = (count_Cu_piP / count_Cu_ele) / (count_LD2_piP / count_LD2_ele);
        error_R_h_A_piP_Cu = R_h_A_piP_Cu * sqrt(1.0 / count_Cu_piP + 1.0 / count_Cu_ele + 1.0 / count_LD2_piP + 1.0 / count_LD2_ele);
        graph_multiplicity_piP_Cu->SetPoint(pointIndex_piP_Cu, binCenter, R_h_A_piP_Cu);
        graph_multiplicity_piP_Cu->SetPointError(pointIndex_piP_Cu, 0, error_R_h_A_piP_Cu);
        pointIndex_piP_Cu++;
    }

    // Calculate multiplicity ratio R_h_A for Sn
    if (count_LD2_ele > 0 && count_Sn_ele > 0 && count_LD2_piP > 0 && count_Sn_piP > 0) {
        R_h_A_piP_Sn = (count_Sn_piP / count_Sn_ele) / (count_LD2_piP / count_LD2_ele);
        error_R_h_A_piP_Sn = R_h_A_piP_Sn * sqrt(1.0 / count_Sn_piP + 1.0 / count_Sn_ele + 1.0 / count_LD2_piP + 1.0 / count_LD2_ele);
        graph_multiplicity_piP_Sn->SetPoint(pointIndex_piP_Sn, binCenter, R_h_A_piP_Sn);
        graph_multiplicity_piP_Sn->SetPointError(pointIndex_piP_Sn, 0, error_R_h_A_piP_Sn);
        pointIndex_piP_Sn++;
    }

    // Calculate multiplicity ratio R_h_A for Carbon
    if (count_LD2_ele > 0 && count_Carbon_ele > 0 && count_LD2_piP > 0 && count_Carbon_piP > 0) {
        R_h_A_piP_Carbon = (count_Carbon_piP / count_Carbon_ele) / (count_LD2_piP / count_LD2_ele);
        error_R_h_A_piP_Carbon = R_h_A_piP_Carbon * sqrt(1.0 / count_Carbon_piP + 1.0 / count_Carbon_ele + 1.0 / count_LD2_piP + 1.0 / count_LD2_ele);
        graph_multiplicity_piP_Carbon->SetPoint(pointIndex_piP_Carbon, binCenter, R_h_A_piP_Carbon);
        graph_multiplicity_piP_Carbon->SetPointError(pointIndex_piP_Carbon, 0, error_R_h_A_piP_Carbon);
        pointIndex_piP_Carbon++;
    }

    // Update the y-axis limits
    double R_h_A_max = std::max({R_h_A_piP_Cu + error_R_h_A_piP_Cu, R_h_A_piP_Sn + error_R_h_A_piP_Sn, R_h_A_piP_Carbon + error_R_h_A_piP_Carbon});
    double R_h_A_min = std::min({R_h_A_piP_Cu - error_R_h_A_piP_Cu, R_h_A_piP_Sn - error_R_h_A_piP_Sn, R_h_A_piP_Carbon - error_R_h_A_piP_Carbon});

    if (R_h_A_max > maxY_piP) maxY_piP = R_h_A_max;
    if (R_h_A_min < minY_piP) minY_piP = R_h_A_min;
}

double margin_piP = 0.1 * (maxY_piP - minY_piP);  // 10% margin
graph_multiplicity_piP_Cu->SetMinimum(minY_piP - margin_piP);
graph_multiplicity_piP_Cu->SetMaximum(maxY_piP + margin_piP);

// Create a canvas for the pi+ graph
TCanvas* canvas_piP = new TCanvas("canvas_multiplicity_piP", "Multiplicity Ratio vs Q^{2} for Carbon, Copper, and Tin", 800, 600);
graph_multiplicity_piP_Cu->SetTitle("Multiplicity Ratio vs Q^{2} for Carbon, Copper, and Tin");
graph_multiplicity_piP_Cu->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
graph_multiplicity_piP_Cu->GetYaxis()->SetTitle("R_{M}^{#pi^{+}}");

// Draw the graph for Cu
graph_multiplicity_piP_Cu->SetMarkerStyle(20);
graph_multiplicity_piP_Cu->SetMarkerColor(kRed);
graph_multiplicity_piP_Cu->SetLineColor(kRed);
graph_multiplicity_piP_Cu->Draw("AP");

// Draw the graph for Sn on top of Cu
graph_multiplicity_piP_Sn->SetMarkerStyle(21);
graph_multiplicity_piP_Sn->SetMarkerColor(kBlue);
graph_multiplicity_piP_Sn->SetLineColor(kBlue);
graph_multiplicity_piP_Sn->Draw("P same");

// Draw the graph for Carbon on top of Sn and Cu
graph_multiplicity_piP_Carbon->SetMarkerStyle(22);
graph_multiplicity_piP_Carbon->SetMarkerColor(kGreen);
graph_multiplicity_piP_Carbon->SetLineColor(kGreen);
graph_multiplicity_piP_Carbon->Draw("P same");

// Add a legend to distinguish between the targets
TLegend* legend = new TLegend(0.11, 0.9, 0.3, 0.8);
legend->AddEntry(graph_multiplicity_piP_Carbon, "C", "p");
legend->AddEntry(graph_multiplicity_piP_Cu, "Cu", "p");
legend->AddEntry(graph_multiplicity_piP_Sn, "Sn", "p");
legend->SetTextSize(0.03);
legend->SetFillStyle(0); 
legend->SetBorderSize(0);
legend->Draw();


// Save the plot
canvas_piP->SaveAs("MultiplicityRatio_Cu_Sn_Carbon_PiPlus.png");


TGraphErrors* graph_multiplicity_piM_Cu = new TGraphErrors();
TGraphErrors* graph_multiplicity_piM_Sn = new TGraphErrors();
TGraphErrors* graph_multiplicity_piM_Carbon = new TGraphErrors();

int pointIndex_piM_Cu = 0, pointIndex_piM_Sn = 0, pointIndex_piM_Carbon = 0;
double maxY_piM = -1e30;  // Initialize with a small value
double minY_piM = 1e30;   // Initialize with a large value
double R_h_A_piM_Cu = 0.0;
double error_R_h_A_piM_Cu = 0.0;
double R_h_A_piM_Sn = 0.0;
double error_R_h_A_piM_Sn = 0.0;
double R_h_A_piM_Carbon = 0.0;
double error_R_h_A_piM_Carbon = 0.0;

for (int bin = 1; bin <= nBins; ++bin) {
    double binCenter = h_LD2_ele->GetBinCenter(bin);
    double count_LD2_ele = h_LD2_ele->GetBinContent(bin);
    double count_LD2_piM = h_LD2_piM->GetBinContent(bin);

    // Get counts for all targets
    double count_Cu_ele = h_Cu_ele->GetBinContent(bin);
    double count_Cu_piM = h_Cu_piM->GetBinContent(bin);
    double count_Sn_ele = h_Sn_ele->GetBinContent(bin);
    double count_Sn_piM = h_Sn_piM->GetBinContent(bin);
    double count_Carbon_ele = h_Carbon_ele->GetBinContent(bin);
    double count_Carbon_piM = h_Carbon_piM->GetBinContent(bin);

    // Calculate multiplicity ratio R_h_A for Cu
    if (count_LD2_ele > 0 && count_Cu_ele > 0 && count_LD2_piM > 0 && count_Cu_piM > 0) {
        R_h_A_piM_Cu = (count_Cu_piM / count_Cu_ele) / (count_LD2_piM / count_LD2_ele);
        error_R_h_A_piM_Cu = R_h_A_piM_Cu * sqrt(1.0 / count_Cu_piM + 1.0 / count_Cu_ele + 1.0 / count_LD2_piM + 1.0 / count_LD2_ele);
        graph_multiplicity_piM_Cu->SetPoint(pointIndex_piM_Cu, binCenter, R_h_A_piM_Cu);
        graph_multiplicity_piM_Cu->SetPointError(pointIndex_piM_Cu, 0, error_R_h_A_piM_Cu);
        pointIndex_piM_Cu++;
    }

    // Calculate multiplicity ratio R_h_A for Sn
    if (count_LD2_ele > 0 && count_Sn_ele > 0 && count_LD2_piM > 0 && count_Sn_piM > 0) {
        R_h_A_piM_Sn = (count_Sn_piM / count_Sn_ele) / (count_LD2_piM / count_LD2_ele);
        error_R_h_A_piM_Sn = R_h_A_piM_Sn * sqrt(1.0 / count_Sn_piM + 1.0 / count_Sn_ele + 1.0 / count_LD2_piM + 1.0 / count_LD2_ele);
        graph_multiplicity_piM_Sn->SetPoint(pointIndex_piM_Sn, binCenter, R_h_A_piM_Sn);
        graph_multiplicity_piM_Sn->SetPointError(pointIndex_piM_Sn, 0, error_R_h_A_piM_Sn);
        pointIndex_piM_Sn++;
    }

    // Calculate multiplicity ratio R_h_A for Carbon
    if (count_LD2_ele > 0 && count_Carbon_ele > 0 && count_LD2_piM > 0 && count_Carbon_piM > 0) {
        R_h_A_piM_Carbon = (count_Carbon_piM / count_Carbon_ele) / (count_LD2_piM / count_LD2_ele);
        error_R_h_A_piM_Carbon = R_h_A_piM_Carbon * sqrt(1.0 / count_Carbon_piM + 1.0 / count_Carbon_ele + 1.0 / count_LD2_piM + 1.0 / count_LD2_ele);
        graph_multiplicity_piM_Carbon->SetPoint(pointIndex_piM_Carbon, binCenter, R_h_A_piM_Carbon);
        graph_multiplicity_piM_Carbon->SetPointError(pointIndex_piM_Carbon, 0, error_R_h_A_piM_Carbon);
        pointIndex_piM_Carbon++;
    }

    // Update the y-axis limits
    double R_h_A_max = std::max({R_h_A_piM_Cu + error_R_h_A_piM_Cu, R_h_A_piM_Sn + error_R_h_A_piM_Sn, R_h_A_piM_Carbon + error_R_h_A_piM_Carbon});
    double R_h_A_min = std::min({R_h_A_piM_Cu - error_R_h_A_piM_Cu, R_h_A_piM_Sn - error_R_h_A_piM_Sn, R_h_A_piM_Carbon - error_R_h_A_piM_Carbon});

    if (R_h_A_max > maxY_piM) maxY_piM = R_h_A_max;
    if (R_h_A_min < minY_piM) minY_piM = R_h_A_min;
}

double margin_piM = 0.1 * (maxY_piM - minY_piM);  // 10% margin
graph_multiplicity_piM_Cu->SetMinimum(minY_piM - margin_piM);
graph_multiplicity_piM_Cu->SetMaximum(maxY_piM + margin_piM);

// Create a canvas for the pi- graph
TCanvas* canvas_piM = new TCanvas("canvas_multiplicity_piM", "Multiplicity Ratio vs Q^{2} for Carbon, Copper, and Tin", 800, 600);
graph_multiplicity_piM_Cu->SetTitle("Multiplicity Ratio vs Q^{2} for Carbon, Copper, and Tin");
graph_multiplicity_piM_Cu->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
graph_multiplicity_piM_Cu->GetYaxis()->SetTitle("R_{M}^{#pi^{-}}");

// Draw the graph for Cu
graph_multiplicity_piM_Cu->SetMarkerStyle(20);
graph_multiplicity_piM_Cu->SetMarkerColor(kRed);
graph_multiplicity_piM_Cu->SetLineColor(kRed);
graph_multiplicity_piM_Cu->Draw("AP");

// Draw the graph for Sn on top of Cu
graph_multiplicity_piM_Sn->SetMarkerStyle(21);
graph_multiplicity_piM_Sn->SetMarkerColor(kBlue);
graph_multiplicity_piM_Sn->SetLineColor(kBlue);
graph_multiplicity_piM_Sn->Draw("P same");

// Draw the graph for Carbon on top of Sn and Cu
graph_multiplicity_piM_Carbon->SetMarkerStyle(22);
graph_multiplicity_piM_Carbon->SetMarkerColor(kGreen);
graph_multiplicity_piM_Carbon->SetLineColor(kGreen);
graph_multiplicity_piM_Carbon->Draw("P same");

// Add a legend to distinguish between the targets
TLegend* legend1 = new TLegend(0.11, 0.9, 0.3, 0.8);
legend1->AddEntry(graph_multiplicity_piM_Carbon, "C", "p");
legend1->AddEntry(graph_multiplicity_piM_Cu, "Cu", "p");
legend1->AddEntry(graph_multiplicity_piM_Sn, "Sn", "p");
legend1->SetTextSize(0.03);
legend1->SetFillStyle(0); 
legend1->SetBorderSize(0);
legend1->Draw();

// Save the plot
canvas_piM->SaveAs("MultiplicityRatio_Cu_Sn_Carbon_PiMinus.png");


}



