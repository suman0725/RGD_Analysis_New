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
    int nBins = 7;
    double nuMin = 1.5;
    double nuMax = 8.5;

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
                h_LD2_ele->Fill(nu);
            }

            // Positive pion selection
            if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pip_temp, p);
                foundPionP = true;
                h_LD2_piP->Fill(nu);
            }

            // Negative pion selection
            if (pid == -211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pim_temp, p);
                foundPionM = true;
                h_LD2_piM->Fill(nu);
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
                h_Carbon_ele->Fill(nu_Carbon);
            }

            if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pip_temp_Carbon, p);
                foundPionP_Carbon = true;
                h_Carbon_piP->Fill(nu_Carbon);
            }

            if (pid == -211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pim_temp_Carbon, p);
                foundPionM_Carbon = true;
                h_Carbon_piM->Fill(nu_Carbon);
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
                    h_Cu_ele->Fill(nu_CuSn);
                }

                if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                    SetLorentzVector(pip_temp_CuSn, p);
                    foundPionP_CuSn = true;
                    h_Cu_piP->Fill(nu_CuSn);
                }

                if (pid == -211 && chi2Pid > -10 && chi2Pid < 10) {
                    SetLorentzVector(pim_temp_CuSn, p);
                    foundPionM_CuSn = true;
                    h_Cu_piM->Fill(nu_CuSn);
                }
            } else if (Vz >= -6.137 && Vz <= 5) { // Sn target
                if (pid == 11 && status < 0 && chi2Pid < 5 && chi2Pid > -5) {
                    SetLorentzVector(el_temp_CuSn, p);
                    foundElectron_CuSn = true;
                    nu_CuSn = beam.Energy() - el_temp_CuSn.Energy();
                    q_CuSn = beam - el_temp_CuSn;
                    h_Sn_ele->Fill(nu_CuSn);
                }

                if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                    SetLorentzVector(pip_temp_CuSn, p);
                    foundPionP_CuSn = true;
                    h_Sn_piP->Fill(nu_CuSn);
                }

                if (pid == -211 && chi2Pid > -10 && chi2Pid < 10) {
                    SetLorentzVector(pim_temp_CuSn, p);
                    foundPionM_CuSn = true;
                    h_Sn_piM->Fill(nu_CuSn);
                }
            }
        }
    }









    // Process multiplicity ratio calculation and graphing for Cu and Sn in a similar way...
    // (Rest of your script follows without significant change)
    // Create TGraphErrors for Cu target
TGraphErrors* graph_multiplicity_Cu_piP = new TGraphErrors();
TGraphErrors* graph_multiplicity_Cu_piM = new TGraphErrors();
int pointIndex_Cu_piP = 0;
int pointIndex_Cu_piM = 0;
// Variables to track the maximum and minimum values
double maxY_Cu = -1e30;  // Very small initial value
double minY_Cu = 1e30;   // Very large initial value

for (int bin = 1; bin <= nBins; ++bin) {
    double binCenter = h_LD2_ele->GetBinCenter(bin);
    double count_LD2_ele = h_LD2_ele->GetBinContent(bin);
    double count_LD2_piP = h_LD2_piP->GetBinContent(bin);
    double count_LD2_piM = h_LD2_piM->GetBinContent(bin);
    double count_Cu_ele = h_Cu_ele->GetBinContent(bin);
    double count_Cu_piP = h_Cu_piP->GetBinContent(bin);
    double count_Cu_piM = h_Cu_piM->GetBinContent(bin);

    double R_h_A_piP = 0.0;
    double R_h_A_piM = 0.0;
    double error_R_h_A_piP = 0.0;
    double error_R_h_A_piM = 0.0;

    // Calculate multiplicity ratio for pi+
    if (count_LD2_ele > 0 && count_Cu_ele > 0) {
        if (count_LD2_piP > 0 && count_Cu_piP > 0) {
            R_h_A_piP = (count_Cu_piP / count_Cu_ele) / (count_LD2_piP / count_LD2_ele);
            error_R_h_A_piP = R_h_A_piP * sqrt(1.0 / count_Cu_piP + 1.0 / count_Cu_ele + 
                                               1.0 / count_LD2_piP + 1.0 / count_LD2_ele);
            graph_multiplicity_Cu_piP->SetPoint(pointIndex_Cu_piP, binCenter, R_h_A_piP);
            graph_multiplicity_Cu_piP->SetPointError(pointIndex_Cu_piP, 0, error_R_h_A_piP);
             // Update maxY and minY for Sn pi+
            if (R_h_A_piP + error_R_h_A_piP > maxY_Cu) maxY_Cu = R_h_A_piP + error_R_h_A_piP;
            if (R_h_A_piP - error_R_h_A_piP < minY_Cu) minY_Cu = R_h_A_piP - error_R_h_A_piP;
            pointIndex_Cu_piP++;
        }

        // Calculate multiplicity ratio for pi-
        if (count_LD2_piM > 0 && count_Cu_piM > 0) {
            R_h_A_piM = (count_Cu_piM / count_Cu_ele) / (count_LD2_piM / count_LD2_ele);
            error_R_h_A_piM = R_h_A_piM * sqrt(1.0 / count_Cu_piM + 1.0 / count_Cu_ele + 
                                               1.0 / count_LD2_piM + 1.0 / count_LD2_ele);
            graph_multiplicity_Cu_piM->SetPoint(pointIndex_Cu_piM, binCenter, R_h_A_piM);
            graph_multiplicity_Cu_piM->SetPointError(pointIndex_Cu_piM, 0, error_R_h_A_piM);

            // Update maxY and minY for Cu pi-
            if (R_h_A_piM + error_R_h_A_piM > maxY_Cu) maxY_Cu = R_h_A_piM + error_R_h_A_piM;
            if (R_h_A_piM - error_R_h_A_piM < minY_Cu) minY_Cu = R_h_A_piM - error_R_h_A_piM;
            pointIndex_Cu_piM++;
        }
    }
}

// Adjust y-axis range based on the calculated max and min values
double margin_Cu = 0.1 * (maxY_Cu - minY_Cu);  // Add 10% margin above and below the data
graph_multiplicity_Cu_piP->SetMinimum(minY_Cu - margin_Cu);    // Extend the y-axis minimum
graph_multiplicity_Cu_piP->SetMaximum(maxY_Cu + margin_Cu);    // Extend the y-axis maximum


// Create a canvas to draw both graphs for Cu
TCanvas* canvas_Cu = new TCanvas("canvas_multiplicity_Cu", "Multiplicity Ratio #pi^{+} and #pi^{-} vs #nu for Cu", 800, 600);
graph_multiplicity_Cu_piP->SetTitle("Multiplicity Ratio #pi^{+} and #pi^{-} vs #nu for Cu;nu;R_{M}^{h}");
graph_multiplicity_Cu_piP->SetMarkerStyle(kFullCircle);
graph_multiplicity_Cu_piP->SetMarkerColor(kBlue);  // Color for pi+
graph_multiplicity_Cu_piP->SetLineColor(kBlue);    // Set error bar color to blue
graph_multiplicity_Cu_piP->Draw("AP");  // Draw pi+ data first

graph_multiplicity_Cu_piM->SetMarkerStyle(kFullSquare);
graph_multiplicity_Cu_piM->SetMarkerColor(kRed);  // Color for pi-
graph_multiplicity_Cu_piM->SetLineColor(kRed);    // Set error bar color to red
graph_multiplicity_Cu_piM->Draw("P");  // Draw pi- data on the same canvas

// Create the legend
TLegend* legend_Cu = new TLegend(0.15, 0.9, 0.3, 0.8);  // Set position (x1, y1, x2, y2)
legend_Cu->AddEntry(graph_multiplicity_Cu_piP, "#pi^{+}", "p");
legend_Cu->AddEntry(graph_multiplicity_Cu_piM, "#pi^{-}", "p");
legend_Cu->SetTextSize(0.03);
legend_Cu->SetFillStyle(0);   
legend_Cu->SetBorderSize(0);

legend_Cu->Draw();  // Draw the legend

canvas_Cu->SaveAs("multiplicity_ratio_Cu.png");


// Create TGraphErrors for Sn target
TGraphErrors* graph_multiplicity_Sn_piP = new TGraphErrors();
TGraphErrors* graph_multiplicity_Sn_piM = new TGraphErrors();
int pointIndex_Sn_piP = 0;
int pointIndex_Sn_piM = 0;

// Variables to track the maximum and minimum values
double maxY_Sn = -1e30;  // Very small initial value
double minY_Sn = 1e30;   // Very large initial value

for (int bin = 1; bin <= nBins; ++bin) {
    double binCenter = h_LD2_ele->GetBinCenter(bin);
    double count_LD2_ele = h_LD2_ele->GetBinContent(bin);
    double count_LD2_piP = h_LD2_piP->GetBinContent(bin);
    double count_LD2_piM = h_LD2_piM->GetBinContent(bin);
    double count_Sn_ele = h_Sn_ele->GetBinContent(bin);
    double count_Sn_piP = h_Sn_piP->GetBinContent(bin);
    double count_Sn_piM = h_Sn_piM->GetBinContent(bin);

    double R_h_A_piP = 0.0;
    double R_h_A_piM = 0.0;
    double error_R_h_A_piP = 0.0;
    double error_R_h_A_piM = 0.0;

    // Calculate multiplicity ratio for pi+
    if (count_LD2_ele > 0 && count_Sn_ele > 0) {
        if (count_LD2_piP > 0 && count_Sn_piP > 0) {
            R_h_A_piP = (count_Sn_piP / count_Sn_ele) / (count_LD2_piP / count_LD2_ele);
            error_R_h_A_piP = R_h_A_piP * sqrt(1.0 / count_Sn_piP + 1.0 / count_Sn_ele + 
                                               1.0 / count_LD2_piP + 1.0 / count_LD2_ele);
            graph_multiplicity_Sn_piP->SetPoint(pointIndex_Sn_piP, binCenter, R_h_A_piP);
            graph_multiplicity_Sn_piP->SetPointError(pointIndex_Sn_piP, 0, error_R_h_A_piP);

            // Update maxY and minY for Sn pi+
            if (R_h_A_piP + error_R_h_A_piP > maxY_Sn) maxY_Sn = R_h_A_piP + error_R_h_A_piP;
            if (R_h_A_piP - error_R_h_A_piP < minY_Sn) minY_Sn = R_h_A_piP - error_R_h_A_piP;

            pointIndex_Sn_piP++;
        }

        // Calculate multiplicity ratio for pi-
        if (count_LD2_piM > 0 && count_Sn_piM > 0) {
            R_h_A_piM = (count_Sn_piM / count_Sn_ele) / (count_LD2_piM / count_LD2_ele);
            error_R_h_A_piM = R_h_A_piM * sqrt(1.0 / count_Sn_piM + 1.0 / count_Sn_ele + 
                                               1.0 / count_LD2_piM + 1.0 / count_LD2_ele);
            graph_multiplicity_Sn_piM->SetPoint(pointIndex_Sn_piM, binCenter, R_h_A_piM);
            graph_multiplicity_Sn_piM->SetPointError(pointIndex_Sn_piM, 0, error_R_h_A_piM);

            // Update maxY and minY for Sn pi-
            if (R_h_A_piM + error_R_h_A_piM > maxY_Sn) maxY_Sn = R_h_A_piM + error_R_h_A_piM;
            if (R_h_A_piM - error_R_h_A_piM < minY_Sn) minY_Sn = R_h_A_piM - error_R_h_A_piM;

            pointIndex_Sn_piM++;
        }
    }
}

// Adjust y-axis range based on the calculated max and min values
double margin_Sn = 0.1 * (maxY_Sn - minY_Sn);  // Add 10% margin above and below the data
graph_multiplicity_Sn_piP->SetMinimum(minY_Sn - margin_Sn);    // Extend the y-axis minimum
graph_multiplicity_Sn_piP->SetMaximum(maxY_Sn + margin_Sn);    // Extend the y-axis maximum

// Create a canvas to draw both graphs for Sn
TCanvas* canvas_Sn = new TCanvas("canvas_multiplicity_Sn", "Multiplicity Ratio #pi^{+} and #pi^{-} vs #nu for Sn", 800, 600);
graph_multiplicity_Sn_piP->SetTitle("Multiplicity Ratio #pi^{+} and #pi^{-} vs #nu for Sn;nu;R_{M}^{h}");
graph_multiplicity_Sn_piP->SetMarkerStyle(kFullCircle);
graph_multiplicity_Sn_piP->SetMarkerColor(kBlue);  // Color for pi+
graph_multiplicity_Sn_piP->SetLineColor(kBlue);    // Set error bar color to blue
graph_multiplicity_Sn_piP->Draw("AP");  // Draw pi+ data first

graph_multiplicity_Sn_piM->SetMarkerStyle(kFullSquare);
graph_multiplicity_Sn_piM->SetMarkerColor(kRed);  // Color for pi-
graph_multiplicity_Sn_piM->SetLineColor(kRed);    // Set error bar color to red
graph_multiplicity_Sn_piM->Draw("P");  // Draw pi- data on the same canvas

// Create the legend
TLegend* legend_Sn = new TLegend(0.15, 0.9, 0.3, 0.8);  // Set position (x1, y1, x2, y2)
legend_Sn->AddEntry(graph_multiplicity_Sn_piP, "#pi^{+}", "p");
legend_Sn->AddEntry(graph_multiplicity_Sn_piM, "#pi^{-}", "p");
legend_Sn->SetTextSize(0.03);
legend_Sn->SetFillStyle(0);   
legend_Sn->SetBorderSize(0);

legend_Sn->Draw();  // Draw the legend

canvas_Sn->SaveAs("multiplicity_ratio_Sn.png");



// Create TGraphErrors for Carbon target
TGraphErrors* graph_multiplicity_Carbon_piP = new TGraphErrors();
TGraphErrors* graph_multiplicity_Carbon_piM = new TGraphErrors();
int pointIndex_Carbon_piP = 0;
int pointIndex_Carbon_piM = 0;

// Variables to track the maximum and minimum values for Carbon pi+ and pi-
double maxY_Carbon = -1e30;  // Very small initial value
double minY_Carbon = 1e30;   // Very large initial value

// Loop over bins to calculate multiplicity ratios for Carbon
for (int bin = 1; bin <= nBins; ++bin) {
    double binCenter = h_LD2_ele->GetBinCenter(bin);
    double count_LD2_ele = h_LD2_ele->GetBinContent(bin);
    double count_LD2_piP = h_LD2_piP->GetBinContent(bin);
    double count_LD2_piM = h_LD2_piM->GetBinContent(bin);
    double count_Carbon_ele = h_Carbon_ele->GetBinContent(bin);
    double count_Carbon_piP = h_Carbon_piP->GetBinContent(bin);
    double count_Carbon_piM = h_Carbon_piM->GetBinContent(bin);
    
    // For π⁺ (positive pion)
    
    if (count_LD2_ele > 0 && count_Carbon_ele > 0) {
        double R_h_A_piP_Carbon = (count_Carbon_piP / count_Carbon_ele) / (count_LD2_piP / count_LD2_ele);
        double error_R_h_A_piP_Carbon = R_h_A_piP_Carbon * sqrt(1.0 / count_Carbon_piP + 1.0 / count_Carbon_ele +
                                                                1.0 / count_LD2_piP + 1.0 / count_LD2_ele);
        graph_multiplicity_Carbon_piP->SetPoint(pointIndex_Carbon_piP, binCenter, R_h_A_piP_Carbon);
        graph_multiplicity_Carbon_piP->SetPointError(pointIndex_Carbon_piP, 0, error_R_h_A_piP_Carbon);

        // Update max and min for pi+
        if (R_h_A_piP_Carbon + error_R_h_A_piP_Carbon > maxY_Carbon) maxY_Carbon = R_h_A_piP_Carbon + error_R_h_A_piP_Carbon;
        if (R_h_A_piP_Carbon - error_R_h_A_piP_Carbon < minY_Carbon) minY_Carbon = R_h_A_piP_Carbon - error_R_h_A_piP_Carbon;

        pointIndex_Carbon_piP++;
    }

    // For π⁻ (negative pion)
   
    if (count_LD2_ele > 0 && count_Carbon_ele > 0) {
        double R_h_A_piM_Carbon = (count_Carbon_piM / count_Carbon_ele) / (count_LD2_piM / count_LD2_ele);
        double error_R_h_A_piM_Carbon = R_h_A_piM_Carbon * sqrt(1.0 / count_Carbon_piM + 1.0 / count_Carbon_ele +
                                                                1.0 / count_LD2_piM + 1.0 / count_LD2_ele);
        graph_multiplicity_Carbon_piM->SetPoint(pointIndex_Carbon_piM, binCenter, R_h_A_piM_Carbon);
        graph_multiplicity_Carbon_piM->SetPointError(pointIndex_Carbon_piM, 0, error_R_h_A_piM_Carbon);

        // Update max and min for pi-
        if (R_h_A_piM_Carbon + error_R_h_A_piM_Carbon > maxY_Carbon) maxY_Carbon = R_h_A_piM_Carbon + error_R_h_A_piM_Carbon;
        if (R_h_A_piM_Carbon - error_R_h_A_piM_Carbon < minY_Carbon) minY_Carbon = R_h_A_piM_Carbon - error_R_h_A_piM_Carbon;

        pointIndex_Carbon_piM++;
    }
}

// Adjust y-axis range based on the calculated max and min values
double margin_Carbon = 0.1 * (maxY_Carbon - minY_Carbon);  // Add 10% margin above and below the data
graph_multiplicity_Carbon_piP->SetMinimum(minY_Carbon - margin_Carbon);    // Extend the y-axis minimum
graph_multiplicity_Carbon_piP->SetMaximum(maxY_Carbon + margin_Carbon);    // Extend the y-axis maximum

// Create the canvas and draw the graphs for Carbon target
auto* canvas_Carbon = new TCanvas("canvas_multiplicity_Carbon", "Multiplicity Ratio #pi^{+} and #pi^{-} vs #nu for Sn", 800, 600);
graph_multiplicity_Carbon_piP->SetTitle("Multiplicity Ratio for #pi^{+} and #pi^{-} in Carbon Vs #nu;#nu;R_{M}^{h}");
graph_multiplicity_Carbon_piP->SetMarkerStyle(20);
graph_multiplicity_Carbon_piP->SetMarkerColor(kBlue);
graph_multiplicity_Carbon_piP->SetLineColor(kBlue);
graph_multiplicity_Carbon_piP->Draw("AP");

graph_multiplicity_Carbon_piM->SetMarkerStyle(21);
graph_multiplicity_Carbon_piM->SetMarkerColor(kRed);
graph_multiplicity_Carbon_piM->SetLineColor(kRed);
graph_multiplicity_Carbon_piM->Draw("P SAME");

// Add legend
auto* legend_Carbon = new TLegend(0.15, 0.9, 0.3, 0.8);
legend_Carbon->AddEntry(graph_multiplicity_Carbon_piP, "#pi^{+}", "p");
legend_Carbon->AddEntry(graph_multiplicity_Carbon_piM, "#pi^{-}", "p");
legend_Carbon->SetTextSize(0.03);
legend_Carbon->SetFillStyle(0); 
legend_Carbon->SetBorderSize(0);
legend_Carbon->Draw();

// Save the plot
canvas_Carbon->SaveAs("MultiplicityRatio_Carbon_PiPlusMinus.png");




}



