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

namespace fs = std::filesystem;
using namespace clas12;

void SetLorentzVector(TLorentzVector &p4, clas12::region_part_ptr rp) {
    double mass = TDatabasePDG::Instance()->GetParticle(rp->par()->getPid())->Mass();
    p4.SetXYZM(rp->par()->getPx(), rp->par()->getPy(), rp->par()->getPz(), mass);
}

double CalculatePt2(TLorentzVector P_h, TLorentzVector q) {

    TVector3 q3 = q.Vect(); 
    TVector3 p_h3 = P_h.Vect(); 

    TVector3 cross_qph = q3.Cross(p_h3);

    // Calculate the magnitudes squared
    double magnitudeCrossProductSquared = cross_qph.Mag2(); // |q × p_h|^2
    double magnitudeQSquared = q3.Mag2();                     // |q|^2

    // Calculate Pt^2 = |q × p_h|^2 / |q|^2
    double Pt2 = magnitudeCrossProductSquared / magnitudeQSquared;

    return Pt2;


}

void multiplicityRatio() {
    auto db = TDatabasePDG::Instance();
    double beamEnergy = 10.532; // GeV
    TLorentzVector beam(0, 0, beamEnergy, beamEnergy); // Assuming the beam is along the z-axis
    TLorentzVector target(0, 0, 0, db->GetParticle(2212)->Mass()); // Proton target
    TLorentzVector el, pip, pim;

    int nBins = 10;
    double pt2Min = 0;
    double pt2Max = 2;
    double LD2_ele_nor = 0; 
    double C_ele_nor = 0;

    // Histograms for π⁺ (positive pion)
    auto* h_LD2_piP = new TH1F("h_LD2_piP", "LD Positive Pions;pt2;Counts", nBins, pt2Min, pt2Max);
    auto* h_C_piP = new TH1F("h_C_piP", "Carbon Positive Pions;pt2;Counts", nBins, pt2Min, pt2Max);

    // Histograms for π⁻ (negative pion)
    auto* h_LD2_piM = new TH1F("h_LD2_piM", "LD Negative Pions;pt2;Counts", nBins, pt2Min, pt2Max);
    auto* h_C_piM = new TH1F("h_C_piM", "Carbon Negative Pions;pt2;Counts", nBins, pt2Min, pt2Max);

    clas12root::HipoChain chainLD2, chainC;
    chainLD2.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstLD2/dst/recon/018309/rec_clas_018309.evio.00035-00039.hipo");
    chainC.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstCxC/dst/recon/018339/rec_clas_018339.evio.00015-00019.hipo");
    chainLD2.db()->turnOffQADB();
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
                nu = beam.Energy() - el_temp.Energy(); // Calculate nu
                q = beam - el_temp; // Calculate q after finding the electron
                Q2 = -q.Mag2();
                LD2_ele_nor++;
            }

            // Positive pion (π⁺) selection
            if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pip_temp, p);
                foundPionP = true;
            }

            // Negative pion (π⁻) selection
            if (pid == -211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pim_temp, p);
                foundPionM = true;
            }
        }

        // Fill histograms if both electron and pion (π⁺ or π⁻) are found
        if (foundElectron && foundPionP) {
            double Pt2 = CalculatePt2(pip_temp, q); 
            h_LD2_piP->Fill(Pt2); 

            
        }
        if (foundElectron && foundPionM) {
            double Pt2 = CalculatePt2(pim_temp, q); 
            h_LD2_piM->Fill(Pt2); 
        }
    }

    // Process Carbon target
    auto& c12_C = chainC.C12ref();
    while (chainC.Next()) {
        auto particles = c12_C->getDetParticles();
        TLorentzVector el_temp_C, pip_temp_C, pim_temp_C;

        bool foundElectron_C = false;
        bool foundPionP_C = false;
        bool foundPionM_C = false;

        double nu_C = 0;
        double Q2_C = 0;
        TLorentzVector q_C;

        for (auto& p : particles) {
            int pid = p->par()->getPid();
            int status = p->par()->getStatus();
            int chi2Pid = p->par()->getChi2Pid();

            if (pid == 11 && status < 0 && chi2Pid < 5 && chi2Pid > -5) {
                SetLorentzVector(el_temp_C, p);
                foundElectron_C = true;
                nu_C = beam.Energy() - el_temp_C.Energy();
                q_C = beam - el_temp_C;
                Q2_C = -q_C.Mag2();
                C_ele_nor++;
            }

            if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pip_temp_C, p);
                foundPionP_C = true;
            }

            if (pid == -211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pim_temp_C, p);
                foundPionM_C = true;
            }
        }

        if (foundElectron_C && foundPionP_C) {
            double Pt2_C = CalculatePt2(pim_temp_C, q_C); 
            h_C_piP->Fill(Pt2_C); ;
        }
        if (foundElectron_C && foundPionM_C) {
            double Pt2_C = CalculatePt2(pim_temp_C, q_C); 
            h_C_piM->Fill(Pt2_C); 
        }
    }

   
    // Create TGraphErrors for π⁺ and π⁻
TGraphErrors* graphP = new TGraphErrors();
TGraphErrors* graphM = new TGraphErrors();
int pointIndexP = 0, pointIndexM = 0;

// Variables to track the maximum and minimum values
double maxY = -1e30;  // Very small initial value
double minY = 1e30;   // Very large initial value

for (int bin = 1; bin <= nBins; ++bin) {
    double binCenter = h_LD2_piP->GetBinCenter(bin);
    
    // For π⁺ (positive pion)
    double count_LD2_piP = h_LD2_piP->GetBinContent(bin);
    double count_C_piP = h_C_piP->GetBinContent(bin);
    if (LD2_ele_nor > 0 && C_ele_nor > 0) {
        double R_h_A_P = (count_C_piP / C_ele_nor) / (count_LD2_piP / LD2_ele_nor);
        double error_R_h_A_P = R_h_A_P * sqrt(1.0 / count_C_piP + 1.0 / C_ele_nor + 
                                              1.0 / count_LD2_piP + 1.0 / LD2_ele_nor);
        graphP->SetPoint(pointIndexP, binCenter, R_h_A_P);
        graphP->SetPointError(pointIndexP, 0, error_R_h_A_P);

        // Update maxY and minY
        if (R_h_A_P + error_R_h_A_P > maxY) maxY = R_h_A_P + error_R_h_A_P;
        if (R_h_A_P - error_R_h_A_P < minY) minY = R_h_A_P - error_R_h_A_P;

        pointIndexP++;
    }

    // For π⁻ (negative pion)
    double count_LD2_piM = h_LD2_piM->GetBinContent(bin);
    double count_C_piM = h_C_piM->GetBinContent(bin);
    if (LD2_ele_nor > 0 && C_ele_nor > 0) {
        double R_h_A_M = (count_C_piM / C_ele_nor) / (count_LD2_piM / LD2_ele_nor);
        double error_R_h_A_M = R_h_A_M * sqrt(1.0 / count_C_piM + 1.0 / C_ele_nor + 
                                              1.0 / count_LD2_piM + 1.0 / LD2_ele_nor);
        graphM->SetPoint(pointIndexM, binCenter, R_h_A_M);
        graphM->SetPointError(pointIndexM, 0, error_R_h_A_M);

        // Update maxY and minY
        if (R_h_A_M + error_R_h_A_M > maxY) maxY = R_h_A_M + error_R_h_A_M;
        if (R_h_A_M - error_R_h_A_M < minY) minY = R_h_A_M - error_R_h_A_M;

        pointIndexM++;
    }
}

// Adjust y-axis range based on the calculated max and min values
double margin = 0.1 * (maxY - minY);  // Add 10% margin above and below the data
graphP->SetMinimum(minY - margin);    // Extend the y-axis minimum
graphP->SetMaximum(maxY + margin);    // Extend the y-axis maximum

// Create the canvas and draw the graphs
auto* canvas = new TCanvas("canvas", "Multiplicity Ratio", 800, 600);
graphP->SetTitle("Multiplicity Ratio for #pi^{+} and #pi^{-} Vs P_{t}^{2};P_{t}^{2};R_{h}^{A}");
graphP->SetMarkerStyle(20);
graphP->SetMarkerColor(kBlue);
graphP->SetLineColor(kBlue);
graphP->Draw("AP");

graphM->SetMarkerStyle(21);
graphM->SetMarkerColor(kRed);
graphM->SetLineColor(kRed);
graphM->Draw("P SAME");

    auto* legend = new TLegend(0.1, 0.7, 0.9, 0.9);
    legend->AddEntry(graphP, "#pi^{+}", "p");
    legend->AddEntry(graphM, "#pi^{-}", "p");
    legend->SetTextSize(0.03);
    legend->SetFillStyle(0); 
    legend->SetBorderSize(0);
    legend->Draw();

    canvas->SaveAs("MultiplicityRatio_PiPlusMinus.png");
}


