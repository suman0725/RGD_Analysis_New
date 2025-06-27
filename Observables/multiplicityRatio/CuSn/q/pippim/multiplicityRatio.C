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
    double Q2Min = 1.0;
    double Q2Max = 8.0;

    // Histograms for LD and Nuclear targets (positive and negative pions)
    auto* h_LD2_ele = new TH1F("h_LD2_ele", "LD Electrons;nu;Counts", nBins, Q2Min, Q2Max);
    auto* h_LD2_piP = new TH1F("h_LD2_piP", "LD Positive Pions;nu;Counts", nBins, Q2Min, Q2Max);
    auto* h_LD2_piM = new TH1F("h_LD2_piM", "LD Negative Pions;nu;Counts", nBins, Q2Min, Q2Max);
    auto* h_C_ele = new TH1F("h_C_ele", "Carbon Electrons;nu;Counts", nBins, Q2Min, Q2Max);
    auto* h_C_piP = new TH1F("h_C_piP", "Carbon Positive Pions;nu;Counts", nBins, Q2Min, Q2Max);
    auto* h_C_piM = new TH1F("h_C_piM", "Carbon Negative Pions;nu;Counts", nBins,Q2Min, Q2Max);

    // Load Hipo files
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
                h_C_ele->Fill(Q2_C);
            }

            if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pip_temp_C, p);
                foundPionP_C = true;
                h_C_piP->Fill(Q2_C);
            }

            if (pid == -211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pim_temp_C, p);
                foundPionM_C = true;
                h_C_piM->Fill(Q2_C);
            }
        }
    }

TGraphErrors* graph_multiplicity_piP = new TGraphErrors();
TGraphErrors* graph_multiplicity_piM = new TGraphErrors();
int pointIndex_piP = 0;
int pointIndex_piM = 0;

for (int bin = 1; bin <= nBins; ++bin) {
    double binCenter = h_LD2_ele->GetBinCenter(bin);
    double count_LD2_ele = h_LD2_ele->GetBinContent(bin);
    double count_LD2_piP = h_LD2_piP->GetBinContent(bin);
    double count_LD2_piM = h_LD2_piM->GetBinContent(bin);
    double count_C_ele = h_C_ele->GetBinContent(bin);
    double count_C_piP = h_C_piP->GetBinContent(bin);
    double count_C_piM = h_C_piM->GetBinContent(bin);

    double R_h_A_piP = 0.0;
    double R_h_A_piM = 0.0;
    double error_R_h_A_piP = 0.0;
    double error_R_h_A_piM = 0.0;

    // Calculate multiplicity ratio for pi+
    if (count_LD2_ele > 0 && count_C_ele > 0) {
        if (count_LD2_piP > 0 && count_C_piP > 0) {
            R_h_A_piP = (count_C_piP / count_C_ele) / (count_LD2_piP / count_LD2_ele);
            error_R_h_A_piP = R_h_A_piP * sqrt(1.0 / count_C_piP + 1.0 / count_C_ele + 
                                               1.0 / count_LD2_piP + 1.0 / count_LD2_ele);
            graph_multiplicity_piP->SetPoint(pointIndex_piP, binCenter, R_h_A_piP);
            graph_multiplicity_piP->SetPointError(pointIndex_piP, 0, error_R_h_A_piP);
            pointIndex_piP++;
        }

        // Calculate multiplicity ratio for pi-
        if (count_LD2_piM > 0 && count_C_piM > 0) {
            R_h_A_piM = (count_C_piM / count_C_ele) / (count_LD2_piM / count_LD2_ele);
            error_R_h_A_piM = R_h_A_piM * sqrt(1.0 / count_C_piM + 1.0 / count_C_ele + 
                                               1.0 / count_LD2_piM + 1.0 / count_LD2_ele);
            graph_multiplicity_piM->SetPoint(pointIndex_piM, binCenter, R_h_A_piM);
            graph_multiplicity_piM->SetPointError(pointIndex_piM, 0, error_R_h_A_piM);
            pointIndex_piM++;
        }
    }
}

// Create a canvas to draw both graphs
TCanvas* canvas = new TCanvas("canvas_multiplicity", "Multiplicity Ratio #pi^{+} and #pi^{-} vs Q^{2}", 800, 600);
graph_multiplicity_piP->SetTitle("Multiplicity Ratio #pi^{+} and #pi^{-} vs Q^{2};Q^{2};R_{h}^{A}");
graph_multiplicity_piP->SetMarkerStyle(kFullCircle);
graph_multiplicity_piP->SetMarkerColor(kBlue);  // Color for pi+
graph_multiplicity_piP->SetLineColor(kBlue);    // Set error bar color to blue
graph_multiplicity_piP->Draw("AP");  // Draw pi+ data first

graph_multiplicity_piM->SetMarkerStyle(kFullSquare);
graph_multiplicity_piM->SetMarkerColor(kRed);  // Color for pi-
graph_multiplicity_piM->SetLineColor(kRed);    // Set error bar color to red
graph_multiplicity_piM->Draw("P");  // Draw pi- data on the same canvas

// Create the legend
TLegend* legend = new TLegend(0.15, 0.9, 0.3, 0.8);  // Set position (x1, y1, x2, y2)
legend->AddEntry(graph_multiplicity_piP, "#pi^{+}", "p");
legend->AddEntry(graph_multiplicity_piM, "#pi^{-}", "p");
legend->SetTextSize(0.03);
legend->SetFillStyle(0);   
legend->SetBorderSize(0);

legend->Draw();  // Draw the legend

canvas->SaveAs("multiplicity_ratio.png");


}

