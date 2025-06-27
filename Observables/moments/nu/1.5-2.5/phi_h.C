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
#include <TProfile.h>
#include <fstream>


namespace fs = std::filesystem;
using namespace clas12;

void SetLorentzVector(TLorentzVector &p4, clas12::region_part_ptr rp) {
    double mass = TDatabasePDG::Instance()->GetParticle(rp->par()->getPid())->Mass();
    p4.SetXYZM(rp->par()->getPx(), rp->par()->getPy(), rp->par()->getPz(), mass);
}


double CalculatePhih(TLorentzVector q, TLorentzVector p, TLorentzVector p_h) {
    // Extract 3-vectors from the Lorentz vectors
    TVector3 q3 = q.Vect();   // q vector
    TVector3 p3 = p.Vect();   // incident electron
    TVector3 p_h3 = p_h.Vect(); // hadron's vector

    // Normalize the q vector
   TVector3 VectorQNorm = q3.Unit();

    // Compute the cross products
    TVector3 cross_qp = VectorQNorm.Cross(p3);      // (q x p)
    TVector3 cross_qph = VectorQNorm.Cross(p_h3);   // (q x p_h)
    TVector3 cross_pph = p3.Cross(p_h3);
   
    // Calculate the angle in degrees
    // Calculate the angle according to the Trento convention
    double phi = TMath::ATan2(cross_pph.Dot(VectorQNorm), cross_qp.Dot(cross_qph)) * TMath::RadToDeg();
    phi = phi + 180.0;

    return phi; 
}


void phi_h() {
    auto db = TDatabasePDG::Instance();
    double beamEnergy = 10.532; // GeV
    TLorentzVector beam(0, 0, beamEnergy, beamEnergy); // Beam vector
    TLorentzVector target(0, 0, 0, db->GetParticle(2212)->Mass()); // Proton target
    TLorentzVector el, pip;

    
    
    auto* h_LD2_Phi = new TH1F("h_LD2_Phi", "LD2 Phi_h distribution;Phi [deg];Counts", 100, 0, 360);
    auto* h_C_Phi = new TH1F("h_C_Phi", "Carbon Phi_h distribution;Phi [deg];Counts", 100, 0, 360);
    auto* h_C_Cos_Phi_average = new TProfile("h_C_Phi_average", "Carbon Phi Average Profile; Phi [deg]; Average Phi", 100, 0, 360);
    auto* h_C_Cos_Phi2_average = new TProfile("h_C_Phi2_average", "Average (Phi)^2 Profile; Phi [deg]; Average (Phi)^2", 100, 0, 130000);
    auto* h_LD2_Cos_Phi_average = new TProfile("h_LD2_Phi_average", "LD2 Phi Average Profile; Phi [deg]; Average Phi", 100, 0, 360);
    auto* h_LD2_Cos_Phi2_average = new TProfile("h_LD2_Phi2_average", "Average (Phi)^2 Profile; Phi [deg]; Average (Phi)^2", 100, 0, 130000);

    // Binning for Phi_h
    int nBins = 72; // Number of bins
    double PhiMin = 0.0; // Minimum value of Phi in degrees
    double PhiMax = 360.0; // Maximum value of Phi in degrees
    double binWidth = (PhiMax - PhiMin) / nBins; // Calculate bin width

    std::vector<double> Phi_D(nBins, 0.0);  // Store <Phi> for LD2
    std::vector<double> Phi_A(nBins, 0.0);  // Store <Phi> for Carbon
    std::vector<double> Phi2_D(nBins, 0.0); // Store <Phi^2> for LD2
    std::vector<double> Phi2_A(nBins, 0.0); // Store <Phi^2> for Carbon


    // Hipo files
    clas12root::HipoChain chainLD2, chainC;
    /* chainLD2.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstLD2/dst/recon/018309/rec_clas_018309.evio.00035-00039.hipo");
    chainC.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstCxC/dst/recon/018339/rec_clas_018339.evio.00015-00019.hipo"); */
    chainLD2.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v6/LD2/skim_run_018307.hipo");
    chainC.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v6/CxC/skim_run_018339.hipo");
    chainLD2.db()->turnOffQADB();
    chainC.db()->turnOffQADB();

    // Process LD2 target
    auto& c12 = chainLD2.C12ref();
    while (chainLD2.Next()) {
        auto particles = c12->getDetParticles();
        TLorentzVector el_temp, pip_temp;
        bool foundElectron = false, foundPion = false;

        for (auto& p : particles) {
            int pid = p->par()->getPid();
            int status = p->par()->getStatus();
            int chi2Pid = p->par()->getChi2Pid();

            // Check for electron
            if (pid == 11 && status < 0 && chi2Pid < 5 && chi2Pid > -5) {
                SetLorentzVector(el_temp, p);
                foundElectron = true;
            }

            // Check for pion
            if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pip_temp, p);
                foundPion = true;
            }
        }

        // Calculate nu, q, and Q^2 if an electron is found
        if (foundElectron) {
            double nu = beam.Energy() - el_temp.Energy(); 
            TLorentzVector q = beam - el_temp; 
            //double Q2 = -q.Mag2();
            if (foundPion) {
                //double z = pip_temp.E() / nu;
                if ( nu < 2.5 && nu >= 1.5 ){
                    double Phi = CalculatePhih(q,el_temp,pip_temp);
                    double Cos_Phi = cos(Phi);
                    double Cos_Phi2 = Cos_Phi * Cos_Phi;
                    h_LD2_Phi->Fill(Phi);
                    h_LD2_Cos_Phi_average->Fill(Cos_Phi, Cos_Phi);
                    h_LD2_Cos_Phi2_average->Fill(Cos_Phi2, Cos_Phi2);
                }
                
            }  
        }
    }

    double weightedSumCos_Phi = 0.0;
    double weightedSumCos_Phi2 = 0.0;
    double N_entries = 0.0;

    for (int i = 1; i <= h_LD2_Cos_Phi_average->GetNbinsX(); ++i) {
        double Cos_Phi_average = h_LD2_Cos_Phi_average->GetBinContent(i); // <Phi> for bin i
        double Cos_Phi2_average = h_LD2_Cos_Phi2_average->GetBinContent(i); // <Phi^2> for bin i
        double N_i = h_LD2_Cos_Phi_average->GetBinEntries(i); // Number of entries in bin i
        

        // Calculate weighted sum of <Phi> and <Phi^2>
        weightedSumCos_Phi += Cos_Phi_average * N_i;
        weightedSumCos_Phi2 += Cos_Phi2_average * N_i;
        N_entries += N_i; // Total number of entries
    }

    // Calculate averages and error
    double average_Cos_Phi = (N_entries > 0) ? (weightedSumCos_Phi / N_entries) : 0.0;
    double average_Cos_Phi2 = (N_entries > 0) ? (weightedSumCos_Phi2 / N_entries) : 0.0;
    double variance = average_Cos_Phi2- average_Cos_Phi * average_Cos_Phi;
    double error_LD2 = (N_entries > 0 && variance >= 0) ? sqrt(variance / N_entries) : 0.0;

    // Output results
    cout << "< Cos Phi_LD2 >= " << average_Cos_Phi << " [degrees]" << std::endl;
    
    cout << "Error < Cos Phi_LD2 >= " << error_LD2 << " [degrees]" << endl;





        
    // Process Carbon target
    auto& c12_C = chainC.C12ref();
    while (chainC.Next()) {
        auto particles = c12_C->getDetParticles();
        TLorentzVector el_temp_C, pip_temp_C;
        bool foundElectron_C = false, foundPion_C = false;

        for (auto& p : particles) {
            int pid = p->par()->getPid();
            int status = p->par()->getStatus();
            int chi2Pid = p->par()->getChi2Pid();

            // Check for electron
            if (pid == 11 && status < 0 && chi2Pid < 5 && chi2Pid > -5) {
                SetLorentzVector(el_temp_C, p);
                foundElectron_C = true;
            }

            // Check for pion
            if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pip_temp_C, p);
                foundPion_C = true;
            }
        }

        // Calculate nu, q, and Q^2 if an electron is found
        if (foundElectron_C) {
            double nu_C = beam.Energy() - el_temp_C.Energy(); 
            TLorentzVector q_C = beam - el_temp_C; 
            //double Q2_C = -q_C.Mag2();

            // If both electron and pion are found, calculate Pt^2
            if (foundPion_C) {
                //double z_C = pip_temp_C.E() / nu_C;   
               // h_C_z->Fill(z_C);     
                   
                if ( nu_C < 2.5 && nu_C >= 1.5 ){
                    double Phi_C = CalculatePhih(q_C, el_temp_C, pip_temp_C);  // Calculate Cos_Phi for Carbon
                    double Cos_Phi_C = cos(Phi_C);
                    double Cos_Phi2_C = Cos_Phi_C * Cos_Phi_C;         // Calculate Cos_Phi^2 for Carbon

                    h_C_Phi->Fill(Phi_C);                                  // Fill histogram for Cos_Phi
                    h_C_Cos_Phi_average->Fill(Cos_Phi_C, Cos_Phi_C);         // Fill profile for <Cos_Phi>
                    h_C_Cos_Phi2_average->Fill(Cos_Phi2_C, Cos_Phi2_C);      // Fill profile for <Cos_Phi^2>
 // Fill profile for <Phi^2>

                }
            }
        }
    }

    double weightedSumCos_Phi_C = 0.0; 
    double weightedSumCos_Phi2_C = 0.0; 
    double N_entries_C = 0.0; 

    for (int i = 1; i <= h_C_Cos_Phi_average->GetNbinsX(); ++i) {
        double Cos_Phi_average_C = h_C_Cos_Phi_average->GetBinContent(i);   // <Cos_Phi> for bin i
        double Cos_Phi2_average_C = h_C_Cos_Phi2_average->GetBinContent(i); // <Cos_Phi^2> for bin i
        double N_i_C = h_C_Cos_Phi_average->GetBinEntries(i);               // Number of entries in bin i

        // Calculate weighted sums
        weightedSumCos_Phi_C += Cos_Phi_average_C * N_i_C; 
        weightedSumCos_Phi2_C += Cos_Phi2_average_C * N_i_C; 
        N_entries_C += N_i_C; // Total number of entries
    }

    double average_Cos_Phi_C = (N_entries_C > 0) ? (weightedSumCos_Phi_C / N_entries_C) : 0.0; // Average <Cos_Phi>
    double average_Cos_Phi2_C = (N_entries_C > 0) ? (weightedSumCos_Phi2_C / N_entries_C) : 0.0; // Average <Cos_Phi^2>
    double variance_C = average_Cos_Phi2_C - average_Cos_Phi_C * average_Cos_Phi_C; // Variance
    double error_C = (N_entries_C > 0 && variance_C >= 0) ? sqrt(variance_C / N_entries_C) : 0.0; // Error

     // Output results
    cout << "< Cos_Phi_C > = " << average_Cos_Phi_C << " [degrees]" << std::endl;
    cout << "Error < Cos_Phi_C > = " << error_C << " [degrees]" << std::endl;


    double Ratio = average_Cos_Phi_C / average_Cos_Phi;
    double Error_Ratio= Ratio * sqrt(pow(error_C / average_Cos_Phi_C, 2) + pow(error_LD2 / average_Cos_Phi, 2));

    cout << " < Cos_Phi_C > / < Cos_Phi_LD2 > = " << Ratio << endl;
    cout << "Error (< Cos_Phi_C > / < Cos_Phi_LD2 >) = " << Error_Ratio << endl;



   
 
   /*  TCanvas* c1 = new TCanvas("c1", "Phi_h Histogram", 800, 600);
    c1->Divide(1, 2);
    c1->cd(1);
    h_C_Phi->Draw(); 
    gPad->SetGrid(); 
    c1->cd(2);
    h_LD2_Phi->Draw(); 
    gPad->SetGrid();
    c1->SaveAs("Phi_histogram_C.png"); */
 






    
}
