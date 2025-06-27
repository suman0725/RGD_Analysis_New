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

double CalculatePt2(TLorentzVector P_h, TLorentzVector q) {
    TVector3 q3 = q.Vect(); 
    TVector3 p_h3 = P_h.Vect(); 
    TVector3 cross_qph = q3.Cross(p_h3);
    double magnitudeCrossProductSquared = cross_qph.Mag2(); // |q × p_h|^2
    double magnitudeQSquared = q3.Mag2();                     // |q|^2
    return magnitudeCrossProductSquared / magnitudeQSquared; // Pt^2
}




/* void addHipoFiles(clas12root::HipoChain& chain, const fs::path& baseDir) {
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
} */

void tmb() {
    auto db = TDatabasePDG::Instance();
    double beamEnergy = 10.532; // GeV
    TLorentzVector beam(0, 0, beamEnergy, beamEnergy); // Beam vector
    TLorentzVector target(0, 0, 0, db->GetParticle(2212)->Mass()); // Proton target
    TLorentzVector el, pip;

    // Histograms for kinematic variables fo LD2
    /* auto* h_LD2_nu = new TH1F("h_LD2_nu", "LD2 nu;nu(GeV);Counts", 200, 1, 11);
    auto* h_LD2_Q2 = new TH1F("h_LD2_Q2", "LD2 Q^{2} distribution;Q^{2} (GeV)^{2};Counts", 200, 0, 8);
    auto* h_LD2_xB = new TH1F("h_LD2_xB", "LD2 Bjorken x distribution;xB;Counts", 100, 0, 1);
    auto* h_LD2_W2 = new TH1F("h_LD2_W2", "LD2 W^{2} distribution;W^{2} (GeV)^{2};Counts", 200, 0, 20);
    auto* h_LD2_Pt2 = new TH1F("h_LD2_Pt2", "LD2 Pt^2 distribution;Pt^2;Counts", 100, 0, 3);
    auto* h_LD2_z = new TH1F("h_LD2_z", "LD2 z distribution;z;Counts", 100, 0, 1); // New histogram for z
    auto* h_LD2_Phi_h = new TH1F("Ph_LD2_hi_h", "LD2 Phi_h distribution;Phi_h (Degree);Counts", 100, 0, 180); // New histogram for z */
    
    // Histograms for kinematic variables for Nuclear Target Carbon
    /* auto* h_C_nu = new TH1F("h_C_nu", "Carbon nu;nu(GeV);Counts", 200, 1, 11);
    auto* h_C_Q2 = new TH1F("h_C_Q2", "Carbon  Q^{2} distribution;Q^{2} (GeV)^{2};Counts", 200, 0, 8);
    auto* h_C_xB = new TH1F("h_C_xB", "Carbon  Bjorken x distribution;xB;Counts", 100, 0, 1);
    auto* h_C_W2 = new TH1F("h_C_W2", "Carbon  W^{2} distribution;W^{2} (GeV)^{2};Counts", 200, 0, 20);
    auto* h_C_Pt2 = new TH1F("h_C_Pt2", "Carbon  Pt^2 distribution;Pt^2;Counts", 100, 0, 3);
    auto* h_C_z = new TH1F("h_C_z", "Carbon  z distribution;z;Counts", 100, 0, 1); // New histogram for z
    auto* h_C_Phi_h = new TH1F("h_C_Phi_h", "Carbon  Phi_h distribution;Phi_h (Degree);Counts", 100, 0, 180); // New histogram for z 
 */

    // Histograms for kinematic variables for Nuclear Target Tin
    /* auto* h_Sn_nu = new TH1F("h_Sn_nu", "Tin LD2_nu;nu(GeV);Counts", 200, 1, 11);
    auto* h_Sn_Q2 = new TH1F("h_Sn_Q2", "Tin Q^{2} distribution;Q^{2} (GeV)^{2};Counts", 200, 0, 8);
    auto* h_Sn_xB = new TH1F("h_Sn_xB", "Tin Bjorken x distribution;xB;Counts", 100, 0, 1);
    auto* h_Sn_W2 = new TH1F("h_Sn_W2", "Tin W^{2} distribution;W^{2} (GeV)^{2};Counts", 200, 0, 20);
    auto* h_Sn_Pt2 = new TH1F("h_Sn_Pt2", "Tin Pt^2 distribution;Pt^2;Counts", 100, 0, 3);
    auto* h_Sn_z = new TH1F("h_Sn_z", "Tin z distribution;z;Counts", 100, 0, 1); // New histogram for z
    auto* h_Sn_Phi_h = new TH1F("h_Sn_Phi_h", "Tin Phi_h distribution;Phi_h (Degree);Counts", 100, 0, 180); // New histogram for z
 */
    /* auto* h_LD2_z = new TH1F("h_LD2_z", "LD2 z distribution;z;Counts", 100, 0, 1); // New histogram for z
    auto* h_C_z = new TH1F("h_C_z", "Carbon  z distribution;z;Counts", 100, 0, 1.5); // New histogram for z */
    auto* h_LD2_pT2 = new TH1F("h_LD2_Pt2", "LD2 Pt^2 distribution;Pt^2;Counts", 100, 0, 2);
    auto* h_C_pT2 = new TH1F("h_C_Pt2", "Carbon  Pt^2 distribution;Pt^2;Counts", 100, 0, 1);
    auto* h_C_pT2_average = new TProfile("h_C_Pt2", "Carbon  Pt^2 distribution;Pt^2;Counts", 100, 0, 1);
    auto* h_C_pT4_average = new TProfile("p_C_pT4_Squared", "Average (p_T^2)^2 Profile; p_T^2 [GeV^2]; Average (p_T^2)^2", 100, 0, 1);
    auto* h_LD2_pT2_average = new TProfile("h_LD2_Pt2", "Carbon  Pt^2 distribution;Pt^2;Counts", 100, 0, 1);
    auto* h_LD2_pT4_average = new TProfile("p_LD2_pT4_Squared", "Average (p_T^2)^2 Profile; p_T^2 [GeV^2]; Average (p_T^2)^2", 100, 0, 1);
    // Histograms for Q^2 for pions
    /* auto* h_LD2_Q2_piP = new TH1F("h_LD2_Q2_pion", "LD2 Q^{2} (Pions);Q^{2} (GeV)^{2};Counts", 200, 0, 8.5);
    auto* h_C_Q2_piP = new TH1F("h_C_Q2_pion", "Carbon Q^{2} (Pions);Q^{2} (GeV)^{2};Counts", 200, 0, 8.5); */
    


    /* // Open output file
    std::ofstream outFile("output.txt");
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open the output file!" << std::endl;
        return;
    } */
    
    // Binning for Q^2
    int nBins = 5;
    double pT2Min = 0.0;
    double pT2Max = 1.5;
    double binWidth = (pT2Max - pT2Min) / nBins;

    std::vector<double> pT2_D(nBins, 0.0); // Store <Pt^2> for LD2
    std::vector<double> pT2_A(nBins, 0.0); // Store <Pt^2> for Carbon
    std::vector<double> pT4_D(nBins, 0.0); // Store <Pt^4> for LD2
    std::vector<double> pT4_A(nBins, 0.0); // Store <Pt^4> for Carbon
    std::vector<int> counts_D(nBins, 0);   // Count events for LD2
    std::vector<int> counts_A(nBins, 0);   // Count events for Carbon
    int totalEntries = 0;
    int totalEntries_C = 0;
    std::vector<double> errors_C(nBins, 0.0);
    std::vector<double> errors_D(nBins, 0.0);
    std::vector<double> errors(nBins, 0.0);  // Standard error for each bin
    std::vector<double> deltapT2(nBins, 0.0);
    


    // Hipo files
    clas12root::HipoChain chainLD2, chainC;
   /*  chainLD2.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstLD2/dst/recon/018309/rec_clas_018309.evio.00035-00039.hipo");
    chainC.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstCxC/dst/recon/018339/rec_clas_018339.evio.00015-00019.hipo"); */
    chainLD2.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v6/LD2/skim_run_018307.hipo");
    chainC.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v6/CxC/skim_run_018339.hipo");
    chainLD2.db()->turnOffQADB();
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

            // If both electron and pion are found, calculate Pt^2
            if (foundPion) {
                //double z = pip_temp.E() / nu;
                if ( nu < 2.5 && nu >= 1.5 ){
                    double pT2_value = CalculatePt2(pip_temp, q);
                    double pT4_value = pT2_value * pT2_value;
                    h_LD2_pT2->Fill(pT2_value);
                    h_LD2_pT2_average->Fill(pT2_value, pT2_value);
                    h_LD2_pT4_average->Fill(pT4_value, pT4_value);
                }
                
            }  
        }
    }

    double wightedSumpT2 = 0.0; 
    double wightedSumpT4 = 0.0; 
    double N_entries = 0.0; 
    for (int i = 1; i <= h_LD2_pT2_average->GetNbinsX(); ++i) {
        double pT2_average = h_LD2_pT2_average->GetBinContent(i);
        double pT4_average = h_LD2_pT4_average->GetBinContent(i);
        double N_i = h_LD2_pT2_average->GetBinEntries(i);

        // Calculate (p_T^2)_i * N_i
        wightedSumpT2 += pT2_average * N_i; 
        wightedSumpT4 += pT4_average * N_i; 
        N_entries += N_i; // Total number of entries
    }
        double average_pT2 = (N_entries > 0) ? (wightedSumpT2 / N_entries) : 0.0;
        double average_pT4 = (N_entries > 0) ? (wightedSumpT4 / N_entries) : 0.0;
        double variance = average_pT4 - average_pT2 * average_pT2;
        double error_D = (N_entries > 0 && variance >= 0) ? sqrt(variance / N_entries) : 0.0;
        cout << "Average p_T^2_LD2= " << average_pT2 << " [GeV^2]" << std::endl;
        cout << "Error_LD2= " << error_D << endl; 




        
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
                    double pT2_value_C = CalculatePt2(pip_temp_C, q_C);
                    double pT4_value_C = pT2_value_C * pT2_value_C;
                    h_C_pT2->Fill(pT2_value_C);
                    h_C_pT2_average->Fill(pT2_value_C, pT2_value_C);
                    h_C_pT4_average->Fill(pT4_value_C, pT4_value_C);
                }
            }
        }
    }

    double wightedSumpT2_C = 0.0; 
    double wightedSumpT4_C = 0.0; 
    double N_entries_C = 0.0; 
    for (int i = 1; i <= h_C_pT2_average->GetNbinsX(); ++i) {
        double pT2_average_C = h_C_pT2_average->GetBinContent(i);
        double pT4_average_C = h_C_pT4_average->GetBinContent(i);
        double N_i_C = h_C_pT2_average->GetBinEntries(i);

        // Calculate (p_T^2)_i * N_i
        wightedSumpT2_C += pT2_average_C * N_i_C; 
        wightedSumpT4_C += pT4_average_C * N_i_C; 
        N_entries_C += N_i_C; // Total number of entries
    }
    double average_pT2_A = (N_entries_C > 0) ? (wightedSumpT2_C / N_entries_C) : 0.0;
    double average_pT4_A = (N_entries_C > 0) ? (wightedSumpT4_C / N_entries_C) : 0.0;
    double variance_C = average_pT4_A - average_pT2_A * average_pT2_A;
    double error_A = (N_entries_C > 0 && variance_C >= 0) ? sqrt(variance_C / N_entries_C) : 0.0;

    double deltapt2 = average_pT2_A - average_pT2; 
    double error = sqrt(error_A * error_A + error_D * error_D); 
    
    cout << "Average p_T^2_C= " << average_pT2_C << " [GeV^2]" << std::endl;
    cout << "Error_C= " << error_C << endl; 

    cout << "deltapT2 is :" << deltapt2 << endl; 
    cout << "Error is :" << error << endl;

    /* outFile << deltapT2 << "\t" << error << std::endl;
    outFile.close(); 
 */

    TCanvas* c1 = new TCanvas("c1", "pT^2 Histogram", 800, 600);
    c1->Divide(1,2);
    c1->cd(1);
    h_C_pT2->Draw(); // Draw the histogram
    
    c1->cd(2);
    h_LD2_pT2->Draw();
    c1->SaveAs("pt_broadening_histogram_C_nu_1.15_2.5.png"); // Save the histogram as an image file



    /* TCanvas* c2 = new TCanvas("c1", "z Histogram", 800, 600);
    h_C_z->Draw(); // Draw the histogram
    c2->SaveAs("z_histogram_C.png");  */
       
    



    
}
