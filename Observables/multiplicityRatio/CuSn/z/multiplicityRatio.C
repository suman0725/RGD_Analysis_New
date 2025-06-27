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



namespace fs = std::filesystem;
using namespace clas12;

void SetLorentzVector(TLorentzVector &p4, clas12::region_part_ptr rp) {
    p4.SetXYZM(rp->par()->getPx(), rp->par()->getPy(), rp->par()->getPz(), p4.M());
}

// Function to calculate Pt^2 for the produced hadron
/* double CalculatePt2(TLorentzVector P_h, TLorentzVector q) {
    // Get the angle between P_h and q
    double theta_a = P_h.Vect().Angle(q.Vect());

    // Calculate Pt^2 using |P_h|^2 * sin^2(theta_a)
    double Pt2 = P_h.P() * P_h.P() * std::sin(theta_a) * std::sin(theta_a);
    
    return Pt2;
}  */
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



// Function to compute φh
/* double CalculatePhih(TLorentzVector q, TLorentzVector p, TLorentzVector p_h) {
    // Extract 3-vectors from the Lorentz vectors
    TVector3 q3 = q.Vect();   // q vector
    TVector3 p3 = p.Vect();   // another particle's 
    TVector3 p_h3 = p_h.Vect(); // hadron's vector

    // Compute the cross products
    TVector3 cross_qp = q3.Cross(p3);      // (q x p)
    TVector3 cross_qph = q3.Cross(p_h3);   // (q x p_h)

    // Compute the dot product of the cross products
    double dotProduct = cross_qp.Dot(cross_qph);

    // Compute the magnitudes of the cross products
    double mag_qp = cross_qp.Mag();
    double mag_qph = cross_qph.Mag();

    // Calculate cos(φh)
    double cos_phi_h = dotProduct / (mag_qp * mag_qph);

    // Return φh by applying the inverse cosine (acos)
    return acos(cos_phi_h) * TMath::RadToDeg(); // result in radians
} */
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

    /*  Compute the dot product of the cross products
    double dotProduct = cross_qp.Dot(cross_qph);

    // Compute the magnitudes of the cross products
    double mag_qp = cross_qp.Mag();
    double mag_qph = cross_qph.Mag();

    // Calculate cos(φh)
    double cos_phi_h = dotProduct / (mag_qp * mag_qph);

    // Return φh by applying the inverse cosine (acos)
    return acos(cos_phi_h) * TMath::RadToDeg(); // result in degrees */
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
    TLorentzVector beam(0, 0, beamEnergy, beamEnergy); // Assuming the beam is along the z-axis
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
    // bin for variable nu 
    int nBins = 8;
    double zMin = 0.1;
    double zMax = 0.9;
    double LD2_ele_nor = 0; 
    double C_ele_nor =0;
   

    // Histograms for LD and Nuclear targets
    //auto* h_LD2_ele = new TH1F("h_LD2_ele", "LD Electrons;nu;Counts", nBins, zMin, zMax);
    auto* h_LD2_piP = new TH1F("h_LD2_piP", "LD Hadrons;z;Counts", nBins, zMin, zMax);
    //auto* h_C_ele = new TH1F("h_C_ele", "Carbon Target Electrons;nu;Counts", nBins, zMin, zMax);
    auto* h_C_piP = new TH1F("h_C_piP", "Carbon Positive Pions;z;Counts", nBins, zMin, zMax);


    /* auto* h_Sn_ele = new TH1F("h_Sn_ele", "Tin Electrons;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Sn_piP = new TH1F("h_Sn_piP", "Tin Positive Pions;nu;Counts", nBins, nuMin, nuMax);
   */

    // including Hipo files : LD2 
    clas12root::HipoChain chainLD2, chainC;
    //string baseDirectoryLD2 = "/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstLD2/dst/recon/018309";
    //addHipoFiles(chainLD2, baseDirectoryLD2);
    //string baseDirectoryC = "/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstCxC/dst/recon/018339"
    chainLD2.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstLD2/dst/recon/018309/rec_clas_018309.evio.00035-00039.hipo");
    chainC.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstCxC/dst/recon/018339/rec_clas_018339.evio.00015-00019.hipo");
    chainLD2.db()->turnOffQADB();
    chainC.db()->turnOffQADB();

    auto& c12 = chainLD2.C12ref();
    while (chainLD2.Next()) {
        auto particles = c12->getDetParticles();
        TLorentzVector el_temp, pip_temp;

        bool foundElectron = false;
        bool foundPion = false;

        double nu = 0; // Declare nu outside of the loop
        double Q2 = 0;
        
        
        TLorentzVector q; // Declare q outside of the loop

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
                double W2 = target.M() * target.M() - Q2 + 2 * target.M() * nu;
                double xB = Q2 / (2 * target.M() * nu);
                double y = nu / beam.Energy();

                

                /* h_LD2_nu->Fill(nu);
                h_LD2_Q2->Fill(Q2);
                h_LD2_xB->Fill(xB);
                h_LD2_W2->Fill(W2);  */

                //h_LD2_ele->Fill(Q2);
                LD2_ele_nor++;
            }

            // Pion selection
            if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pip_temp, p);
                foundPion = true;

                
            }
        }

        // Only calculate and fill histograms if both electron and pion are found
        if (foundElectron && foundPion) {
            
            /* double Pt2 = CalculatePt2(pip_temp, q); 
            h_LD2_Pt2->Fill(Pt2); */
            

            double z = pip_temp.E() / nu; // z = E_h / nu
            //h_LD2_z->Fill(z); // Fill the z histogram
            

            /* double Phi_ppi =  CalculatePhih(q,beam,pip_temp);
            h_LD2_Phi_h->Fill(Phi_ppi); */

            h_LD2_piP->Fill(z);


            
        } 
       
       /*  if (foundElectron && foundPion) {

        double Pt2 = CalculatePt2(pip_temp, q); 
            h_LD2_Pt2->Fill(Pt2);

            double Phi_ppi =  CalculatePhih(q,el_temp,pip_temp);
                    h_LD2_Phi_h->Fill(Phi_ppi);
       } */
    }
 
    
    // Process Carbon target
    auto& c12_C = chainC.C12ref();
    while (chainC.Next()) {
        auto particles = c12_C->getDetParticles();
        TLorentzVector el_temp_C, pip_temp_C;

        bool foundElectron_C = false;
        bool foundPion_C = false;

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
                double W2_C = target.M() * target.M() - Q2_C + 2 * target.M() * nu_C;
                double xB_C = Q2_C / (2 * target.M() * nu_C);

                /* h_C_nu->Fill(nu_C);
                h_C_Q2->Fill(Q2_C);
                h_C_xB->Fill(xB_C);
                h_C_W2->Fill(W2_C);
 */
                //h_C_ele->Fill(Q2_C);
                C_ele_nor ++; 
            }

            if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pip_temp_C, p);
                foundPion_C = true;
                 
            }
        }

        if (foundElectron_C && foundPion_C) {
            /* double Pt2_C = CalculatePt2(pip_temp_C, q_C);
            h_C_Pt2->Fill(Pt2_C); */
            double z_C = pip_temp_C.E() / nu_C;
            //h_C_z->Fill(z_C);
            /* double Phi_ppi_C = CalculatePhih(q_C, beam, pip_temp_C);
            h_C_Phi_h->Fill(Phi_ppi_C); */

            h_C_piP->Fill(z_C);

            
        } 
    } 


    TGraphErrors* graph = new TGraphErrors();
    int pointIndex = 0;

    for (int bin = 1; bin <= nBins; ++bin) {
        double binCenter = h_LD2_piP->GetBinCenter(bin);
        double count_LD2_ele =  LD2_ele_nor;        //h_LD2_ele->GetBinContent(bin);
        double count_LD2_piP = h_LD2_piP->GetBinContent(bin);
        double count_C_ele =  C_ele_nor;         //h_C_ele->GetBinContent(bin);
        double count_C_piP = h_C_piP->GetBinContent(bin);

        std::cout << "Bin Center: " << binCenter << ", Counts - "
          << "count_C_piP: " << count_C_piP << ", "
          << "count_C_ele: " << count_C_ele << ", "
          << "count_LD2_piP: " << count_LD2_piP << ", "
          << "count_LD2_ele: " << count_LD2_ele << std::endl;


        double R_h_A = 0.0;
        double error_R_h_A = 0.0;
        if (count_LD2_ele > 0 && count_C_ele > 0) {
            R_h_A = (count_C_piP / count_C_ele) / (count_LD2_piP / count_LD2_ele);

            error_R_h_A = R_h_A * sqrt(1.0 / count_C_piP + 1.0 / count_C_ele + 
                                        1.0 / count_LD2_piP + 1.0 / count_LD2_ele);
        
        graph->SetPoint(pointIndex, binCenter, R_h_A);
        graph->SetPointError(pointIndex, 0, error_R_h_A); // Set the error

        // Create a text object to display the values
        /* TLatex latex;
        latex.SetTextSize(0.03);  // Set text size
        latex.SetTextAlign(12);  // Align text left and center vertically
        latex.DrawLatex(binCenter, R_h_A, Form("R = %.2f, Error = %.2f", R_h_A, error_R_h_A)); */

        // Print the values
        std::cout << "Bin Center: " << binCenter
                  << ", R_h_A: " << R_h_A
                  << ", Error: " << error_R_h_A << std::endl;

        pointIndex++;
        }
        
    }


    // Create a canvas to draw the graph
    TCanvas* canvas = new TCanvas("canvas", "Multiplicity Ratio vs z", 800, 600);
    
    // Set graph titles and styles
    graph->SetTitle("Multiplicity Ratio vs z;z;R_h^A");
    graph->SetMarkerStyle(kFullCircle);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue); // Optional: set the line color

    // Draw the graph with error bars
    graph->Draw("AP");

    // Save the canvas as a PNG file
    canvas->SaveAs("multiplicity_ratio_z.png");


    
    // Drawing histograms for kinematic variables and saving as PDF for LD2
    /* TCanvas* can1 = new TCanvas();
    h_LD2_nu->Draw();
    can1->Print("Kinematic_Variables_LD2_018309.pdf(", "pdf");

    TCanvas* can2 = new TCanvas();
    h_LD2_Q2->Draw();
    can2->Print("Kinematic_Variables_LD2_018309.pdf", "pdf");

    TCanvas* can3 = new TCanvas();
    h_LD2_xB->Draw();
    can3->Print("Kinematic_Variables_LD2_018309.pdf", "pdf");

    TCanvas* can4 = new TCanvas();
    h_LD2_W2->Draw();
    can4->Print("Kinematic_Variables_LD2_018309.pdf", "pdf");

    // Old Pt2 histogram (log scale and x-axis range)
    TCanvas* can5 = new TCanvas();
    can5->SetLogy();            // Set y-axis to log scale
    h_LD2_Pt2->GetXaxis()->SetRangeUser(0, 3);  // Set x-axis range
    h_LD2_Pt2->Draw();
    can5->Print("Kinematic_Variables_LD2_018309.pdf", "pdf");

    TCanvas* can6 = new TCanvas();
    h_LD2_Phi_h->Draw();
    can6->Print("Kinematic_Variables_LD2_018309.pdf", "pdf");

    // New z histogram (log scale and x-axis range)
    TCanvas* can7 = new TCanvas();
    can7->SetLogy();            // Set y-axis to log scale
    h_LD2_z->GetXaxis()->SetRangeUser(0, 1); // Set x-axis range
    h_LD2_z->Draw();
    can7->Print("Kinematic_Variables_LD2_018309.pdf)", "pdf"); 


    // Drawing histograms for kinematic variables and saving as PDF for C (Carbon)
    TCanvas* can8 = new TCanvas();
    h_C_nu->Draw();
    can8->Print("Kinematic_Variables_C_018339.pdf(", "pdf");

    TCanvas* can9 = new TCanvas();
    h_C_Q2->Draw();
    can9->Print("Kinematic_Variables_C_018339.pdf", "pdf");

    TCanvas* can10 = new TCanvas();
    h_C_xB->Draw();
    can10->Print("Kinematic_Variables_C_018339.pdf", "pdf");

    TCanvas* can11 = new TCanvas();
    h_C_W2->Draw();
    can11->Print("Kinematic_Variables_C_018339.pdf", "pdf");

    TCanvas* can12 = new TCanvas();
    h_C_Pt2->Draw();
    can12->Print("Kinematic_Variables_C_018339.pdf", "pdf");

    TCanvas* can13 = new TCanvas();
    h_C_Phi_h->Draw();
    can13->Print("Kinematic_Variables_C_018339.pdf", "pdf");

    TCanvas* can14 = new TCanvas();
    h_C_z->Draw();
    can14->Print("Kinematic_Variables_C_018339.pdf)", "pdf");   */
}
