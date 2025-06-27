/*

Total four targets : Liquid deuterium (LD2), Carbon Carbon (CxC), Copper Tin (CuSn)
#NOTE: But, RGD used CuSn target for same beam , so  events from the 
both targets ( Cu and Sn) are stored in the same hipo files. Electron Vz needs to be applied as a cut to seperate them for the each independent target's events. 

***********************************************************************************************************************************
Multiplicity Ratio (R) Extraction for all three nuclear targets ( Cu, Sn, CxC) comparing with the liquid deuterium (LD2) target
***********************************************************************************************************************************
 
First, we have to load the hipo files for all targets from the desired directories. 
Second, looping over events for each independent target (for Cu and Sn we have to apply vz electron cut)
Base directory is: /lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v11
Inside base directory, three different directories according to different target : CuSn  CxC  LD2. Hipo files related to each target 
are stored in each associated target directories. But, CuSn is used for both target at the same time. 

Target configuration: Inbending (Torus = -1) and outbending (Torus = +1)
Each run has been triggered with the specific configuration and saved in the variable name "torus" inside the bank name "RUN::config". 
We have to analyze data from both configuration independently. 

Code should go to base directory, look the needed target directory, inside loop each hipo file, open each hipo file, loop through each events in each hipo file, 
open each event, load necessary bank ( REC::Particle, RUN::config for R #NOTE: I believe each run runnumber is saved as one hipo file, so I think checking "torus" 
from first event would be enough to separate that hipo file is inbending or outbening: 
#include "reader.h" 

Once we find inbending or outbending, (we will use variable form REC::Particles banks only onwards) first we have to find all the electrons from all events from all files in the target directory. Simply taking pid ==11 from REC::Particle bank
Then, selecting the trigger electron : 
1. status < 0, first row, -5< chi2pid < 5, vz cut, for forward detector abs(status)/1000==2, ( i think it would be better 
if we start tracking counting the electron, that means we have to track how many electrons are cutting off after each cut, litterally after each cut)
only the vz cut is different : 
These cuts are for outbending configuration electron vz
Target                            Vz range  (cm)
            (electron)              (+ve pions)                 (-ve pions)
Cu       [-10.5023, -6.4529]        [-9.76024, -5.1676]         [-10.9553, -5.8338]
Sn       [-6.0086, 5.0000]          [-4.6413, 5.0000]           [-5.4629, 5.0000]
LD2      [-20.0000, 5.0000]         [-20.0000, 5.0000]          [-20.0000,5.0000]
CxC      [-10.5319, 5.0000]         [-10.5319, 5.0000]          [-10.3501, 5.0000]

I still need to get specific cuts for inbending configuration for electron, +ve pions, -ve pions vz 

for +ve pion 

After having trigger electron: calculating all the associated kinematic variables : nu, q, Q^2, y, W
for kinematical variable cuts: Q2 > 1 && y > 0.25 && y < 0.85 && W > 2 , we can refine  trigger electrons after applying these cuts

Now, we have to choose only  event which has trigger electron ( what we expect is each event has one trigger electron): 
And inside event with the trigger electron, select the pions (i have to do for +ve and -ve): while selecting the pions,
i think for tracking purpose, we can count all +ve/-ve pions before using trigger electron cut, count after using trigger electron cut, then apply vz and chi2pid cut and 
keep track of pions. After applying all the cuts to select +ve/-ve pions, now calculate kinematica variable related to hadron: z, Pt
Kinematic varialbes specific to the hadron for TMDs: 0.3 < z < 0.7 and p_t^2 < 1.2 GeV^2 


***********************************************************************************************************************************
Definition of R
***********************************************************************************************************************************
R_M (z, nu, Q^2, p_t^2, phi) = [ N_h ( z, nu, Q^2, pt^2, phi )/N_DIS (nu, Q^2) ]_A / [ N_h ( z, nu, Q^2, pt^2, phi )/N_DIS (nu, Q^2) ]_D
// A-> C, Cu, Sn ; D-> LD2
1. R_M (nu) = [N_h(nu)/ N_e(nu)]_A/[N_h(nu)/ N_e(nu)]_D {R_M dependance with variable nu, we have to }
For each bin of nu, counting the number of hadrons and electron  and get R_M in each bin  by taking the ratio
2. R_M (Q^2) = [N_h(Q^2)/ N_e(Q^2)]_A/[N_h(Q^2)/N_e(Q^2)]_D {R_M dependance with variable Q^2}
For each bin of Q^2, counting the number of hadrons and electron  and get R_M in each bin by taking the ratio 
3. R_M (z) = [N_h(z)/ N_e]_A/[N_h(z)/ N_e]_D
For each bin of Z, counting the number of hadrons and here N_e is the normalized electron not z dependend , and get R_M in each bin by taking the ratio 
4. R_M (p_t^2) = [N_h(p_t^2)/ N_e]_A/[N_h(p_t^2)/ N_e]_D
For each bin of p_t^2, counting the number of hadrons and here N_e is the normalized electron not p_t^2 dependend , and get R_M in each bin by taking the ratio 

[N_h, N_e]_A: Number of SIDS hadron h and DIS electrons on a target A
[N_h, N_e]_D: Number of SIDS hadron h and DIS electrons on a target LD2 



*/


// for example here is the example for 

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





//we can define required four momentum vector for incoming electron, target proton, scattered electron, produced +ve pions , produced -ve pions
 auto db = TDatabasePDG::Instance();
// for eg. we can get mass of any particle in pdd with pid like this: db->GetParticle(2212)->Mass()
//we can use beam energy 
double beamEnergy = 10.532; // GeV
// Bin for variables nu, Q^2, z, p_t^2 
//for Q^2 ; nBins = 10, min 1, max = 9 
//for nu, int nBins = 10;min= 0.1 , max= 1
//for pt_2 int nBins = 10;min = 0, max 1.5 
// just make function to create bins before main function and apply for all variable 

    // Histograms for LD2, Cu, Sn, carbon targets (positive and negative pions): 
    auto* h_LD2_ele = new TH1F("h_LD2_ele", "LD Electrons;nu;Counts", nBins, nuMin, nuMax);
    auto* h_LD2_piP = new TH1F("h_LD2_piP", "LD Positive Pions;nu;Counts", nBins, nuMin, nuMax);
    auto* h_LD2_piM = new TH1F("h_LD2_piM", "LD Negative Pions;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Cu_ele = new TH1F("h_Cu_ele", "Copper Electrons;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Cu_piP = new TH1F("h_Cu_piP", "Copper Positive Pions;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Cu_piM = new TH1F("h_Cu_piM", "Copper Negative Pions;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Sn_ele = new TH1F("h_Sn_ele", "Tin Electrons;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Sn_piP = new TH1F("h_Sn_piP", "Tin Positive Pions;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Sn_piM = new TH1F("h_Sn_piM", "Tin Negative Pions;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Carbon_ele = new TH1F("h_Carbon_ele", "Carbon Electrons;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Carbon_piP = new TH1F("h_Carbon_piP", "Carbon Positive Pions;nu;Counts", nBins, nuMin, nuMax);
    auto* h_Carbon_piM = new TH1F("h_Carbon_piM", "Carbon Negative Pions;nu;Counts", nBins, nuMin, nuMax);

    //this is for nu variables but we have to create for rest of three variables 

    // Load Hipo files
    //we have to  find the strategy load hipo files, if the script beccomes complex then we can make four different script for four target, suggest me here
    

    //for event loop
     for (const auto& file : hipoFiles) {
        cout << "Processing file: " << file << endl;
        hipo::reader reader;
        reader.open(file.c_str());
        hipo::dictionary factory;
        reader.readDictionary(factory);

        while (reader.next() && event_count < maxEvents) {
            event_count++;
            hipo::event event;
            reader.read(event);

            hipo::bank PART(factory.getSchema("REC::Particle"));
              hipo::bank RUN(factory.getSchema("RUN::config"));
              // i think first event is enough to see inbending or outbending

        // to extract all the information from rec::particle bank we can make
        ParticleData getParticleData(int partIdx, hipo::bank& PART,  hipo::bank& RUN,
                             ) {
            ParticleData pd;
            pd.px = PART.getFloat("px", partIdx);
            pd.py = PART.getFloat("py", partIdx);
            pd.pz = PART.getFloat("pz", partIdx);
            pd.p = sqrt(pd.px * pd.px + pd.py * pd.py + pd.pz * pd.pz);
            pd.vz = PART.getFloat("vz", partIdx);
            pd.status = PART.getShort("status", partIdx);
            pd.pid = PART.getInt("pid", partIdx);
            pd.chi2pid = PART.getFloat("chi2pid", partIdx); }
        // one common function before main function once we use that function 
         //something like this 
        // Extract particle data from HIPO banks

        // all information stores when call from inside for all targets 
       
        // Extract particles and count electrons/positrons
            vector<ParticleData> particles(PART.getRows());
            for (int i = 0; i < PART.getRows(); ++i) {
                particles[i] = getParticleData(i, PART, RUN);
            }

            for (const auto& p : particles) {

                // Electron selection
                if (pid == 11 ){
                    //total electrons
                
                    if (status < 0/* and first row of REC*/ && chi2Pid < 5 && chi2Pid > -5 /* aaply vz cut according to targe*/) {
                
                    
                        nu = beam.Energy() - el_temp.Energy();
                        q = beam - el_temp;
                        Q2 = -q.Mag2(); 
                        double xB = Q2 / (2 * target.M() * nu);
                        //add y, W double W = sqrt(target.M() * target.M() + 2 * target.M() * nu - Q2);
                            //double y = nu / beam.E();
                        //apply those  cuts Q2 > 1 && y > 0.25 && y < 0.85 && W > 2 >> hastrigger
                        //fill ( nu) to get [N_e(nu)]_A if C, Cu, Sn, [N_e(nu)]_D if LD2
                        //fill ( Q^2) to get [N_e(Q^2)]_A if C, Cu, Sn, [N_e(Q^2)]_D if LD2
                        //if z dependence just count the elctron here [N_e]_A, [N_e]_D
                        //ifp_t^2 dependence just count the elctron here [N_e]_A, [N_e]_D
                        //fill ny for nu variable 
                        //we have to for all target DIS
                    }

                }
            }

   
    if (!hasTrigger) continue;
            for (const auto& p : particles) {
                if (abs(p.status) / 1000 != 2) continue;
                //collect +ve pions and -ve pions  using 
                if (pid == 211){

                    if ( /*vz cut and chi2pid cut */){
                        //calculate z and p_t^2 
                        //double z = pip_temp.E() / nu; // z = E_h / nu
                        //I think we can add the function to calculate the p_t^2 before main function 
                        /*double CalculatePt2(TLorentzVector P_h, TLorentzVector q) {

                            TVector3 q3 = q.Vect(); 
                            TVector3 p_h3 = P_h.Vect(); 

                            TVector3 cross_qph = q3.Cross(p_h3);

                            // Calculate the magnitudes squared
                            double magnitudeCrossProductSquared = cross_qph.Mag2(); // |q × p_h|^2
                            double magnitudeQSquared = q3.Mag2();                     // |q|^2

                            // Calculate Pt^2 = |q × p_h|^2 / |q|^2
                            double Pt2 = magnitudeCrossProductSquared / magnitudeQSquared;

                            return Pt2;
                            apply the condition 0.3 < z < 0.7 and p_t^2 < 1.2 GeV^2 not at the same time but according to necessity 
                            Fill (Q^2) here if dependence with Q^2 for [N_h(Q^2)]_A and [N_h(Q^2)]_D
                            Fill (nu) here if dependence with nu for [N_h(nu)]_A and [N_h(nu)}_D

                            Fill ( z) for [N_h(z)]_A and [N_h(z)]_D
                            fill(pt^2) for [N_h(p_t^2)]_A and [N_h(p_t^2)]_D 


                        }*/

                    }
                }
   
        // Process only for Carbon target (assume Carbon target Vz range if necessary)
        //We have to finish upto this like for all target ok 
        //But count for all the target which are filled will be counted for each bin in next part if 
        // this is not the correct approach feel free to modify ok 
        // for those electrons counting for z and pt^2 we get from before part 
        //ultimately we should just run the code and get the result
        //and code should be very simple to understand 
           


    // Process multiplicity ratio calculation and graphing for Cu and Sn in a similar way...
   
    // Create TGraphErrors for Cu target
    //from this below this only for nu variables so, 
    //but we have to for rest three Q^2, Z, p_t^2 
    //before getting R, we should have the counting for each bin ( and just total electrons number if z and pt ^2  dependance )
    //once we have all the counting for A and D then just take the ratio to get R
    // and include the error on R
    // and put in the same TGraph for C, Cu, and Sn 
    // we have to do for inbending and outbending 
    //we have to do for +ve pions and -ve pions
    //That means should have 4 mulitplicty plot for +ve pion for all target 
   // and  have 4 mulitplicty plot for -ve pion for all target 
   // for 8 for each configuration right?
   //what should be the correct strategy so that script should looks simple but working all task what do you suggest 
/* TGraphErrors* graph_multiplicity_piP_Cu = new TGraphErrors();
TGraphErrors* graph_multiplicity_piP_Sn = new TGraphErrors();
TGraphErrors* graph_multiplicity_piP_Carbon = new TGraphErrors();


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


} */




