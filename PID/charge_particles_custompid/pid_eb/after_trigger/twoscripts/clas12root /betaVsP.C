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
#include <TGraph.h>
#include "clas12reader.h"
#include "HipoChain.h"
#include <TLine.h>
using namespace clas12;
using namespace std;

void betaVsP() {     

// Process LD Target Files
    clas12root::HipoChain chainLD;
    chainLD.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v9/LD2/skim_run_018536.hipo");
    chainLD.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v9/LD2/skim_run_018537.hipo");
    chainLD.db()->turnOffQADB();
    auto& c12LD = chainLD.C12ref();

    int status; 
    float vz, phi, momentum, chi2pid; 
    bool found_trig_electron = false; 

    int event_count = 0;

	TH2F *h_vz_vs_phi = new TH2F("h_vz_vs_phi", "Vz vs phi; Momentum p (GeV); beta", 200, 0, 10, 200, 0.8, 1.05);  // Y-axis range for Epcal
	
	
	

     

    while (chainLD.Next()) {

    //if (event_count >= 10000) break; 
    if (event_count >= 10000000) break;         
    auto particles = c12LD->getDetParticles();
        for (auto& p : particles) { 
            // Check if the particle is an electron
            if (p->par()->getPid() == 11) { 
                int status = p->par()->getStatus(); 
                vz = p->par()->getVz(); 
                chi2pid = p->par()->getVz(); 
                momentum = p->par()->getP(); 
			//Apply forward detector status cut only 
			if ((abs(status)/2000)==1 && chi2pid > -5 && chi2pid < 5 && vz >= -20 && vz <= 5) {
                found_trig_electron = true; 
                break; 

            }				
        }

        if (found_trig_electron){
            for (auto& p : particles) { 
                int status = p->par()->getStatus(); 
                float mom = p->par()->getP(); 
                float beta = p->par()->getBeta(); 
                if ((abs(status)/2000) != 1) continue;
                //cout << "mom: " << mom << " and " << "beta: " << beta; 
                if (p->par()->getPid() == 211) { 
                    h_vz_vs_phi->Fill(mom,beta); 

                }


            }


        }
	}	
    event_count++; 
    }
	TCanvas* canvas_betaVsPP = new TCanvas("canvas_betaVsPP", "beta Vs Momentum for EB 211", 1200, 800);
    h_vz_vs_phi->Draw("COLZ");
    gPad->SetLogz(); 
	canvas_betaVsPP->SaveAs("beta_vs_p_EBpion_Forward.pdf");
 
}
