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
int getSector(double phi) {
    // Assuming phi is already in the expected range of [-180, 180]
    if (phi >= -30 && phi < 30) return 1;
    if (phi >= 30 && phi < 90) return 2;
    if (phi >= 90 && phi < 150) return 3;
    if ((phi >= 150 && phi <= 180) || (phi >= -180 && phi < -150)) return 4;
    if (phi >= -150 && phi < -90) return 5;
    if (phi >= -90 && phi < -30) return 6;

    return -1; // Default if not in any sector
}

void pid() {     

// Process LD Target Files
    clas12root::HipoChain chainLD;
    chainLD.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/prod/v4ob_aideLD2/dst/recon/018528/rec_clas_018528.evio.00310-00314.hipo");
    chainLD.db()->turnOffQADB();
    auto& c12LD = chainLD.C12ref();

    int status; 
    int event_count = 0;
    std::vector<TH1F*> h_nPhe_negative(6);
    for (int i = 0; i < 6; ++i) {
        h_nPhe_negative[i] = new TH1F(Form("h_nPhe_negative_sector%d", i+1), 
                                      Form("Number of Photoelectrons for Sector %d", i+1), 
                                      100, 0, 50);  // Adjust bins and range as needed
    }

    while (chainLD.Next()) {        
    auto particles = c12LD->getDetParticles();
        for (auto& p : particles) { 
            // Check if the particle is negatively charged
            if (p->par()->getCharge() == -1) { 
		int status = p->par()->getStatus(); 

		//Apply forward detector status cut only 
		if (status > -4000 && status <= 2000) {
                // Get the number of photoelectrons from HTCC
                int nPhe = p->che(HTCC)->getNphe(); // Use HTCC to get the number of photoelectrons
		double phi = p->getPhi()* TMath::RadToDeg(); 
		int sector = getSector(phi);
		cout << "Phi: " << phi << " Sector: " << sector << endl;
			if (sector > 0 && sector <= 6) {
	                    h_nPhe_negative[sector - 1]->Fill(nPhe);
	                }
		}
                // Get the status of the particle
                status = p->par()->getStatus();
                
                // Apply forward detector cut
                if (status > -4000 && status <= 2000) {
                    // Particle is measured in the forward detector
                    // Do something with the valid particle, if needed
                }
            }
        }

        // Example of printing the number of negatively charged particles found
       // cout << "Processed an event." << endl;
    }

     // Create a canvas and divide it into six pads (2 rows x 3 columns)
    TCanvas* c = new TCanvas("c", "Number of Photoelectrons by Sector", 1200, 800);
    c->Divide(3, 2);  // Divide canvas into 6 parts (2 rows, 3 columns)


     for (int i = 0; i < 6; ++i) {
        c->cd(i+1);  // Move to the i-th pad
        h_nPhe_negative[i]->Draw();
	// Set log scale on y-axis
        gPad->SetLogy();

        // Set y-axis range (10^-6 to 1e6)
        h_nPhe_negative[i]->GetYaxis()->SetRangeUser(1e0, 1e6);

        // Customize marker style
        h_nPhe_negative[i]->SetMarkerStyle(20);  // Dots
        h_nPhe_negative[i]->SetMarkerSize(1);    // Adjust size as needed
        h_nPhe_negative[i]->Draw("P");           // "P" option for points

        // Draw a vertical red line at nPhe = 2
        TLine* line = new TLine(2, 1e0, 2, 1e6);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);  // Adjust line width as needed
        line->Draw("same");
    }

    // Save the canvas as a PNG file
    c->SaveAs("Photoelectrons_negative_particle_LD2_all_sectors_FD.png");

}
