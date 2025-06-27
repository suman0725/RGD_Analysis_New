#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include "clas12reader.h"
#include <TApplication.h>
#include <iostream>
#include <cstdlib>
#include <filesystem>
#include <TChain.h>
#include <chrono>
#include "HipoChain.h"
namespace fs = std::filesystem;

using namespace clas12;
using namespace std; 

void readEvents() {

    auto start = std::chrono::high_resolution_clock::now();


    clas12root::HipoChain chain;
   
    chain.Add("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/Test/skim_run_018329.hipo");
    /* TFile *outFile = new TFile("output_histograms4.root", "RECREATE");
    TTree *tree = new TTree("ElectronData", "Complete Electron Data Tree");
 */
    // Variables to be stored
    //float vz, vx, vy, phi, status;

    // Branches for the tree
    /* tree->Branch("vz", &vz, "vz/F");
    tree->Branch("vx", &vx, "vx/F");
    tree->Branch("vy", &vy, "vy/F");
    tree->Branch("status", &status, "status/I");
    tree->Branch("phi", &phi, "phi/F"); */
     
                
    chain.db()->turnOffQADB();
    auto& c12 = chain.C12ref();
    int counter = 0; 

    while (chain.Next()) {
        if (counter == 0 ) break; 
        auto particles = c12->getDetParticles(); 
        for (auto& p : particles) {
            int charge = p->par()->getCharge(); 
            cout << "charge is = " << charge; 
            if (charge <= 0 ) continue;
           /*  float mom = p->par()->getP(); 
            float protonbeta = mom / sqrt(mom * mom + 0.93827 * 0.93827); */
            float path = p->sci(FTOF1B)->getPath(); 
            float time = p->sci(FTOF1B)->getTime(); 
            float beta = path / time;

            cout << "Beta from SCIN = " << beta << endl; 

            //tree->Fill();
        }
        counter++; 
    }

    /* outFile->Write();
    outFile->Close(); */

    

   auto end = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = end - start;
   std::cout << "Elapsed time: " << elapsed.count() << " seconds." << std::endl;
}
