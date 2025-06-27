#include <cstdlib>
#include <string>
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

using namespace clas12;

void SetLorentzVector(TLorentzVector &p4, clas12::region_part_ptr rp) {
    double mass = TDatabasePDG::Instance()->GetParticle(rp->par()->getPid())->Mass();
    p4.SetXYZM(rp->par()->getPx(), rp->par()->getPy(), rp->par()->getPz(), mass);
}


void mass() {
    auto db = TDatabasePDG::Instance();
    double beamEnergy = 10.532; // GeV
    TLorentzVector beam(0, 0, beamEnergy, beamEnergy); // Assuming the beam is along the z-axis
    TLorentzVector target(0, 0, 0, db->GetParticle(2212)->Mass()); // Proton target

    clas12root::HipoChain chainLD2; 
    chainLD2.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstLD2/dst/recon/018309/rec_clas_018309.evio.00035-00039.hipo");
    //chainLD2.Add("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/Test/test.hipo");
    chainLD2.db()->turnOffQADB();

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
                cout<< "electron momentum: " << "Px:"<< el_temp.Px() << endl; 
                std::cout << "Electron mass: " << el_temp.M() << " GeV/c^2" << std::endl; // Print the mass of the electron
                foundElectron = true;
                nu = beam.Energy() - el_temp.Energy(); // Calculate nu
                q = beam - el_temp; // Calculate q after finding the electron
                Q2 = -q.Mag2();
                double W2 = target.M() * target.M() - Q2 + 2 * target.M() * nu;
                double xB = Q2 / (2 * target.M() * nu);
                double y = nu / beam.Energy();
            }

            // Pion selection
            if (pid == 211 && chi2Pid > -10 && chi2Pid < 10) {
                SetLorentzVector(pip_temp, p);
                std::cout << "Mass of selected pion: " << pip_temp.M() << " GeV/c^2" << std::endl;
                foundPion = true;
            }
        }        
        if (foundElectron && foundPion) {
            double z = pip_temp.E() / nu; 
            // Additional calculations can be done here
        } 
    }
}
