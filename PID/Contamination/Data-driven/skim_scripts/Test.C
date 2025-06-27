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
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "HipoChain.h"

using namespace clas12;

void SetLorentzVector(TLorentzVector &p4, clas12::region_part_ptr rp) {
    p4.SetXYZM(rp->par()->getPx(), rp->par()->getPy(),
               rp->par()->getPz(), p4.M());
}

void Test() {
    gBenchmark->Start("timer");
    int counter = 0;
    int total_electron_count = 0;
    int total_positron_count = 0;
    const int max_event = 10000; // Process up to 10,000 events

    // Initialize HipoChain with the specific file
    clas12root::HipoChain chain;
    //chain.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v11/CuSn/skim_run_019066.hipo");
    chain.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v11/LD2/skim_run_018873.hipo");

    auto config_c12 = chain.GetC12Reader();
    auto& c12 = chain.C12ref();

    // Event loop
    while (chain.Next() && counter < max_event) {
        // Debug output for first 10 events
        if (counter < 10) {
            int particle_count = 0;
            std::cout << "Event " << counter << ": ";
            for (auto& p : c12->getDetParticles()) {
                int pid = p->par()->getPid();
                std::cout << pid << " ";
                particle_count++;
                if (pid == 11) total_electron_count++;
                if (pid == -11) total_positron_count++;
            }
            std::cout << ", " << particle_count << " particles" << std::endl;
        } else {
            for (auto& p : c12->getDetParticles()) {
                int pid = p->par()->getPid();
                if (pid == 11) total_electron_count++;
                if (pid == -11) total_positron_count++;
            }
        }

        counter++;
    }

    // Print results
    auto bcharge = chain.TotalBeamCharge();
    std::cout << "Number of Events = " << counter << " total charge = " << bcharge << std::endl;
    std::cout << "Total EB electrons = " << total_electron_count << std::endl;
    std::cout << "Total EB positrons = " << total_positron_count << std::endl;

    gBenchmark->Stop("timer");
    gBenchmark->Print("timer");
}