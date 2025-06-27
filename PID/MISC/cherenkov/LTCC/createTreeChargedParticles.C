
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
#include <filesystem>

namespace fs = std::filesystem;
using namespace clas12;
using namespace std;

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
void addHipoFiles(clas12root::HipoChain& chain, const std::vector<std::string>& baseDirs) {
    // A vector to store the paths of the .hipo files
    std::vector<fs::path> files;

    // Iterate through all provided base directories
    for (const auto& baseDir : baseDirs) {
        for (const auto& entry : fs::recursive_directory_iterator(baseDir)) {
            if (entry.path().extension() == ".hipo") {
                files.push_back(entry.path());  // Add to the list
            }
        }
    }

    // Sort the files by their paths (alphabetically)
    std::sort(files.begin(), files.end());

    // Add each file to the HipoChain
    for (const auto& file : files) {
        chain.Add(file.string());  // Add the file to the chain
        std::cout << "Added file: " << file << std::endl;
    }
}



void createTreeChargedParticles() {
    clas12root::HipoChain chainLD;
    std::vector<std::string> baseDirectories = {
        "/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstLD2/dst/recon/018431", //
        "/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstLD2/dst/recon/018433",
        "/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstLD2/dst/recon/018528",
        "/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstLD2/dst/recon/018644",
        "/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstLD2/dst/recon/018851",
        "/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstLD2/dst/recon/018873"
    };
    addHipoFiles(chainLD, baseDirectories);    
    chainLD.db()->turnOffQADB();
    auto& c12LD = chainLD.C12ref();

    // Create a ROOT file and TTree
    TFile* outFile = new TFile("charged_particles_LTCC_HTCC.root", "RECREATE");
    TTree* tree = new TTree("charged_particle", "Tree for charged particles");

    // Branch variables
    int pid;
    float px, py, pz, vx, vy, vz, vt, beta, chi2pid, startTime;
    int charge, detector, sector;
    float nphe, time, path, x, y, z, dtheta, dphi;
    short status;

    // Create branches
    tree->Branch("pid", &pid, "pid/I");
    tree->Branch("px", &px, "px/F");
    tree->Branch("py", &py, "py/F");
    tree->Branch("pz", &pz, "pz/F");
    tree->Branch("vx", &vx, "vx/F");
    tree->Branch("vy", &vy, "vy/F");
    tree->Branch("vz", &vz, "vz/F");
    tree->Branch("vt", &vt, "vt/F");
    tree->Branch("charge", &charge, "charge/I");
    tree->Branch("beta", &beta, "beta/F");
    tree->Branch("chi2pid", &chi2pid, "chi2pid/F");
    tree->Branch("status", &status, "status/S");

    // Additional branches for Cherenkov variables
    tree->Branch("detector", &detector, "detector/I");
    tree->Branch("sector", &sector, "sector/I");
    tree->Branch("nphe", &nphe, "nphe/F");
    tree->Branch("time", &time, "time/F");
    tree->Branch("path", &path, "path/F");
    tree->Branch("x", &x, "x/F");
    tree->Branch("y", &y, "y/F");
    tree->Branch("z", &z, "z/F");
    tree->Branch("dtheta", &dtheta, "dtheta/F");
    tree->Branch("dphi", &dphi, "dphi/F");

   
    int counter = 0; 
    // Loop through events
    while (chainLD.Next()) {  

       // if (counter>=100000) break; 

        auto particles = c12LD->getDetParticles();
        for (auto& p : particles) { 
            charge = p->par()->getCharge();
            if (charge != 0) { // Process only charged particles
                pid = p->par()->getPid();
                px = p->par()->getPx();
                py = p->par()->getPy();
                pz = p->par()->getPz();
                vx = p->par()->getVx();
                vy = p->par()->getVy();
                vz = p->par()->getVz();
                vt = p->par()->getVt();
                beta = p->par()->getBeta();
                chi2pid = p->par()->getChi2Pid();
                status = p->par()->getStatus();

                // Check for LTCC
                if (p->che(LTCC)) {
                    detector = p->che(LTCC)->getDetector();
                    sector = p->che(LTCC)->getSector();
                    nphe = p->che(LTCC)->getNphe();
                    time = p->che(LTCC)->getTime();
                    path = p->che(LTCC)->getPath();
                    x = p->che(LTCC)->getX();
                    y = p->che(LTCC)->getY();
                    z = p->che(LTCC)->getZ();
                    dtheta = p->che(LTCC)->getDtheta();
                    dphi = p->che(LTCC)->getDPhi();
                    
                    // Fill the tree
                    tree->Fill();
                } 

                // Check for HTCC
                if (p->che(HTCC)) {
                    detector = p->che(HTCC)->getDetector();
                    sector = p->che(HTCC)->getSector();
                    nphe = p->che(HTCC)->getNphe();
                    time = p->che(HTCC)->getTime();
                    path = p->che(HTCC)->getPath();
                    x = p->che(HTCC)->getX();
                    y = p->che(HTCC)->getY();
                    z = p->che(HTCC)->getZ();
                    dtheta = p->che(HTCC)->getDtheta();
                    dphi = p->che(HTCC)->getDPhi();
                
                    // Fill the tree
                    tree->Fill();
                } 
            }
        }
    counter ++;
    }

    // Write and close the file
    outFile->Write();
    outFile->Close();

    std::cout << "Tree with charged particles and Cherenkov data saved to charged_particles.root" << std::endl;
}
