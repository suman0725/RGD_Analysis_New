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

int determineSector(double phi) {
    phi *= TMath::RadToDeg(); // Convert phi from radians to degrees
    if (phi >= -30 && phi < 30) return 0;
    else if (phi >= 30 && phi < 90) return 1;
    else if (phi >= 90 && phi < 150) return 2;
    else if (phi >= 150 || phi < -150) return 3;
    else if (phi >= -150 && phi < -90) return 4;
    else if (phi >= -90 && phi < -30) return 5;
    return -1; // Should not happen
}

void vertex3() {
    auto start = std::chrono::high_resolution_clock::now();

    clas12root::HipoChain chain;
    std::string baseDirectory = "/lustre19/expphy/volatile/clas12/rg-d/production/prod/v3_ob_CuSn/dst/recon";
    addHipoFiles(chain, baseDirectory);

    TFile *outFile = new TFile("vertex3.root", "RECREATE");
    TTree *tree = new TTree("ElectronData", "Complete Electron Data Tree");

    // Variables to be stored
    float vz, vx, vy, px, py, pz, charge,  phi, status;
    int sector;

    // Branches for the tree
    tree->Branch("vz", &vz, "vz/F");
    tree->Branch("vx", &vx, "vx/F");
    tree->Branch("vy", &vy, "vy/F");
    tree->Branch("phi", &phi, "phi/F");
    tree->Branch("status", &status, "status/I");
    tree->Branch("sector", &sector, "sector/I");

    chain.db()->turnOffQADB();
    auto& c12 = chain.C12ref();
    int c=0;
    while (chain.Next()) {
            //if(c>1000000) break; 
        auto electrons = c12->getByID(11);
        for (auto& electron : electrons) {
            phi = electron->getPhi()TMath::RadToDeg();
            sector = determineSector(phi); // Determine the sector based on phi
            if (sector != -1) {
                vz = electron->par()->getVz();
                vx = electron->par()->getVx();
                vy = electron->par()->getVy();
                status = electron->par()->getStatus();
                //int sectag=sector;
                tree->Fill();
             c++;
            }
	}
    }

    outFile->Write();
    outFile->Close();
    delete outFile; 
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds." << std::endl;
}
