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
#include "TStyle.h"
#include "TColor.h"

namespace fs = std::filesystem;

using namespace clas12;
using namespace std;

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


void betavsp() {     


    clas12root::HipoChain chainLD;
    /* std::vector<std::string> baseDirectories = {
    "/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstLD2/dst/recon/018431"
    };  */
    chainLD.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstLD2/dst/recon/018433/rec_clas_018433.evio.00035-00039.hipo");
   // addHipoFiles(chainLD, baseDirectories);    
    chainLD.db()->turnOffQADB();
    auto& c12LD = chainLD.C12ref();

    int status; 
    int sector; 
    
    int event_count = 0;
    TH2F *h_beta_vs_p = new TH2F("h_beta_vs_p", "Beta vs Momentum;Momentum (GeV/c);Beta",
                                 100, 0, 5, 100, 0, 1.2);
 
 
    while (chainLD.Next()) {   

        auto particles = c12LD->getDetParticles();
        for (auto& p : particles) { 
            int charge = p->par()->getCharge();
        	int status = p->par()->getStatus();
            double beta = p->par()->getBeta();  
            double mom = p->par()->getP(); 
            int detector = p->sci(FTOF1A)->getDetector(); 
            int layer = p->sci(FTOF1A)->getLayer(); 
            if (charge > 0 && abs(status / 2000) == 1 && detector==12 && layer==1) {
                h_beta_vs_p->Fill(mom, beta);
                                
            }

        }
    }
    TCanvas* canvas = new TCanvas("canvas", "chi2Pid distribution", 800, 600);
    // Define a custom palette (Reversed order: orange to dark blue)
    Int_t nColors = 5; // Number of colors in the gradient
    Double_t stops[]  = {0.00, 0.25, 0.50, 0.75, 1.00}; // Gradient stops
    Double_t red[]    = {0.00, 0.00, 0.00, 1.00, 1.00}; // Red increases at the higher end (orange)
    Double_t green[]  = {0.00, 0.00, 1.00, 0.00, 0.50}; // Green peaks in the middle
    Double_t blue[]   = {0.50, 1.00, 0.00, 0.00, 0.00}; // Blue dominates at the lower end (dark blue)
    Int_t palette = TColor::CreateGradientColorTable(nColors, stops, red, green, blue, 99);
    gStyle->SetPalette(99, palette);

    h_beta_vs_p->SetStats(0); // Hide stats box
    h_beta_vs_p->SetOption("COLZ"); // Color map
    h_beta_vs_p->SetContour(99); // Smooth color gradient
    canvas->SetLogz();
    h_beta_vs_p->Draw("COLZ");

    canvas->SaveAs("beta_vs_p.png");

}
