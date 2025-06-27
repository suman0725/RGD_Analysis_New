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

using namespace clas12;
using namespace std;

// Function to set the Lorentz vector for a given particle
void SetLorentzVector(TLorentzVector &p4, clas12::region_part_ptr rp) {
    p4.SetXYZM(rp->par()->getPx(), rp->par()->getPy(), rp->par()->getPz(), p4.M());
}

void multiplicityRatio2_new_new_new() {
    auto db = TDatabasePDG::Instance();
    int nBins = 7;
    double zMin = 0.1;
    double zMax = 0.8;
    const double m_e = db->GetParticle(11)->Mass(); // Mass of electron
    const double m_pion = db->GetParticle(211)->Mass(); // Mass of pion

    TLorentzVector beam(0, 0, 10.532, 10.532); // Beam 4-vector (assuming 10.532 GeV)
    TLorentzVector target(0, 0, 0, db->GetParticle(2212)->Mass()); // Target 4-vector (assuming proton target)

    // Histograms for LD and Nuclear targets
    auto* h_LD_hadrons = new TH1F("h_LD_hadrons", "LD Hadrons;z;Counts", nBins, zMin, zMax);
    auto* h_Nuc_hadrons = new TH1F("h_Nuc_hadrons", "Nuclear Target Hadrons;z;Counts", nBins, zMin, zMax);

    // Process LD Target Files
    clas12root::HipoChain chainLD;
    chainLD.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/prod/v4ob_aideLD2/dst/recon/018528/rec_clas_018528.evio.00310-00314.hipo");
    chainLD.db()->turnOffQADB();
    auto& c12LD = chainLD.C12ref();

    double total_LD_electrons = 0;

    while (chainLD.Next()) {
        auto particles = c12LD->getDetParticles();
        double nu = 0.0;
        double max_E = 0.0;
        TLorentzVector selectedElectron;
        bool eventPassesCuts = false;
        bool hasElectron = false;
        bool hasHadron = false;

        // Find the leading electron (highest energy)
        for (auto& p : particles) {
            int pid = p->par()->getPid();

            if (pid == 11) { // Electron
                SetLorentzVector(selectedElectron, p);
                double E = selectedElectron.E();
                if (E > max_E) {
                    max_E = E;
                    nu = beam.E() - E; // Update nu for the selected electron
                    hasElectron = true; // Mark that we have found an electron
                }
            } else if (abs(pid) == 211) { // Hadron (e.g., pion)
                hasHadron = true; // Mark that we have found a hadron
            }
        }

        if (hasElectron && hasHadron && max_E > 0) {
            // Apply cuts based on Q², W, y
            TLorentzVector q = beam - selectedElectron;
            double Q2 = -q.Mag2();
            double W = sqrt(target.M() * target.M() + 2 * target.M() * nu - Q2);
            double y = nu / beam.E();

            if (Q2 > 1 && y > 0.25 && y < 0.85 && W > 2) {
                eventPassesCuts = true; // Mark event as passing the cuts
            }
        }

        // If an electron was selected and cuts were passed, process hadrons
        if (eventPassesCuts) {
            total_LD_electrons += 1; // Count selected electron
            for (auto& p : particles) {
                int pid = p->par()->getPid();
                if (abs(pid) == 211) { // Hadron (e.g., pions)
                    TLorentzVector pion;
                    pion.SetM(m_pion); // Set mass to pion mass
                    SetLorentzVector(pion, p);
                    double z = pion.E() / nu;
                    h_LD_hadrons->Fill(z);
                }
            }
        }
    }

    // Process Nuclear Target Files
    clas12root::HipoChain chainNuc;
    chainNuc.Add("/lustre24/expphy/cache/hallb/scratch/rg-d/production/prod/v4ob_aideCxC/dst/recon/018451/rec_clas_018451.evio.00310-00314.hipo");
    chainNuc.db()->turnOffQADB();
    auto& c12Nuc = chainNuc.C12ref();

    double total_Nuc_electrons = 0;

    while (chainNuc.Next()) {
        auto particles = c12Nuc->getDetParticles();
        double nu_c = 0.0;
        double max_E = 0.0;
        TLorentzVector selectedElectron;
        bool eventPassesCuts = false;
        bool hasElectron = false;
        bool hasHadron = false;

        // Find the leading electron (highest energy)
        for (auto& p : particles) {
            int pid = p->par()->getPid();

            if (pid == 11) { // Electron
                SetLorentzVector(selectedElectron, p);
                double E = selectedElectron.E();
                if (E > max_E) {
                    max_E = E;
                    nu_c = beam.E() - E; // Update nu_c for the selected electron
                    hasElectron = true; // Mark that we have found an electron
                }
            } else if (abs(pid) == 211) { // Hadron (e.g., pion)
                hasHadron = true; // Mark that we have found a hadron
            }
        }

        if (hasElectron && hasHadron && max_E > 0) {
            // Apply cuts based on Q², W, y
            TLorentzVector q = beam - selectedElectron;
            double Q2 = -q.Mag2();
            double W = sqrt(target.M() * target.M() + 2 * target.M() * nu_c - Q2);
            double y = nu_c / beam.E();

            if (Q2 > 1 && y > 0.25 && y < 0.85 && W > 2) {
                eventPassesCuts = true; // Mark event as passing the cuts
            }
        }

        // If an electron was selected and cuts were passed, process hadrons
        if (eventPassesCuts) {
            total_Nuc_electrons += 1; // Count selected electron
            for (auto& p : particles) {
                int pid = p->par()->getPid();
                if (abs(pid) == 211) { // Hadron (e.g., pions)
                    TLorentzVector pion;
                    pion.SetM(m_pion); // Set mass to pion mass
                    SetLorentzVector(pion, p);
                    double z = pion.E() / nu_c;
                    h_Nuc_hadrons->Fill(z);
                }
            }
        }
    }

    // Apply double ratio calculation
    TGraph* graph = new TGraph();
    int pointIndex = 0;

    for (int bin = 1; bin <= nBins; ++bin) {
        double binCenter = h_LD_hadrons->GetBinCenter(bin);
        double count_LD_hadrons = h_LD_hadrons->GetBinContent(bin);
        double count_Nuc_hadrons = h_Nuc_hadrons->GetBinContent(bin);

        double ratio = 0.0;
        if (total_LD_electrons > 0 && total_Nuc_electrons > 0) {
            ratio = (count_Nuc_hadrons / total_Nuc_electrons) / (count_LD_hadrons / total_LD_electrons);
        }

        graph->SetPoint(pointIndex++, binCenter, ratio);
    }

    // Save histograms and scatter plot
    TFile* outFile = new TFile("multiplicity_ratio_z.root", "RECREATE");
    h_LD_hadrons->Write();
    h_Nuc_hadrons->Write();
    graph->Write();
    outFile->Close();

    // Optionally display scatter plot in a canvas
    TCanvas* canvas = new TCanvas("canvas", "Multiplicity Ratio vs z", 800, 600);
    graph->SetTitle("Multiplicity Ratio vs z;z;R");
    graph->SetMarkerStyle(kFullCircle);
    graph->Draw("AP");

    canvas->SaveAs("multiplicity_ratio_vs_z_4.png");
}
