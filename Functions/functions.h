#include <cstdlib>
#include <iostream>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <cmath>




// calculating phi_h
double CalculatePhih(TLorentzVector q, TLorentzVector p, TLorentzVector p_h) {
    // Extract 3-vectors from the Lorentz vectors
    TVector3 q3 = q.Vect();   // q vector
    TVector3 p3 = p.Vect();   // incident electron
    TVector3 p_h3 = p_h.Vect(); // hadron's vector

    // Normalize the q vector
   TVector3 VectorQNorm = q3.Unit();

    // Compute the cross products
    TVector3 cross_qp = VectorQNorm.Cross(p3);      // (q x p)
    TVector3 cross_qph = VectorQNorm.Cross(p_h3);   // (q x p_h)
    TVector3 cross_pph = p3.Cross(p_h3);
   
    // Calculate the angle in degrees
    // Calculate the angle according to the Trento convention
    double phi = TMath::ATan2(cross_pph.Dot(VectorQNorm), cross_qp.Dot(cross_qph)) * TMath::RadToDeg();
    phi = phi + 180.0;

    return phi; 
}

double CalculatePt2(TLorentzVector P_h, TLorentzVector q) {
    TVector3 q3 = q.Vect(); 
    TVector3 p_h3 = P_h.Vect(); 
    TVector3 cross_qph = q3.Cross(p_h3);
    double magnitudeCrossProductSquared = cross_qph.Mag2(); // |q × p_h|^2
    double magnitudeQSquared = q3.Mag2();                     // |q|^2
    return magnitudeCrossProductSquared / magnitudeQSquared; // Pt^2
}


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

// nu 
double nu = incidentelectron.Energy() - scatteredelectron.Energy(); 

TLorentzVector q = beam - el_temp; 

double Q2 = -q.Mag2();

double xB = Q2 / (2 * target.M() * nu); // Calculate Bjorken x

double W2 = target.M()*target.M() + 2*target.M()*nu - Q2; // Calculate W^2

double y = nu / beamEnergy;



double z = pip. Energy()/beamEnergy; // calculate z for the pi+


// Function to calculate expected beta for a given particle, momentum from Drift Chamber and mass from PDG
double beta(double p, double mass) {
    return p / sqrt(p*p + mass*mass);
}

//reverse indexing 
IndexMap loadMapByIndex(hipo::bank& fromBank, const char* idxVarName) {
    IndexMap map;
    if (fromBank.getRows() > 0) {
        for (int iFrom = 0; iFrom < fromBank.getRows(); ++iFrom) {
            int iTo = fromBank.getInt(idxVarName, iFrom);
            map[iTo].push_back(iFrom);
        }
    }
    return map;
}


double phi = p->getPhi()* TMath::RadToDeg();