#include <TFile.h>
#include <TH1F.h>
#include <TSystem.h>
#include <ROOT/RDataFrame.hxx>
#include <iostream>
#include <vector>

using namespace std;
using namespace ROOT;

// chi2 PIC cut functions (momentum-dependent)
double getdtCutNeg(double p) {
    return -0.232 + (-6.405) * exp(-8.123 * p); 
}

double getdtCutPos(double p) {
    return 0.225 + 5.245 * exp(-7.298 * p); 
}

// Calculate theoretical beta (for reference, not used in saving)
double getTheoreticalBeta(double p, double mass) {
    return p / sqrt(p * p + mass * mass);
}

int main() {
    ROOT::EnableImplicitMT(16);

    TFile* file = new TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/Skim/pkptreeCxC_9_test_modified.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file pkptreeCxC_9_test_modified.root\n";
        ROOT::DisableImplicitMT();
        return 1;
    }
    ROOT::RDataFrame df_pions("EB_pos_pion_assumed", file);

    vector<double> pBins;
    double pMin = 1.0, pMax = 7.0;
    int nBins = static_cast<int>((pMax - pMin) / 0.3);
    for (int i = 0; i <= nBins; i++) {
        pBins.push_back(pMin + i * 0.3);
    }

    TFile* outfile = new TFile("beta_histograms.root", "RECREATE");

    for (int i = 0; i < nBins; i++) {
        double pLow = pBins[i], pHigh = pBins[i + 1];

        // Unfiltered (before cut) data
        auto df_pions_before = df_pions.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"});
        ROOT::RDF::TH1DModel modelBetaBefore(
            TString::Format("beta_before_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #beta; Counts", pLow, pHigh),
            100, 0.96, 1.03
        );
        auto histo_before = df_pions_before.Histo1D(modelBetaBefore, "beta");
        histo_before->Write();

        // Filtered (after cut) data
        auto filtered_df_pions = df_pions.Filter([pLow, pHigh](float p) { return p >= pLow && p < pHigh; }, {"p"})
                                            .Filter([](float p, float dt) {
                                                double dtMin = getdtCutNeg(p);
                                                double dtMax = getdtCutPos(p);
                                                return dt > dtMin && dt < dtMax;
                                            }, {"p", "dt"});
        ROOT::RDF::TH1DModel modelBetaAfter(
            TString::Format("beta_after_%d", i),
            TString::Format("p: [%.1f-%.1f) GeV/c; #beta; Counts", pLow, pHigh),
            100, 0.96, 1.03
        );
        auto histo_after = filtered_df_pions.Histo1D(modelBetaAfter, "beta");
        histo_after->Write();
    }

    outfile->Close();
    file->Close();
    delete file;
    delete outfile;
    ROOT::DisableImplicitMT();

    cout << "Histograms saved to beta_histograms.root" << endl;
    return 0;
}