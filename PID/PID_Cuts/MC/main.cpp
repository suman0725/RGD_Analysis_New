#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <cmath>
#include "reader.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TAxis.h"

namespace fs = std::filesystem;
using namespace std;

void matching() {
    // Define momentum bins (0.1–7 GeV, 20 bins, 0.345 GeV width)
    const int nBins = 30;
    const double pMin = 0.1, pMax = 9.0;
    double pBins[nBins + 1];
    for (int i = 0; i <= nBins; ++i) {
        pBins[i] = pMin + i * (pMax - pMin) / nBins;
    }

    // Counters for pions and kaons (accumulated across files)
    vector<int> pionCounts(nBins, 0); // Reconstructed π⁺
    vector<int> kaonCounts(nBins, 0); // π⁺ matched to MC K⁺

    // Directory containing HIPO files
    string hipoDir = "/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/DATA/MC/temp_filtered";
    vector<string> hipoFiles;
    for (const auto& entry : fs::directory_iterator(hipoDir)) {
        if (entry.is_regular_file() && entry.path().extension() == ".hipo") {
            hipoFiles.push_back(entry.path().string());
        }
    }

    if (hipoFiles.empty()) return;

    // Process each HIPO file
    int totalEvents = 0;
    for (const auto& hipoFile : hipoFiles) {
        hipo::reader reader;
        reader.open(hipoFile.c_str());

        // Read dictionary and initialize banks
        hipo::dictionary factory;
        reader.readDictionary(factory);
        hipo::event event;
        hipo::bank REC_PART(factory.getSchema("REC::Particle"));
        hipo::bank MC_PART(factory.getSchema("MC::Particle"));

        if (!factory.hasSchema("REC::Particle") || !factory.hasSchema("MC::Particle")) continue;

        // Event loop
        while (reader.next()) {
            totalEvents++;
            if (totalEvents % 100000 == 0) {
                cout << "Processed " << totalEvents << " events" << endl;
            }

            reader.read(event);
            event.getStructure(REC_PART);
            event.getStructure(MC_PART);

            int recRows = REC_PART.getRows();
            int mcRows = MC_PART.getRows();

            // Loop over reconstructed particles
            for (int i = 0; i < recRows; ++i) {
                if (REC_PART.getInt("pid", i) != 211) continue; // Select π⁺

                // Get reconstructed kinematics
                float px = REC_PART.getFloat("px", i);
                float py = REC_PART.getFloat("py", i);
                float pz = REC_PART.getFloat("pz", i);
                float pRec = sqrt(px * px + py * py + pz * pz);
                float thetaRec = atan2(sqrt(px * px + py * py), pz) * 180.0 / M_PI;
                float phiRec = atan2(py, px) * 180.0 / M_PI;

                if (std::isnan(pRec) || std::isnan(thetaRec) || std::isnan(phiRec) || pRec <= 0) continue;

                // Find momentum bin
                int bin = -1;
                for (int j = 0; j < nBins; ++j) {
                    if (pRec >= pBins[j] && pRec < pBins[j + 1]) {
                        bin = j;
                        break;
                    }
                }
                if (bin < 0) continue;
                pionCounts[bin]++;

                // Loop over MC particles
                double minDeltaR = 1e6;
                int bestMCIndex = -1;
                for (int j = 0; j < mcRows; ++j) {
                    float mcPx = MC_PART.getFloat("px", j);
                    float mcPy = MC_PART.getFloat("py", j);
                    float mcPz = MC_PART.getFloat("pz", j);
                    float thetaMC = atan2(sqrt(mcPx * mcPx + mcPy * mcPy), mcPz) * 180.0 / M_PI;
                    float phiMC = atan2(mcPy, mcPx) * 180.0 / M_PI;

                    if (std::isnan(thetaMC) || std::isnan(phiMC)) continue;

                    double deltaTheta = std::abs(thetaRec - thetaMC);
                    double deltaPhi = std::abs(phiRec - phiMC);
                    if (deltaPhi > 180.0) deltaPhi = 360.0 - deltaPhi;

                    if (deltaTheta < 1.0 && deltaPhi < 3.0) {
                        double deltaR = std::sqrt(deltaTheta * deltaTheta + deltaPhi * deltaPhi);
                        if (deltaR < minDeltaR) {
                            minDeltaR = deltaR;
                            bestMCIndex = j;
                        }
                    }
                }

                // Check if matched to K⁺ (PID = 321)
                if (bestMCIndex >= 0 && MC_PART.getInt("pid", bestMCIndex) == 321) {
                    kaonCounts[bin]++;
                }
            }
        }
    }

    // Prepare data for TGraph
    vector<double> pMidPoints(nBins);
    vector<double> contamination(nBins);
    for (int i = 0; i < nBins; ++i) {
        pMidPoints[i] = (pBins[i] + pBins[i + 1]) / 2.0;
        contamination[i] = pionCounts[i] > 0 ? static_cast<double>(kaonCounts[i]) / pionCounts[i] : 0.0;
    }

    // Save results to text file
    ofstream outFile("kaon_contamination.txt");
    if (outFile.is_open()) {
        outFile << fixed << setprecision(4);
        outFile << "# K⁺ contamination in reconstructed π⁺: (Reconstructed π⁺ matched to MC K⁺) / (Total reconstructed π⁺)\n";
        outFile << "MomentumBin(GeV)\tContamination\tPionCount\tKaonCount\n";
        for (int i = 0; i < nBins; ++i) {
            outFile << pMidPoints[i] << "\t" << contamination[i] << "\t" << pionCounts[i] << "\t" << kaonCounts[i] << "\n";
        }
        outFile.close();
    }

    // Create TGraph for PDF output
    TGraph* gCont = new TGraph(nBins, pMidPoints.data(), contamination.data());
    gCont->SetTitle("K^{+} Contamination in #pi^{+};p (GeV/c);Contamination");
    gCont->SetMarkerStyle(20); // Solid circles
    gCont->SetMarkerColor(kBlack); // Black points

    TCanvas c1("c1", "Kaon Contamination", 800, 600);
    gStyle->SetGridColor(kGray); // Set grid color to gray
    c1.SetGridx(); // Enable x grid
    c1.SetGridy(); // Enable y grid
    gCont->Draw("AP"); // Draw points only
    c1.SaveAs("kaon_contamination.pdf");
    delete gCont;
}

int main() {
    matching();
    return 0;
}