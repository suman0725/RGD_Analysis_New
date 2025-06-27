#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include <filesystem>
#include <algorithm>
#include "clas12reader.h"
#include "TCanvas.h"
#include "TH2F.h"

namespace fs = std::filesystem;

using namespace clas12;
using namespace std;

typedef std::map<int, std::vector<int>> IndexMap;

// Function to load map by index
IndexMap loadMapByIndex(hipo::bank& fromBank, const char* idxVarName) {
    IndexMap map;
    if (fromBank.getRows() > 0) {
        for (int iFrom = 0; iFrom < fromBank.getRows(); ++iFrom) {
            int iTo = fromBank.getInt(idxVarName, iFrom);
            if (map.find(iTo) == map.end()) {
                map[iTo] = std::vector<int>();
            }
            map[iTo].push_back(iFrom);
        }
    }
    return map;
}

// Function to add .hipo files from directories to vector
void getHipoFiles(std::vector<std::string>& filePaths, const std::vector<std::string>& baseDirs) {
    // Iterate through all provided base directories
    for (const auto& baseDir : baseDirs) {
        for (const auto& entry : fs::recursive_directory_iterator(baseDir)) {
            if (entry.path().extension() == ".hipo") {
                filePaths.push_back(entry.path().string());  // Add file path to vector
            }
        }
    }

    // Sort the file paths alphabetically
    std::sort(filePaths.begin(), filePaths.end());
}

void cherenkovAngleVsMomentum_1() {
    // Example directories where your .hipo files are stored
    std::vector<std::string> baseDirs = {
        "/lustre24/expphy/cache/hallb/scratch/rg-d/production/Bspot/v5dstLD2/dst/recon/018433/"
        // Add more directories as needed
    };

    // Store the .hipo file paths
    std::vector<std::string> filePaths;
    getHipoFiles(filePaths, baseDirs);

    // Prepare the reader, event, and bank structure
    hipo::reader reader;
    hipo::dictionary factory;
    reader.readDictionary(factory);
    hipo::event event;

    // Initialize 2D histogram for Cherenkov angle (in mrad) vs momentum
    TH2F* hCherAngleVsMomentum = new TH2F("hCherAngleVsMomentum", "Cherenkov Angle vs Momentum (in mrad)", 100, 0, 5, 100, 0, 1000);

    // Iterate over the files and process events
    int counter = 0;
    for (const auto& filePath : filePaths) {
        // Open the .hipo file
        reader.open(filePath.c_str());  
        std::cout << "Processing file: " << filePath << std::endl;

        // Read events from the file
        while (reader.next()) {
            reader.read(event);
            // Get the particle and RICH banks
            hipo::bank PART(factory.getSchema("REC::Particle"));
            hipo::bank RICH(factory.getSchema("RICH::Ring"));
            event.getStructure(PART);
            event.getStructure(RICH);

            // Load the map for RICH::Ring based on particle index
            IndexMap richMap = loadMapByIndex(RICH, "pindex");

            int nrows = PART.getRows();
            for (int i = 0; i < nrows; i++) {
                int charge = PART.getInt("charge", i); // Particle charge
                if (charge > 0) { // Only positively charged particles
                    float px = PART.getFloat("px", i);
                    float py = PART.getFloat("py", i);
                    float pz = PART.getFloat("pz", i);

                    // Calculate momentum
                    float momentum = std::sqrt(px * px + py * py + pz * pz);

                    // Apply momentum filter (3-8 GeV)
                    if (momentum < 3 || momentum > 8) continue;

                    // Check RICH::Ring data linked to this particle
                    for (int richIndex : richMap[i]) {
                        float etaC = RICH.getFloat("etaC", richIndex); // Cherenkov angle in radians

                        // Convert etaC (radians) to milliradians (mrad)
                        float cherAngleInMrad = etaC * 1000;

                        // Fill the histogram
                        hCherAngleVsMomentum->Fill(momentum, cherAngleInMrad);
                    }
                }
            }
            counter++;
        }
        
    }

    std::cout << "Number of events processed: " << counter << std::endl;

    // Draw the histogram
    TCanvas* c1 = new TCanvas("c1", "Cherenkov Angle vs Momentum (in mrad)");
    hCherAngleVsMomentum->Draw("COLZ");

    // Set log scale for Z-axis using gPad
    gPad->SetLogz();

    // Save the plot
    c1->Print("CherAngleVsMomentum_mrad_LogZ_MomentumFilter.png");
}
