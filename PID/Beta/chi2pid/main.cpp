#include <iostream>
#include <chrono>
#include <filesystem>
#include <reader.h>
#include <cmath>
using namespace std;
namespace fs = std::filesystem;


const int HTCC_DETECTOR = 15;
const int LTCC_DETECTOR = 18; // Added for LTCC handling
const int ECAL_DETECTOR = 7;
const int FTOF_DETECTOR = 12;
const double m_pi = 0.13957;
const double m_k  = 0.49367;
const float c = 29.9792458; // speed of light in cm/ns

const vector<pair<int, vector<int>>> chargedBetaDetectors = {
    {FTOF_DETECTOR, {2, 1, 3}}, {CTOF_DETECTOR, {1}}, {ECAL_DETECTOR, {1, 4, 7}}
};

// Function to calculate expected beta for a given particle, momentum from Drift Chamber and mass from PDG
double beta(double p, double mass) {
    return p / sqrt(p*p + mass*mass);
}


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


float computeDeltaT(const ParticleData& p, int pid) {
    float delta_t = 99999.0f;
    bool found = false;

    for (const auto& det : chargedBetaDetectors) {
        for (int layer : det.second) {
            if (hasHit(p, det.first, layer)) {
                for (const auto& hit : p.hits) {
                    if (hit.detector == det.first && hit.layer == layer) {
                        float mass = PDG_MASS[pid];
                        float beta_theory = p.p / sqrt(p.p * p.p + mass * mass);
                        if (beta_theory <= 0) return 99999.0f;
                        float vt = hit.time - hit.path / (SPEED_OF_LIGHT * beta_theory);
                        delta_t = vt - p.start_time;
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
            if (found) break;
        }
        if (found) break;
    }
    return delta_t;
}



int main (){


    // Start the timer
    auto start = chrono::high_resolution_clock::now();


    ifstream inputFile("directories.txt");
    if (!inputFile.is_open()) {
        cerr << "Error: Could not open directories.txt" << endl;
        return 1;
    }
    vector<string> directories;
    string dir;
    while (getline(inputFile, dir)) {
        if (!dir.empty()) directories.push_back(dir);
    }
    inputFile.close();
    if (directories.empty()) {
        cerr << "Error: No directories found" << endl;
        return 1;
    }


    // Loop through each directory
    for (const auto& dir : directories) {
        cout << "Processing directory: " << dir << endl;
        vector<string> hipoFiles;
        for (const auto& entry : fs::directory_iterator(dir)) {
            if (entry.is_regular_file() && entry.path().extension() == ".hipo") {
                hipoFiles.push_back(entry.path().string());
            }
        }
        if (hipoFiles.empty()) {
            cout << "No .hipo files found in directory: " << dir << endl;
            continue;
        }

        const int maxEvents = 100000;
        int event_count = 0;
        int total_events_processed = 0;
        int directory_count = 0;
        int hipo_file_count = 0;

        // Process each HIPO file
        for (const auto& file : hipoFiles) {
            cout << "Processing file: " << file << endl;
            // Add your HIPO file processing code here
            hipo::reader reader;
            reader.open(file.c_str());
            hipo::dictionary factory;
            reader.readDictionary(factory);

            while (reader.next() && event_count < maxEvents) {
                event_count++;

                hipo::event event;
                reader.read(event);
                hipo::bank PART(factory.getSchema("REC::Particle"));
                hipo::bank EVENT(factory.getSchema("REC::Event"));
                hipo::bank SCIN(factory.getSchema("REC::Scintillator"));
                hipo::bank CHER(factory.getSchema("REC::Cherenkov"));
                hipo::bank CALO(factory.getSchema("REC::Calorimeter"));
                event.getStructure(PART);
                event.getStructure(EVENT);
                event.getStructure(SCIN);
                event.getStructure(CHER);
                event.getStructure(CALO);
            }
        }
    }




    double beta_pi = beta(p, m_pi);
    double beta_k  = beta(p, m_k);
    double beta_mid = 0.5 * (beta_pi + beta_k);





    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Execution time: " << duration.count() << " milliseconds" << endl;
}