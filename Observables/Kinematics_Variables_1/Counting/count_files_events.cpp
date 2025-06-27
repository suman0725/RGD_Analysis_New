#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <chrono>
#include "reader.h"

using namespace std;
using namespace std::chrono;
namespace fs = std::filesystem;

// Load HIPO files for a target directory
vector<string> loadHipoFiles(const string& baseDir, const string& targetDir) {
    vector<string> files;
    string dir = baseDir + "/" + targetDir;
    if (!fs::exists(dir)) {
        cerr << "Directory " << dir << " does not exist!" << endl;
        return files;
    }
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (entry.path().extension() == ".hipo") {
            files.push_back(entry.path().string());
        }
    }
    return files;
}

// Check configuration and extract run number (IN: torus = -1.0, OB: torus = 1.0)
pair<bool, int> isConfigMatch(const string& file, const string& config) {
    hipo::reader reader;
    reader.open(file.c_str());
    hipo::dictionary factory;
    reader.readDictionary(factory);
    if (!factory.hasSchema("RUN::config")) {
        cerr << "Error: No RUN::config bank in " << file << endl;
        return {false, 0};
    }
    hipo::bank RUN(factory.getSchema("RUN::config"));
    if (reader.gotoEvent(0)) {
        hipo::event event;
        reader.read(event);
        event.getStructure(RUN);
        if (RUN.getRows() == 0) {
            cerr << "Error: Empty RUN::config bank in " << file << endl;
            return {false, 0};
        }
        float torus = RUN.getFloat("torus", 0);
        int run = RUN.getInt("run", 0);
        bool match = (config == "IN" && torus == -1.0) || (config == "OB" && torus == 1.0);
        return {match, run};
    }
    cerr << "Error: Failed to read event 0 in " << file << endl;
    return {false, 0};
}

// Count events in a single HIPO file and return processing time
pair<long, double> countEventsInFile(const string& file) {
    auto start = high_resolution_clock::now();
    hipo::reader reader;
    reader.open(file.c_str());
    long eventCount = 0;
    while (reader.next()) {
        eventCount++;
    }
    auto end = high_resolution_clock::now();
    double duration = duration_cast<microseconds>(end - start).count() / 1e6; // Seconds
    return {eventCount, duration};
}

// Process a target directory and count files and events
void processTarget(const string& baseDir, const string& targetDir, const vector<string>& logicalTargets,
                  map<string, long>& totalFiles, map<string, map<string, long>>& configFiles,
                  map<string, map<string, vector<pair<string, long>>>>& eventsPerFile,
                  map<string, long>& totalEvents, map<string, map<string, set<int>>>& runNumbers,
                  double& totalTime, long& processedFiles, long totalFilesCount) {
    auto files = loadHipoFiles(baseDir, targetDir);
    vector<string> configs = {"IN", "OB"};

    // Map targetDir to logical targets
    vector<string> targets;
    if (targetDir == "CxC") {
        targets = {"C"};
    } else if (targetDir == "LD2") {
        targets = {"LD2"};
    } else if (targetDir == "CuSn") {
        targets = {"Cu", "Sn"};
    }

    // Initialize maps for this targetDir
    for (const auto& target : targets) {
        totalFiles[target] += files.size();
        for (const auto& config : configs) {
            configFiles[target][config] = 0;
            eventsPerFile[target][config] = {};
            runNumbers[target][config] = {};
        }
        totalEvents[target] = 0;
    }

    // Count files, events, and run numbers per configuration
    for (const auto& file : files) {
        for (const auto& config : configs) {
            auto [match, run] = isConfigMatch(file, config);
            if (match) {
                for (const auto& target : targets) {
                    configFiles[target][config]++;
                    auto [eventCount, duration] = countEventsInFile(file);
                    eventsPerFile[target][config].emplace_back(file, eventCount);
                    totalEvents[target] += eventCount;
                    runNumbers[target][config].insert(run);
                    totalTime += duration;
                    processedFiles++;
                }
            }
        }
        // Estimate remaining time
        if (processedFiles > 0 && processedFiles < totalFilesCount) {
            double avgTimePerFile = totalTime / processedFiles;
            long remainingFiles = totalFilesCount - processedFiles;
            double remainingTime = avgTimePerFile * remainingFiles;
            int minutes = static_cast<int>(remainingTime / 60);
            int seconds = static_cast<int>(remainingTime) % 60;
            cout << "Processed " << processedFiles << "/" << totalFilesCount << " files. "
                 << "Estimated time remaining: " << minutes << "m " << seconds << "s" << endl;
        }
    }
}

// Main function
int main() {
    string baseDir = "/lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v11";
    vector<string> targetDirs = {"LD2", "CxC", "CuSn"};
    vector<string> logicalTargets = {"LD2", "C", "Cu", "Sn"};
    vector<string> configs = {"IN", "OB"};

    // Maps to store counts
    map<string, long> totalFiles; // Total files per target
    map<string, map<string, long>> configFiles; // Files per config per target
    map<string, map<string, vector<pair<string, long>>>> eventsPerFile; // Events per file per config per target
    map<string, long> totalEvents; // Total events per target
    map<string, map<string, set<int>>> runNumbers; // Run numbers per config per target

    // Initialize maps
    for (const auto& target : logicalTargets) {
        totalFiles[target] = 0;
        for (const auto& config : configs) {
            configFiles[target][config] = 0;
            eventsPerFile[target][config] = {};
            runNumbers[target][config] = {};
        }
        totalEvents[target] = 0;
    }

    // Count total files for time estimation
    long totalFilesCount = 0;
    for (const auto& dir : targetDirs) {
        totalFilesCount += loadHipoFiles(baseDir, dir).size();
    }

    // Process each target directory
    double totalTime = 0.0;
    long processedFiles = 0;
    for (const auto& dir : targetDirs) {
        processTarget(baseDir, dir, logicalTargets, totalFiles, configFiles, eventsPerFile, totalEvents,
                      runNumbers, totalTime, processedFiles, totalFilesCount);
    }

    // Prepare output
    ostringstream output;
    output << "=== File and Event Counts ===\n\n";
    for (const auto& target : logicalTargets) {
        output << "Target: " << target << "\n";
        output << "  Total Files: " << totalFiles[target] << "\n";
        for (const auto& config : configs) {
            output << "  " << config << " Files: " << configFiles[target][config] << "\n";
            output << "  " << config << " Run Numbers: ";
            if (runNumbers[target][config].empty()) {
                output << "None";
            } else {
                bool first = true;
                for (const auto& run : runNumbers[target][config]) {
                    if (!first) output << ", ";
                    output << run;
                    first = false;
                }
            }
            output << "\n";
            output << "  " << config << " Events Per File:\n";
            for (const auto& [file, count] : eventsPerFile[target][config]) {
                output << "    " << file << ": " << count << " events\n";
            }
        }
        output << "  Total Events: " << totalEvents[target] << "\n\n";
    }

    // Print to console
    cout << output.str();

    // Save to file
    ofstream outFile("counts.txt");
    if (outFile.is_open()) {
        outFile << output.str();
        outFile.close();
        cout << "Results saved to counts.txt" << endl;
    } else {
        cerr << "Error: Could not open counts.txt for writing!" << endl;
    }

    return 0;
}