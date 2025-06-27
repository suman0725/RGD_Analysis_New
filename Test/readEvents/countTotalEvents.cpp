#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>  // For automated file detection
#include "reader.h"    // HIPO library header

using namespace std;
namespace fs = std::filesystem;

int main() {
    // Path to the text file containing directory paths
    string dirListFile = "directories.txt";
    ifstream inputFile(dirListFile);

    if (!inputFile.is_open()) {
        cerr << "Error: Could not open " << dirListFile << endl;
        return 1;
    }

    // Read directories from the file
    vector<string> directories;
    string dir;
    while (getline(inputFile, dir)) {
        if (!dir.empty()) {
            directories.push_back(dir);
        }
    }
    inputFile.close();

    if (directories.empty()) {
        cerr << "Error: No directories found in " << dirListFile << endl;
        return 1;
    }

    const int EVENT_LIMIT = 10000000;  // 10 million event limit
    int totalEventCount = 0;
    bool limitReached = false;

    // Iterate over directories
    for (const auto& dir : directories) {
        cout << "Processing directory: " << dir << endl;

        // Find all .hipo files in the directory
        vector<string> hipoFiles;
        for (const auto& entry : fs::directory_iterator(dir)) {
            if (entry.is_regular_file() && entry.path().extension() == ".hipo") {
                hipoFiles.push_back(entry.path().string());
            }
        }

        if (hipoFiles.empty()) {
            cout << "  No .hipo files found in directory: " << dir << endl;
            continue;
        }

        // Iterate over files in the directory
        for (const auto& file : hipoFiles) {
            cout << "  Opening file: " << file << endl;

            // Open HIPO file
            hipo::reader reader;
            reader.open(file.c_str());

            // Process events in the file
            hipo::event event;
            while (reader.next() && !limitReached) {
                totalEventCount++;
                if (totalEventCount >= EVENT_LIMIT) {
                    limitReached = true;
                    cout << "Reached 10 million events at:" << endl;
                    cout << "  Directory: " << dir << endl;
                    cout << "  File: " << file << endl;
                    cout << "  Total events processed: " << totalEventCount << endl;
                    break;
                }
            }

            
            if (limitReached) break;
        }

        if (limitReached) break;
    }

    if (!limitReached) {
        cout << "Processed all files, but less than 10 million events in total." << endl;
        cout << "Total events processed: " << totalEventCount << endl;
    }

    return 0;
}
