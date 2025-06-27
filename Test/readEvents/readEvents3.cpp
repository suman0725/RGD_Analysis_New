#include <cstdlib>
#include <iostream>
#include "reader.h"
#include <map>
#include <fstream>


using namespace std; 


// Typedef for the index map
typedef std::map<int, std::vector<int>> IndexMap;

// Function to create a map from index values
IndexMap loadMapByIndex(hipo::bank& fromBank, const char* idxVarName) {
    IndexMap map;
    if (fromBank.getRows() > 0) {
        for (int iFrom = 0; iFrom < fromBank.getRows(); ++iFrom) {
            int iTo = fromBank.getInt(idxVarName, iFrom);
            if (map.find(iTo) == map.end()) {
                map[iTo] = vector<int>();
            }
            map[iTo].push_back(iFrom);
        }
    }
    return map;
}



int main(int argc, char** argv) {

    std::cout << "Reading file example program (HIPO) " << __cplusplus << std::endl;

    char inputFile[256];

    if (argc > 1) {
        sprintf(inputFile, "%s", argv[1]);
    } else {
        std::cout << " *** Please provide a file name..." << std::endl;
        exit(1);
    }


    ofstream out("output.txt");
    if (!out.is_open()) {
        std::cerr << "Error opening file for writing." << std::endl;
        return -1; // Exit if the file cannot be opened
    }

    std::cout << "Opening file: " << inputFile << std::endl;
    hipo::reader reader;
    reader.open(inputFile);
    std::cout << "File opened successfully." << std::endl;

    hipo::dictionary factory;
    reader.readDictionary(factory);
    std::cout << "Dictionary loaded, showing schemas:" << std::endl;
    factory.show();

    hipo::event event;
    int counter = 0;
   

    hipo::bank PART(factory.getSchema("REC::Particle"));
    hipo::bank EVENT(factory. getSchema("REC::Event"));
    hipo::bank SCIN (factory.getSchema("REC::Scintillator"));


   //if (reader.gotoEvent(0)) { // Go to the first event
   while (reader.next() == true){
    if (counter >= 20) break; 
        reader.read(event);      // Read the event data
        event.getStructure(PART);
        event.getStructure(EVENT);
        event.getStructure(SCIN);
         IndexMap scinMap = loadMapByIndex(SCIN, "pindex");
        PART.show();  
        SCIN.show();               // Show REC::Particle bank data
        //EVENT.show();            // Show REC::Event bank data

        int nerows = EVENT.getRows();
        //std::cout << "REC::Event contains " << nerows << " event rows" << std::endl;

        float startTime = EVENT.getFloat("startTime", 0);
        int nrows = PART.getRows();
        // cout << "Total rows in REC:: Particle banks are: " << nrows << endl;
        for (int i = 0; i < nrows; i++) {
            // Get particle properties
            int pid = PART.getInt("pid", i);
            float px = PART.getFloat("px", i);
            float py = PART.getFloat("py", i);
            float pz = PART.getFloat("pz", i);
            int charge = PART.getInt("charge", i);

            // Skip neutral or negatively charged particles
            // if (charge <= 0) continue;

            // Get beta from REC::Particle for positively charged particles
            float beta = PART.getFloat("beta", i);
            

            // Process associated SCIN rows for this particle
            int scirows = SCIN.getRows();
            //cout << "Total rows for SCIN are: " << scirows << endl;



            if (scinMap.find(i) != scinMap.end()) {
                for (int j : scinMap[i]) {
                    int detector = SCIN.getByte("detector",j);
                    //cout << "detector = " << detector << endl;
                    if ( detector == 12) { 

                        // Get path and time for beta calculation
                        float path = SCIN.getFloat("path", j);
                        float time = SCIN.getFloat("time", j);

                        out << "Beta measured from TOF = " << beta 
                            << " with pid value = " << pid 
                            << " for row index REC::Particle " << i << std::endl;

                        float beta_from_scin = path / ((time - startTime) * 30);
                        out << "Beta from SCIN = " << beta_from_scin 
                            << " with row index REC::Particle = " << i 
                            << " time = " << time 
                            << " path = " << path 
                            << " startTime = " << startTime << std::endl;

                    }
                }
            }

        }
        counter++;
    }

    out.close();

}