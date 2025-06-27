#include <cstdlib>
#include <iostream>
#include "reader.h"
#include <map>
#include <cmath>

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


   //while (reader.hasNext()==true) { 
    if (reader.gotoEvent(0)){// Go to the first event
    reader.read(event);      // Read the event data
    event.getStructure(PART);
    event.getStructure(EVENT);
    event.getStructure(SCIN);


    // Load the map between SCIN and Particle banks using the `pindex` variable
    IndexMap scinMap = loadMapByIndex(SCIN, "pindex");
    
     PART.show();  
   // SCIN.show();                // Show REC::Particle bank data
    //EVENT.show();            // Show REC::Event bank data

    int nerows = EVENT.getRows();
    //std::cout << "REC::Event contains " << nerows << " event rows" << std::endl;

    float startTime = EVENT.getFloat("startTime", 0);
    //if (startTime <= 0) continue;
    
    int nrows = PART.getRows();

    //cout << "Total rows in REC:: Particle banks are: " << nrows << endl;
    for (int i = 0; i < nrows; i++) {
    // Get particle properties
    int pid = PART.getInt("pid", i);
    float px = PART.getFloat("px", i);
    float py = PART.getFloat("py", i);
    float pz = PART.getFloat("pz", i);
    int charge = PART.getInt("charge", i);
    float p = sqrt(px * px + py * py + pz * pz); 
    float beta = 0; 

    
     beta = PART.getFloat("beta",i);

    float beta_p = p / sqrt(p * p + 0.13957 * 0.13957);

    cout << beta << " " << beta_p << endl; 
    


 
    /* if (charge <= 0) continue;

   
    float beta = PART.getFloat("beta", i);
    cout << "Beta measured from TOF = " << beta 
         << " with pid value = " << pid 
         << " for row index " << i << endl;

   
    int scirows = SCIN.getRows();
    cout << "Total rows for SCIN are: " << scirows << endl;

    if (scinMap.find(i) != scinMap.end()) {
                for (int iScinRow : scinMap[i]) {

        // Get path and time for beta calculation
        float path = SCIN.getFloat("path", iScinRow);
        float time = SCIN.getFloat("time", iScinRow);

      
            float beta_from_scin = path/time ;
            cout << "Beta from SCIN = " << beta_from_scin 
                 << " with row index = " << i << endl;

                }
    } */
}

}
}