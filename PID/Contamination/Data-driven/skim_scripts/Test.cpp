#include <cstdlib>
#include <iostream>
#include <hipo4/reader.h> // Use full path for clarity

int main() {
    std::cout << "Counting electrons and positrons in HIPO file" << std::endl;

    // Hardcode the input file
    const char* inputFile = "/lustre24/expphy/cache/hallb/scratch/rg-d/production/skim_pass0v11/CuSn/skim_run_019066.hipo";
    std::cout << "Opening file: " << inputFile << std::endl;

    // Open the HIPO file
    hipo::reader reader;
    reader.open(inputFile);
    std::cout << "File opened successfully." << std::endl;

    // Load dictionary
    hipo::dictionary factory;
    reader.readDictionary(factory);
    std::cout << "Dictionary loaded." << std::endl;

    // Initialize event and bank
    hipo::event event;
    hipo::bank PART(factory.getSchema("REC::Particle"));

    // Counters
    int counter = 0;
    int total_electron_count = 0;
    int total_positron_count = 0;
    const int maxEvents = 10000; // Process up to 10,000 events

    // Event loop
    while (reader.next() && counter < maxEvents) {
        reader.read(event);
        event.getStructure(PART);

        // Debug output for first 10 events
        if (counter < 10) {
            std::cout << "Event " << counter << ": " << PART.getRows() << " particles, PIDs = ";
            for (int i = 0; i < PART.getRows(); i++) {
                int pid = PART.getInt("pid", i);
                std::cout << pid << " ";
            }
            std::cout << std::endl;
        }

        // Count electrons and positrons
        int nrows = PART.getRows();
        for (int i = 0; i < nrows; i++) {
            int pid = PART.getInt("pid", i);
            if (pid == 11) {
                total_electron_count++;
            } else if (pid == -11) {
                total_positron_count++;
            }
        }

        counter++;
    }

    // Print results
    std::cout << "Processed events = " << counter << std::endl;
    std::cout << "Total EB electrons = " << total_electron_count << std::endl;
    std::cout << "Total EB positrons = " << total_positron_count << std::endl;

    return 0;
}