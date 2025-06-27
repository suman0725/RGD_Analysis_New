#include <iostream>
#include <vector>
#include "CCDB/Calibration.h"
#include "CCDB/CalibrationGenerator.h"

int main() {
    // Create calibration object using CalibrationGenerator
    ccdb::Calibration *calib = ccdb::CalibrationGenerator::CreateCalibration("mysql://clas12reader@clasdb.jlab.org/clas12");

    if (!calib) {
        std::cerr << "Failed to create Calibration object!" << std::endl;
        return 1;
    }

    // Retrieve calibration data
    std::vector<std::vector<double>> values;
    if (!calib->GetCalib(values, "/calibration/ftof/tres")) {
        std::cerr << "Failed to get calibration data!" << std::endl;
        return 1;
    }

    // Print retrieved calibration values
    for (const auto &row : values) {
        for (double val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}

