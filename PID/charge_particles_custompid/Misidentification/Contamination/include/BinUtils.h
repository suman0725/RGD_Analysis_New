#ifndef BIN_UTILS_H
#define BIN_UTILS_H

#include <vector>
#include <iostream>

std::vector<double> createBins(double min, double max, int nBins) {
    // Validate inputs
    if (min >= max || nBins < 1) {
        std::cerr << "Error: Invalid bin parameters (min >= max or nBins < 1)." << std::endl;
        return {};
    }

    // Generate linear bins
    std::vector<double> bins(nBins + 1);
    double step = (max - min) / nBins;
    for (int i = 0; i <= nBins; ++i) {
        bins[i] = min + i * step;
    }

    return bins;
}

#endif // BIN_UTILS_H