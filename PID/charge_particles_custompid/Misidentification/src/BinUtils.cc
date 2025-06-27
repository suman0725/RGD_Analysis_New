#include "BinUtils.h"
#include <vector>

std::vector<Bin> createBins(double min, double max, double step, const std::string& varName, const std::string& unit) {
    std::vector<Bin> bins;
    for (double i = min; i < max; i += step) {
        double binMax = (i + step >= max) ? max : i + step;
        if (binMax <= i) continue; // Skip zero-width bins
        bins.emplace_back(i, binMax, varName, unit);
    }
    return bins;
}
