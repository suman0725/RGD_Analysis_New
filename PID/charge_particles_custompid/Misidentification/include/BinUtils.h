#ifndef BIN_UTILS_H
#define BIN_UTILS_H

#include <string>
#include <vector> // Added for std::vector

struct Bin {
    double min, max;
    std::string varName, unit;

    Bin(double m, double M, const std::string& v = "", const std::string& u = "")
        : min(m), max(M), varName(v), unit(u) {}
};

std::vector<Bin> createBins(double min, double max, double step, const std::string& varName, const std::string& unit);

#endif