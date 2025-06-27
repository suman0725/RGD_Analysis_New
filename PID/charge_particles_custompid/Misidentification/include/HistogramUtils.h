#ifndef HISTOGRAM_UTILS_H
#define HISTOGRAM_UTILS_H

#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <string>
#include <vector>
#include <map>
#include "BinUtils.h"

struct HistoConfig {
    enum class Type { OneD, TwoD };
    Type type;
    std::string name, var;
    double xMin, xMax;
    int xBins;
    std::string yVar;
    double yMin, yMax;
    int yBins;

    HistoConfig(Type t, const std::string& n, const std::string& v, double xm, double xM, int xb,
                const std::string& yv = "", double ym = 0, double yM = 0, int yb = 0)
        : type(t), name(n), var(v), xMin(xm), xMax(xM), xBins(xb),
          yVar(yv), yMin(ym), yMax(yM), yBins(yb) {}
};

std::map<std::string, TH1*> makeHisto(
    std::vector<TTree*> trees, // Changed from const std::vector<TTree*>& to std::vector<TTree*>
    const std::vector<std::string>& particleNames,
    const HistoConfig& config,
    const Bin* bin = nullptr
);

#endif