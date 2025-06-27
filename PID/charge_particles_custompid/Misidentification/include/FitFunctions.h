#ifndef FIT_FUNCTIONS_H
#define FIT_FUNCTIONS_H

#include <TH1.h>
#include <TCanvas.h>
#include <string>

struct FitResult {
    double gaussConst;   // Gaussian constant
    double gaussMean;    // Gaussian mean
    double gaussSigma;   // Gaussian sigma
    double cbConst;      // Crystal Ball constant
    double mean;         // Crystal Ball mean
    double sigma;        // Crystal Ball sigma
    double alpha;        // Crystal Ball alpha
    double n;            // Crystal Ball n
    int status;          // Fit status
};

FitResult fitCrystalBall(TH1F* hist, const std::string& name, double fitMin, double fitMax, TCanvas* canvas, int color, const std::string& particle);

#endif