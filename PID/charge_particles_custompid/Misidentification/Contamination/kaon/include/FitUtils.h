#ifndef FIT_UTILS_H
#define FIT_UTILS_H

#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <map>
#include <iostream>
#include <string>
#include <iomanip>

// Custom double Crystal Ball function
Double_t doubleCrystalBall(Double_t* x, Double_t* par) {
    Double_t t = (x[0] - par[1]) / par[2]; // t = (x - mean) / sigma
    Double_t absAlphaLeft = TMath::Abs(par[3]);
    Double_t absAlphaRight = TMath::Abs(par[6]);

    // Left tail (power-law for t < -|alpha_left|)
    Double_t left = (t < -absAlphaLeft) ? 
        par[0] * TMath::Power(par[4] / absAlphaLeft, par[4]) * TMath::Exp(-0.5 * absAlphaLeft * absAlphaLeft) /
        TMath::Power(par[4] / absAlphaLeft - absAlphaLeft - t, par[4]) : 
        par[0] * TMath::Exp(-0.5 * t * t);

    // Right tail (power-law for t > |alpha_right|)
    Double_t right = (t > absAlphaRight) ? 
        par[5] * TMath::Power(par[7] / absAlphaRight, par[7]) * TMath::Exp(-0.5 * absAlphaRight * absAlphaRight) /
        TMath::Power(par[7] / absAlphaRight - absAlphaRight + t, par[7]) : 
        par[5] * TMath::Exp(-0.5 * t * t);

    return left + right; // Sum of both CB contributions
}

std::map<std::string, double> fitHistogram(TH1F* hist, const char* fitName, double fitMin, double fitMax, const char* fitType = "gaus") {
    std::map<std::string, double> fitParams;

    if (!hist || hist->GetEntries() < 50) {
        std::cerr << "Error: Invalid histogram or too few entries (" << (hist ? hist->GetEntries() : 0) << ")." << std::endl;
        return fitParams;
    }

    // Dynamically determine fit range for both tails
    double mean = hist->GetMean();
    double sigma = hist->GetRMS();
    double maxVal = hist->GetMaximum();
    int maxBin = hist->GetMaximumBin();
    double maxX = hist->GetBinCenter(maxBin);
    double leftTailLimit = maxX - 4.0 * sigma; // Extended left tail
    double rightTailLimit = maxX + 4.0 * sigma; // Extended right tail
    double min = (fitMin < leftTailLimit) ? fitMin : leftTailLimit;
    double max = (fitMax > rightTailLimit) ? fitMax : rightTailLimit;

    // Gaussian fit with standardized name "gausFit"
    TF1* gausFit = new TF1("gausFit", "gaus", min, max);
    gausFit->SetParameters(maxVal, mean, sigma);
    gausFit->SetParLimits(2, 1e-6, sigma * 3.0); // Enforce minimum sigma
    hist->Fit(gausFit, "QR+"); // Quiet, add to histogram
    fitParams["gaus_constant"] = gausFit->GetParameter(0); // Add the Constant parameter
    fitParams["gaus_mean"] = gausFit->GetParameter(1);
    fitParams["gaus_sigma"] = std::max(gausFit->GetParameter(2), 1e-6); // Prevent zero sigma
    fitParams["gaus_chi2"] = gausFit->GetNDF() > 0 ? gausFit->GetChisquare() / gausFit->GetNDF() : 0.0;
    std::cout << "Gaussian fit for " << fitName << ": mean=" << fitParams["gaus_mean"]
              << ", sigma=" << std::fixed << std::setprecision(6) << fitParams["gaus_sigma"] 
              << ", chi2/NDF=" << fitParams["gaus_chi2"] << std::endl;

    // Apply Crystal Ball fit only if fitType is "both"
    if (std::string(fitType) == "both") {
        // Estimate tail asymmetry to determine initial alpha
        double leftTailHeight = hist->GetBinContent(hist->FindBin(maxX - 2.0 * sigma));
        double rightTailHeight = hist->GetBinContent(hist->FindBin(maxX + 2.0 * sigma));
        double tailAsymmetry = (rightTailHeight - leftTailHeight) / (rightTailHeight + leftTailHeight + 1e-6);
        double initialAlphaLeft = (tailAsymmetry < -0.1) ? 1.5 : 0.5; // Stronger left tail
        double initialAlphaRight = (tailAsymmetry > 0.1) ? 1.5 : 0.5; // Stronger right tail

        // Double Crystal Ball fit using custom function
        TF1* cbFit = new TF1("cbFit", doubleCrystalBall, min, max, 8);
        cbFit->SetParameters(
            gausFit->GetParameter(0) * 0.5, // Constant for left CB
            gausFit->GetParameter(1),       // Mean
            gausFit->GetParameter(2),       // Sigma
            initialAlphaLeft,               // Alpha left
            7.0,                            // n left
            gausFit->GetParameter(0) * 0.5, // Constant for right CB
            initialAlphaRight,              // Alpha right
            7.0                             // n right
        );
        cbFit->SetParLimits(0, 0.0, maxVal);        // Constant left
        cbFit->SetParLimits(1, mean - 2.0 * sigma, mean + 2.0 * sigma);
        cbFit->SetParLimits(2, 0.0, sigma * 3.0);
        cbFit->SetParLimits(3, 0.5, 5.0);           // Alpha left
        cbFit->SetParLimits(4, 1.0, 30.0);          // n left
        cbFit->SetParLimits(5, 0.0, maxVal);        // Constant right
        cbFit->SetParLimits(6, 0.5, 5.0);           // Alpha right
        cbFit->SetParLimits(7, 1.0, 30.0);          // n right
        hist->Fit(cbFit, "QR+R");                   // Quiet, robust fitting
        fitParams["cb_mean"] = cbFit->GetParameter(1);
        fitParams["cb_sigma"] = cbFit->GetParameter(2);
        fitParams["cb_alpha_left"] = cbFit->GetParameter(3);
        fitParams["cb_n_left"] = cbFit->GetParameter(4);
        fitParams["cb_alpha_right"] = cbFit->GetParameter(6);
        fitParams["cb_n_right"] = cbFit->GetParameter(7);
        fitParams["cb_chi2"] = cbFit->GetNDF() > 0 ? cbFit->GetChisquare() / cbFit->GetNDF() : 0.0;
        std::cout << "Double Crystal Ball fit for " << fitName << ": mean=" << fitParams["cb_mean"]
                  << ", sigma=" << std::fixed << std::setprecision(6) << fitParams["cb_sigma"]
                  << ", alpha_left=" << fitParams["cb_alpha_left"]
                  << ", n_left=" << fitParams["cb_n_left"]
                  << ", alpha_right=" << fitParams["cb_alpha_right"]
                  << ", n_right=" << fitParams["cb_n_right"]
                  << ", chi2/NDF=" << fitParams["cb_chi2"] << std::endl;
    }

    return fitParams;
}

// Exponential function implementation
Double_t exponential(Double_t* x, Double_t* par) {
    return par[0] * TMath::Exp(-par[1] * x[0]);
}

// Function to fit a histogram with exponential function for a specific bin range
void fitExponentialInBinRange(TH1F* hist, double fitMin, double fitMax) {
    // Check if the histogram is valid
    if (!hist || hist->GetEntries() < 50) {
        std::cerr << "Error: Invalid histogram or too few entries (" << (hist ? hist->GetEntries() : 0) << ")." << std::endl;
        return;
    }

    // Create an exponential fit function with initial parameters
    TF1* expFit = new TF1("expFit", exponential, fitMin, fitMax, 2); // 2 parameters for exponential: amplitude and decay constant
    expFit->SetParameters(1.0, 0.1); // Initial guess for the parameters (amplitude = 1, lambda = 0.1)

    // Fit the histogram in the specified bin range
    hist->Fit(expFit, "QR+", "", fitMin, fitMax); // QR+ suppresses output and adds the fit to the histogram

    // Output fit parameters
    std::cout << "Exponential fit for range [" << fitMin << ", " << fitMax << "]:" << std::endl;
    std::cout << "Amplitude: " << expFit->GetParameter(0) << ", Lambda: " << expFit->GetParameter(1) << std::endl;
    std::cout << "Chi2/NDF: " << (expFit->GetNDF() > 0 ? expFit->GetChisquare() / expFit->GetNDF() : 0.0) << std::endl;
}

#endif // FIT_UTILS_H