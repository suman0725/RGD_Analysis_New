#include "FitFunctions.h"
#include <TF1.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TLine.h>
#include <iostream>

FitResult fitCrystalBall(TH1F* hist, const std::string& name, double fitMin, double fitMax, TCanvas* canvas, int color, const std::string& particle) {
    FitResult result;
    if (!hist || hist->GetEntries() == 0) {
        std::cerr << "Error: Invalid histogram for fitting: " << name << std::endl;
        result.status = -1;
        return result;
    }

    // Initial Gaussian fit to estimate parameters
    TF1* gauss = new TF1("gaussFit", "gaus", fitMin, fitMax);
    hist->Fit(gauss, "RQ"); // R: fit range, Q: quiet
    result.gaussConst = gauss->GetParameter(0);
    result.gaussMean = gauss->GetParameter(1);
    result.gaussSigma = gauss->GetParameter(2);
    double norm = result.gaussConst;
    double alpha = 1.0; // Tail parameter
    double n = 2.0;     // Power-law exponent

    // Create Crystal Ball function with Gaussian-derived initial parameters
    TF1* cb = new TF1(name.c_str(), "crystalball", fitMin, fitMax);
    cb->SetParameters(norm, result.gaussMean, result.gaussSigma, alpha, n);
    cb->SetParLimits(0, 0, norm * 2.0); // Wider norm limit
    cb->SetParLimits(1, result.gaussMean - 5 * result.gaussSigma, result.gaussMean + 5 * result.gaussSigma); // Wider mean range
    cb->SetParLimits(2, result.gaussSigma * 0.3, result.gaussSigma * 3.0); // Wider sigma range
    cb->SetParLimits(3, 0.1, 10.0); // Wider alpha range
    cb->SetParLimits(4, 0.5, 20.0); // Wider n range

    // Fit the histogram and get fit status
    TFitResultPtr fitResult = hist->Fit(cb, "RQS"); // R: fit range, Q: quiet, S: store result
    result.status = fitResult->Status(); // Fit status (0 = success)

    // Extract Crystal Ball fit parameters
    result.cbConst = cb->GetParameter(0);
    result.mean = cb->GetParameter(1);
    result.sigma = cb->GetParameter(2);
    result.alpha = cb->GetParameter(3);
    result.n = cb->GetParameter(4);

    // Draw the fit if canvas is provided
    if (canvas) {
        cb->SetLineColor(color);
        cb->SetLineStyle(1);
        cb->Draw("SAME");
        canvas->Update();

        // Add vertical lines at mean and mean ± 3 * sigma
        TLine* lineMean = new TLine(result.mean, 0, result.mean, hist->GetMaximum());
        TLine* lineLower = new TLine(result.mean - 3 * result.sigma, 0, result.mean - 3 * result.sigma, hist->GetMaximum());
        TLine* lineUpper = new TLine(result.mean + 3 * result.sigma, 0, result.mean + 3 * result.sigma, hist->GetMaximum());
        lineMean->SetLineColor(kRed);
        lineLower->SetLineColor(kRed);
        lineUpper->SetLineColor(kRed);
        lineMean->SetLineStyle(1);  // Solid line for mean
        lineLower->SetLineStyle(2); // Dashed line for ±3 sigma
        lineUpper->SetLineStyle(2);
        lineMean->Draw("SAME");
        lineLower->Draw("SAME");
        lineUpper->Draw("SAME");
        canvas->Update();
    }

    delete gauss; // Clean up Gaussian fit function
    return result;
}