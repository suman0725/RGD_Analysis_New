#ifndef CONTAMINATION_UTILS_H
#define CONTAMINATION_UTILS_H

#include <TMath.h>
#include <map>
#include <string>
#include <cmath>
#include <TF1.h>
#include <iostream>

/**
 * Calculate contamination by finding the kaon area in the pion region
 * @param fitParams1 Fit parameters for the first Gaussian (e.g., pions)
 * @param fitParams2 Fit parameters for the second Gaussian (e.g., kaons)
 * @param c1 Reference to store the first intersection point
 * @param c2 Reference to store the second intersection point
 * @return Contamination percentage
 */
double calculateOverlapContamination(const std::map<std::string, double>& fitParams1, const std::map<std::string, double>& fitParams2, 
                                     double& c1, double& c2) {
    // Extract parameters for the first Gaussian (pions)
    if (fitParams1.find("gaus_constant") == fitParams1.end() ||
        fitParams1.find("gaus_mean") == fitParams1.end() ||
        fitParams1.find("gaus_sigma") == fitParams1.end()) {
        return -1.0; // Invalid parameters
    }
    double C1 = fitParams1.at("gaus_constant");
    double mu1 = fitParams1.at("gaus_mean");
    double sigma1 = fitParams1.at("gaus_sigma");

    // Extract parameters for the second Gaussian (kaons)
    if (fitParams2.find("gaus_constant") == fitParams2.end() ||
        fitParams2.find("gaus_mean") == fitParams2.end() ||
        fitParams2.find("gaus_sigma") == fitParams2.end()) {
        return -1.0; // Invalid parameters
    }
    double C2 = fitParams2.at("gaus_constant");
    double mu2 = fitParams2.at("gaus_mean");
    double sigma2 = fitParams2.at("gaus_sigma");

    // Prevent division by zero or negative values
    if (sigma1 <= 0 || sigma2 <= 0 || C1 <= 0 || C2 <= 0) {
        return -1.0;
    }

    // Ensure mu1 > mu2 (pions have higher beta)
    if (mu1 < mu2) {
        std::swap(mu1, mu2);
        std::swap(sigma1, sigma2);
        std::swap(C1, C2);
    }

    // Create TF1 objects for integration
    TF1* pionFit = new TF1("pionFit", "gaus", 0.96, 1.0);
    pionFit->SetParameters(C1, mu1, sigma1);
    TF1* kaonFit = new TF1("kaonFit", "gaus", 0.96, 1.0);
    kaonFit->SetParameters(C2, mu2, sigma2);

    // Calculate pion endpoints (where Gaussian drops to 10^-5 of peak height or < 0.1, capped between 0.96 and 1.0)
    double pionLeftEndpoint = mu1;
    double threshold = std::min(0.1, C1 * 1e-5); // Use 10^-5 of peak or 0.1, whichever is smaller
    while (pionFit->Eval(pionLeftEndpoint) > threshold && pionLeftEndpoint > 0.96) {
        pionLeftEndpoint -= 0.001 * sigma1; // Finer step size
    }
    pionLeftEndpoint = std::max(0.96, pionLeftEndpoint);

    double pionRightEndpoint = mu1;
    while (pionFit->Eval(pionRightEndpoint) > threshold && pionRightEndpoint < 1.03) {
        pionRightEndpoint += 0.001 * sigma1;
    }
    pionRightEndpoint = std::min(1.03, pionRightEndpoint);

    // Calculate kaon right endpoint
    double kaonRightEndpoint = mu2;
    threshold = std::min(0.1, C2 * 1e-5); // Use 10^-5 of peak or 0.1, whichever is smaller
    while (kaonFit->Eval(kaonRightEndpoint) > threshold && kaonRightEndpoint < 1.0) {
        kaonRightEndpoint += 0.001 * sigma2;
    }
    kaonRightEndpoint = std::min(1.0, kaonRightEndpoint);

    // Compute the total pion area (from left to right endpoint)
    //double total_pion_area = pionFit->Integral(pionLeftEndpoint, pionRightEndpoint);
    double total_pion_area = pionFit->Integral(0.9880, 1.0120);
    if (total_pion_area <= 0) {
        delete pionFit;
        delete kaonFit;
        return -1.0;
    }

    // Compute the intersection points using the quadratic equation
    double sigma1_sq = sigma1 * sigma1;
    double sigma2_sq = sigma2 * sigma2;
    double a = 1.0 / (2 * sigma2_sq) - 1.0 / (2 * sigma1_sq);
    double b = mu1 / sigma1_sq - mu2 / sigma2_sq;
    double c_coeff = (mu2 * mu2) / (2 * sigma2_sq) - (mu1 * mu1) / (2 * sigma1_sq) - TMath::Log(C2 / C1);
    double discriminant = b * b - 4 * a * c_coeff;

    c1 = 0.0;
    c2 = 0.0;
    if (discriminant >= 0) {
        double sqrt_disc = TMath::Sqrt(discriminant);
        c1 = (-b + sqrt_disc) / (2 * a);
        c2 = (-b - sqrt_disc) / (2 * a);
        std::cout << "Intersection points: c1 = " << c1 << ", c2 = " << c2 << std::endl;
    } else {
        std::cout << "No valid intersection" << std::endl;
        delete pionFit;
        delete kaonFit;
        return -1.0;
    }

   
    // Select the intersection point between mu2 and mu1
    double c = -1.0;
    if (c1 >= mu2 && c1 <= mu1) {
        c = c1;
        std::cout << "Using c1 = " << c1 << " as the intersection point (between mu2 = " << mu2 << " and mu1 = " << mu1 << ")" << std::endl;
    } else if (c2 >= mu2 && c2 <= mu1) {
        c = c2;
        std::cout << "Using c2 = " << c2 << " as the intersection point (between mu2 = " << mu2 << " and mu1 = " << mu1 << ")" << std::endl;
    } else {
        std::cout << "No intersection point between mu2 = " << mu2 << " and mu1 = " << mu1 << std::endl;
        delete pionFit;
        delete kaonFit;
        return -1.0;
    }

    // Ensure c is within the valid range
    c = std::max(0.96, std::min(c, 0.998));

    // Compute the kaon area in the pion region (from c to kaonRightEndpoint)
    //double kaon_area = kaonFit->Integral(c, kaonRightEndpoint);
    double kaon_area = kaonFit->Integral(c,1);

    // Compute contamination as a percentage
    double contamination = 100.0 * kaon_area / total_pion_area;

    // Clean up
    delete pionFit;
    delete kaonFit;

    return contamination >= 0 ? contamination : 0.0;
}

#endif // CONTAMINATION_UTILS_H