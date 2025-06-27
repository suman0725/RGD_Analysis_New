#ifndef CONTAMINATION_UTILS_H
#define CONTAMINATION_UTILS_H
#include <TMath.h>
#include <map>
#include <string>
#include <iostream>
#include <algorithm>
#include <ostream>
#include <TF1.h>
/**
 * Calculate contamination using user-defined pion range and kaon right boundary
 * @param fitParams1 Fit parameters for the first Gaussian (e.g., pions)
 * @param fitParams2 Fit parameters for the second Gaussian (e.g., kaons)
 * @param pionLeft Left boundary for pion integration
 * @param pionRight Right boundary for pion integration
 * @param kaonRight Right boundary for kaon integration
 * @param c1 Reference to store the first intersection point
 * @param c2 Reference to store the second intersection point
 * @return Contamination percentage
 */
double calculateOverlapContamination(const std::map<std::string, double>& fitParams1, 
    const std::map<std::string, double>& fitParams2, 
    double pionLeft, double pionRight, 
    double kaonRight, 
    double& c1, double& c2) {
    // Check if parameters for the first Gaussian (pions) are valid
    if (fitParams1.find("gaus_constant") == fitParams1.end() ||
    fitParams1.find("gaus_mean") == fitParams1.end() ||
    fitParams1.find("gaus_sigma") == fitParams1.end() ||
    fitParams1.at("gaus_constant") <= 0 ||
    fitParams1.at("gaus_sigma") <= 0) {
    std::cout << "Invalid pion parameters\n";
    return -1.0; // Invalid parameters
    }
    double C1 = fitParams1.at("gaus_constant");
    double mu1 = fitParams1.at("gaus_mean");
    double sigma1 = fitParams1.at("gaus_sigma");

    // Check if parameters for the second Gaussian (kaons) are valid
    bool kaonsValid = true;
    if (fitParams2.find("gaus_constant") == fitParams2.end() ||
    fitParams2.find("gaus_mean") == fitParams2.end() ||
    fitParams2.find("gaus_sigma") == fitParams2.end() ||
    fitParams2.at("gaus_constant") <= 0 ||
    fitParams2.at("gaus_sigma") <= 0) {
    kaonsValid = false;
    }
    double C2 = kaonsValid ? fitParams2.at("gaus_constant") : 0.0;
    double mu2 = kaonsValid ? fitParams2.at("gaus_mean") : 0.0;
    double sigma2 = kaonsValid ? fitParams2.at("gaus_sigma") : 0.0;

    // If kaons are missing, return 0 contamination (no overlap)
    if (!kaonsValid) {
    c1 = 0.0;
    c2 = 0.0;
    std::cout << "Kaon data missing, setting contamination to 0\n";
    return 0.0;
    }

    // Ensure mu1 > mu2 (pions have higher beta)
    if (mu1 < mu2) {
    std::swap(mu1, mu2);
    std::swap(sigma1, sigma2);
    std::swap(C1, C2);
    }

    // Compute intersection points using the quadratic equation
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
    std::cout << "No valid intersection\n";
    return -1.0;
    }

    // Select intersection point between mu2 (kaons) and mu1 (pions)
    double kaonLeft = -1.0;
    if (c1 >= mu2 && c1 <= mu1) {
    kaonLeft = c1;
    } else if (c2 >= mu2 && c2 <= mu1) {
    kaonLeft = c2;
    } else {
        // If no intersection lies between mu2 and mu1, choose the point closest to kaon peak (mu2)
        kaonLeft = (std::abs(c1 - mu2) < std::abs(c2 - mu2)) ? c1 : c2;
        std::cout << "No intersection between mu2 = " << mu2 << " and mu1 = " << mu1 
                  << ", selecting point closest to kaon peak: " << kaonLeft << std::endl;
    }

    // Create TF1 objects for integration
    TF1* pionFit = new TF1("pionFit", "gaus", pionLeft, pionRight);
    pionFit->SetParameters(C1, mu1, sigma1);
    TF1* kaonFit = new TF1("kaonFit", "gaus", kaonLeft, kaonRight);
    kaonFit->SetParameters(C2, mu2, sigma2);

    // Compute total pion area in the user-defined pion range
    double total_pion_area = pionFit->Integral(pionLeft, pionRight);
    if (total_pion_area <= 0) {
    std::cout << "Invalid pion area: " << total_pion_area << std::endl;
    delete pionFit;
    delete kaonFit;
    return -1.0;
    }

    // Compute kaon area in the user-defined kaon range (overlap with pion region)
    double kaon_area = kaonFit->Integral(kaonLeft, kaonRight);
    if (kaon_area < 0) kaon_area = 0.0;

    // Compute contamination as a percentage
    double contamination = 100.0 * kaon_area / total_pion_area;
    if (contamination < 0) contamination = 0.0;

    // Clean up
    delete pionFit;
    delete kaonFit;

    std::cout << "Pion range: [" << pionLeft << ", " << pionRight << "], Kaon range: [" 
    << kaonLeft << ", " << kaonRight << "], Contamination: " << contamination << "%\n";

    return contamination;
}

#endif // CONTAMINATION_UTILS_H