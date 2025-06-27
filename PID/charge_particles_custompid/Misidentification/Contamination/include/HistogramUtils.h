#ifndef HISTOGRAM_UTILS_H
#define HISTOGRAM_UTILS_H

#include <TH1F.h>
#include <TTree.h>
#include <TString.h>
#include <TBranch.h>
#include <vector>
#include <iostream>
#include <limits>

TH1F* createHistogram(TTree* tree, const char* branchName, const std::vector<double>& bins = {},
                      const char* title = "Histogram", const char* xLabel = "Value",
                      const char* yLabel = "Frequency", int lineColor = kBlue,
                      int fillColor = kBlue, int fillStyle = 3001, const char* histName = "hist",
                      double minVal = 0, double maxVal = 0, const char* drawStyle = "HIST",
                      double markerSize = 1.0) {
    // Validate inputs
    if (!tree || !branchName || tree->GetEntries() == 0) {
        std::cerr << "Error: Invalid tree or branch name, or empty tree." << std::endl;
        return nullptr;
    }
    if (minVal != 0 && maxVal != 0 && minVal >= maxVal) {
        std::cerr << "Error: Invalid range (minVal >= maxVal)." << std::endl;
        return nullptr;
    }

    // Create histogram
    TH1F* hist;
    if (!bins.empty()) {
        // Use provided bins
        hist = new TH1F(histName, title, bins.size() - 1, bins.data());
    } else if (minVal != 0 && maxVal != 0) {
        // Use specified range with fixed number of bins
        hist = new TH1F(histName, title, 100, minVal, maxVal);
    } else {
        // Auto-determine range from tree data
        TBranch* branch = tree->GetBranch(branchName);
        if (!branch) {
            std::cerr << "Error: Branch " << branchName << " not found." << std::endl;
            return nullptr;
        }

        // Check branch type
        TClass* expectedClass = nullptr;
        EDataType expectedType;
        branch->GetExpectedType(expectedClass, expectedType);
        double autoMin = std::numeric_limits<double>::max();
        double autoMax = std::numeric_limits<double>::lowest();

        if (expectedType == kFloat_t) {
            Float_t value;
            tree->SetBranchAddress(branchName, &value);
            for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
                tree->GetEntry(i);
                if (value < autoMin) autoMin = value;
                if (value > autoMax) autoMax = value;
            }
        } else if (expectedType == kDouble_t) {
            Double_t value;
            tree->SetBranchAddress(branchName, &value);
            for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
                tree->GetEntry(i);
                if (value < autoMin) autoMin = value;
                if (value > autoMax) autoMax = value;
            }
        } else {
            std::cerr << "Error: Unsupported branch type for " << branchName << std::endl;
            return nullptr;
        }

        // Add small padding to ensure all data is included
        double padding = (autoMax - autoMin) * 0.05;
        if (padding == 0) padding = 0.1; // Avoid zero padding for flat distributions
        hist = new TH1F(histName, title, 30, autoMin - padding, autoMax + padding);
        tree->ResetBranchAddresses();
    }

    // Fill histogram
    TString drawCmd = TString::Format("%s>>%s", branchName, histName);
    tree->Draw(drawCmd, "", "goff");

    // Set histogram style
    hist->SetXTitle(xLabel);
    hist->SetYTitle(yLabel);
    hist->SetLineColor(lineColor);
    if (TString(drawStyle) == "HIST") {
        hist->SetFillColor(fillColor);
        hist->SetFillStyle(fillStyle);
    } else if (TString(drawStyle) == "P") {
        hist->SetMarkerStyle(kFullCircle); // Circle markers for points
        hist->SetMarkerSize(markerSize);   // Set marker size
        hist->SetMarkerColor(lineColor);   // Match marker color to line color
        hist->SetFillStyle(0);             // No fill for points
    } else {
        std::cerr << "Warning: Unknown drawStyle '" << drawStyle << "'. Defaulting to HIST." << std::endl;
        hist->SetFillColor(fillColor);
        hist->SetFillStyle(fillStyle);
    }
    hist->SetStats(0);

    return hist;
}

#endif // HISTOGRAM_UTILS_H