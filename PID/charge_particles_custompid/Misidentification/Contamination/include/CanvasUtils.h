#ifndef CANVAS_UTILS_H
#define CANVAS_UTILS_H

#include <TCanvas.h>
#include <TH1F.h>
#include <vector>
#include <iostream>
#include <map>

TCanvas* createCanvas(const std::vector<TH1F*>& histograms = {}, int rows = 1, int cols = 1,
                      const std::vector<int>& subplotIndices = {1}, const char* canvasName = "canvas",
                      int canvasWidth = 1600, int canvasHeight = 800) {
    // Validate inputs
    if (rows < 1 || cols < 1) {
        std::cerr << "Error: Rows and columns must be positive." << std::endl;
        return nullptr;
    }
    if (subplotIndices.empty() && !histograms.empty()) {
        std::cerr << "Error: Subplot indices must be provided if histograms are given." << std::endl;
        return nullptr;
    }
    for (int idx : subplotIndices) {
        if (idx < 1 || idx > rows * cols) {
            std::cerr << "Error: Invalid subplot index " << idx << "." << std::endl;
            return nullptr;
        }
    }

    // Create canvas
    TCanvas* canvas = new TCanvas(canvasName, canvasName, canvasWidth, canvasHeight);
    canvas->Divide(cols, rows);

    // If no histograms, return empty divided canvas
    if (histograms.empty()) {
        canvas->Update();
        return canvas;
    }

    // Track which pads have already been drawn to handle overlaying with "SAME"
    std::map<int, bool> padDrawn;

    // Plot histograms
    for (size_t i = 0; i < histograms.size() && i < subplotIndices.size(); ++i) {
        if (!histograms[i]) {
            std::cerr << "Warning: Histogram " << i << " is null, skipping." << std::endl;
            continue;
        }
        canvas->cd(subplotIndices[i]);

        // Determine draw option based on histogram style
        TString drawOption = (histograms[i]->GetFillStyle() == 0) ? "P" : "HIST";

        // If this pad has already been drawn, use "SAME" to overlay
        if (padDrawn[subplotIndices[i]]) {
            drawOption += " SAME";
        }

        histograms[i]->Draw(drawOption);
        padDrawn[subplotIndices[i]] = true;

        // Draw any associated fit function
        if (TF1* fit = histograms[i]->GetFunction("fit")) {
            fit->Draw("SAME");
        }
    }

    canvas->Update();
    return canvas;
}

#endif // CANVAS_UTILS_H