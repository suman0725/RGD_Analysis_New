#ifndef CANVAS_MANAGER_H
#define CANVAS_MANAGER_H

#include <TLegend.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TLatex.h>
#include <string>
#include <vector>
#include <functional>
#include "FitFunctions.h"
#include "BinUtils.h" // Added for Bin type

class CanvasManager {
private:
    // Reordered members to match constructor initialization list
    std::string name_;
    std::string title_;
    int width_;
    int height_;
    std::string baseDir_;
    std::string projectName_;
    std::string outputFile_;
    bool useBeforeAfter_;
    std::vector<int> colors_;
    bool useLegend_;
    std::string binVariable_;

public:
    CanvasManager(const std::string& name, const std::string& title, int width, int height,
                  const std::string& baseDir, const std::string& projectName,
                  const std::string& outputFile, bool useBeforeAfter,
                  const std::vector<int>& colors, bool useLegend,
                  const std::string& binVariable);

    bool validateHistogramPair(const std::pair<std::string, TH1*>& histPair);

    void drawSingleHistogram(
        TCanvas* canvas,
        int padNum,
        const std::pair<std::string, TH1*>& histPair,
        TLegend* legend,
        FitResult& fitResult,
        const std::string& binLabel
    );

    TCanvas* plotSideBySide(
        const std::vector<std::pair<std::pair<std::string, TH1*>, std::pair<std::string, TH1*>>>& histsPerBin,
        const std::vector<Bin*>& bins,
        std::function<void(TCanvas*, size_t, TH1*, TH1*, FitResult&, FitResult&)> fitFunc
    );

    TCanvas* plotOverlapping(
        const std::vector<std::pair<std::map<std::string, TH1*>, std::map<std::string, TH1*>>>& histsPerBin,
        const std::vector<Bin*>& bins,
        bool beforeAfter = false,
        std::function<void(TCanvas*, size_t, const std::map<std::string, TH1*>&, const std::map<std::string, TH1*>)> contaminationFunc = nullptr
    );
};

#endif