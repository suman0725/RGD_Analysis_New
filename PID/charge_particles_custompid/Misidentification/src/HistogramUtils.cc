#include "HistogramUtils.h"
#include <iostream>

std::map<std::string, TH1*> makeHisto(
    std::vector<TTree*> trees,
    const std::vector<std::string>& particleNames,
    const HistoConfig& config,
    const Bin* bin
) {
    std::map<std::string, TH1*> histMap;

    // Validate inputs
    if (trees.empty() || particleNames.empty()) {
        std::cerr << "Error: Empty trees or particleNames in makeHisto" << std::endl;
        return histMap;
    }

    // Construct the full cut string
    std::string fullCut = "chi2pid != 99999.0"; // Hardcoded base cut
    if (bin) {
        // Add bin range cut (assuming the bin variable is 'p' for momentum)
        std::string binCut = Form("p >= %.2f && p < %.2f", bin->min, bin->max);
        fullCut += " && " + binCut;
    }
    std::cout << "Applying cut: " << (fullCut.empty() ? "none" : fullCut) << std::endl;

    // Create histograms for each particle
    for (size_t i = 0; i < particleNames.size(); ++i) {
        const std::string& particle = particleNames[i];
        if (particle.empty()) {
            std::cerr << "Error: Empty particle name at index " << i << std::endl;
            continue;
        }

        std::string histName;
        if (bin) {
            histName = Form("%s_%s_%.2f_%.2f", config.name.c_str(), particle.c_str(), bin->min, bin->max);
        } else {
            histName = Form("%s_%s", config.name.c_str(), particle.c_str());
        }

        // Create histogram based on type
        TH1* hist = nullptr;
        if (config.type == HistoConfig::Type::OneD) {
            hist = new TH1F(histName.c_str(), histName.c_str(), config.xBins, config.xMin, config.xMax);
        } else { // TwoD
            hist = new TH2F(histName.c_str(), histName.c_str(),
                            config.xBins, config.xMin, config.xMax,
                            config.yBins, config.yMin, config.yMax);
        }

        // Fill histogram from trees
        for (auto* tree : trees) {
            if (!tree) {
                std::cerr << "Warning: Null tree pointer in makeHisto" << std::endl;
                continue;
            }
            std::string drawCmd;
            if (config.type == HistoConfig::Type::OneD) {
                drawCmd = Form("%s>>+%s", config.var.c_str(), histName.c_str());
            } else { // TwoD
                drawCmd = Form("%s:%s>>+%s", config.yVar.c_str(), config.var.c_str(), histName.c_str());
            }
            tree->Draw(drawCmd.c_str(), fullCut.c_str(), "goff");
        }

        // Debug output
        std::cout << particle << " histogram " << histName << ": entries=" << hist->GetEntries() << std::endl;
        histMap[particle] = hist;
    }

    return histMap;
}