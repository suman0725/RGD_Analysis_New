#include "CanvasManager.h"
#include <TSystem.h>
#include <iostream>

CanvasManager::CanvasManager(
    const std::string& name, const std::string& title, int width, int height,
    const std::string& baseDir, const std::string& projectName,
    const std::string& outputFile, bool useBeforeAfter,
    const std::vector<int>& colors, bool useLegend,
    const std::string& binVariable
) : name_(name), title_(title), width_(width), height_(height),
    baseDir_(baseDir), projectName_(projectName), outputFile_(outputFile),
    useBeforeAfter_(useBeforeAfter), colors_(colors), useLegend_(useLegend),
    binVariable_(binVariable) {}

bool CanvasManager::validateHistogramPair(const std::pair<std::string, TH1*>& histPair) {
    if (!histPair.second) {
        std::cout << "Validation failed for " << histPair.first << ": Null histogram" << std::endl;
        return false;
    }
    if (histPair.second->GetEntries() == 0) {
        std::cout << "Validation failed for " << histPair.first << ": Zero entries" << std::endl;
        return false;
    }
    std::cout << "Validation passed for " << histPair.first << ": entries=" << histPair.second->GetEntries()
              << ", NbinsX=" << histPair.second->GetNbinsX() << std::endl;
    return true;
}

void CanvasManager::drawSingleHistogram(
    TCanvas* canvas, int padNum, const std::pair<std::string, TH1*>& histPair,
    TLegend* legend, FitResult& fitResult, const std::string& binLabel
) {
    canvas->cd(padNum);
    TVirtualPad* vpad = canvas->GetPad(padNum);
    TPad* pad = dynamic_cast<TPad*>(vpad);
    if (pad) {
        pad->Clear();
        pad->SetLogy(0);
    }
    if (histPair.second) {
        histPair.second->Draw("HIST");
        std::cout << "Adding legend entries for " << histPair.first << std::endl;
        if (legend) {
            legend->AddEntry(histPair.second, histPair.first.c_str(), "l");
        }
        TLatex* label = new TLatex(0.1, 0.95, binLabel.c_str());
        label->SetNDC();
        label->Draw();
    }
}

TCanvas* CanvasManager::plotSideBySide(
    const std::vector<std::pair<std::pair<std::string, TH1*>, std::pair<std::string, TH1*>>>& histsPerBin,
    const std::vector<Bin*>& bins,
    std::function<void(TCanvas*, size_t, TH1*, TH1*, FitResult&, FitResult&)> fitFunc
) {
    std::string pdfDir = baseDir_ + projectName_ + "/pdf/";
    gSystem->Exec(Form("mkdir -p %s", pdfDir.c_str()));
    TCanvas* canvas = new TCanvas(name_.c_str(), title_.c_str(), width_, height_);
    bool firstPage = true;
    int pageCount = 0;

    for (size_t binIndex = 0; binIndex < histsPerBin.size(); ++binIndex) {
        const auto& [pionHistPair, kaonHistPair] = histsPerBin[binIndex];
        std::cout << "Processing bin " << binIndex << " (histsPerBin.size() = " << histsPerBin.size() << ", bins.size() = " << bins.size() << ")" << std::endl;

        if (!validateHistogramPair(pionHistPair) || !validateHistogramPair(kaonHistPair)) {
            std::cout << "Skipping bin " << binIndex << ": invalid histograms" << std::endl;
            continue;
        }

        canvas->Clear();
        canvas->Divide(2, 1);

        std::string binLabel = Form("%.2f #leq %s < %.2f", bins[binIndex]->min, binVariable_.c_str(), bins[binIndex]->max);
        FitResult pionFit, kaonFit;
        TLegend* legendPions = useLegend_ ? new TLegend(0.6, 0.7, 0.8, 0.9) : nullptr;
        TLegend* legendKaons = useLegend_ ? new TLegend(0.6, 0.7, 0.8, 0.9) : nullptr;

        if (legendPions) std::cout << "Legend created for Pions in bin " << binIndex << std::endl;
        if (legendKaons) std::cout << "Legend created for Kaons in bin " << binIndex << std::endl;

        canvas->cd(1);
        drawSingleHistogram(canvas, 1, pionHistPair, legendPions, pionFit, binLabel);
        if (legendPions) {
            legendPions->Draw();
        }

        canvas->cd(2);
        drawSingleHistogram(canvas, 2, kaonHistPair, legendKaons, kaonFit, binLabel);
        if (legendKaons) {
            legendKaons->Draw();
        }

        if (fitFunc) {
            fitFunc(canvas, binIndex, pionHistPair.second, kaonHistPair.second, pionFit, kaonFit);
        }

        canvas->cd(0);
        canvas->Update();

        if (firstPage) {
            std::cout << "Opening PDF for first valid bin: " << binIndex << std::endl;
            canvas->Print(Form("%ssidebyside_%s(", pdfDir.c_str(), outputFile_.c_str()));
            firstPage = false;
        } else {
            std::cout << "Adding page for bin " << binIndex << std::endl;
        }

        std::string binName = Form("%.2f_%.2f", bins[binIndex]->min, bins[binIndex]->max);
        std::cout << "Bin " << binIndex << ": min=" << bins[binIndex]->min << ", max=" << bins[binIndex]->max << std::endl;
        std::cout << "Printing canvas for bin " << binIndex << " (Page " << ++pageCount << ")" << std::endl;
        canvas->Print(Form("%ssidebyside_%s", pdfDir.c_str(), outputFile_.c_str()));

        delete legendPions;
        delete legendKaons;
    }

    if (pageCount > 0) {
        std::cout << "Closing PDF (Page " << pageCount << ")" << std::endl;
        canvas->Print(Form("%ssidebyside_%s)", pdfDir.c_str(), outputFile_.c_str()));
    }

    return canvas;
}

TCanvas* CanvasManager::plotOverlapping(
    const std::vector<std::pair<std::map<std::string, TH1*>, std::map<std::string, TH1*>>>& histsPerBin,
    const std::vector<Bin*>& bins,
    bool beforeAfter,
    std::function<void(TCanvas*, size_t, const std::map<std::string, TH1*>&, const std::map<std::string, TH1*>)> contaminationFunc
) {
    std::string pdfDir = baseDir_ + projectName_ + "/pdf/";
    gSystem->Exec(Form("mkdir -p %s", pdfDir.c_str()));
    TCanvas* canvas = new TCanvas(name_.c_str(), title_.c_str(), width_, height_);
    bool firstPage = true;
    int pageCount = 0;

    for (size_t binIndex = 0; binIndex < histsPerBin.size(); ++binIndex) {
        const auto& [beforeHists, afterHists] = histsPerBin[binIndex];
        std::cout << "Processing bin " << binIndex << " (histsPerBin.size() = " << histsPerBin.size() << ", bins.size() = " << bins.size() << ")" << std::endl;

        // Validate histograms (customize based on your needs)
        bool valid = true;
        for (const auto& pair : beforeHists) {
            if (!validateHistogramPair({pair.first, pair.second})) {
                valid = false;
                break;
            }
        }
        for (const auto& pair : afterHists) {
            if (!validateHistogramPair({pair.first, pair.second})) {
                valid = false;
                break;
            }
        }
        if (!valid) {
            std::cout << "Skipping bin " << binIndex << ": invalid histograms" << std::endl;
            continue;
        }

        canvas->Clear();
        canvas->Divide(1, 2); // Divide into two pads for before/after

        std::string binLabel = Form("%.2f #leq %s < %.2f", bins[binIndex]->min, binVariable_.c_str(), bins[binIndex]->max);
        TLegend* legendBefore = useLegend_ ? new TLegend(0.6, 0.7, 0.8, 0.9) : nullptr;
        TLegend* legendAfter = useLegend_ ? new TLegend(0.6, 0.7, 0.8, 0.9) : nullptr;

        canvas->cd(1);
        for (const auto& pair : beforeHists) {
            drawSingleHistogram(canvas, 1, {pair.first, pair.second}, legendBefore, FitResult(), binLabel);
        }
        if (legendBefore) {
            legendBefore->Draw();
        }

        canvas->cd(2);
        for (const auto& pair : afterHists) {
            drawSingleHistogram(canvas, 2, {pair.first, pair.second}, legendAfter, FitResult(), binLabel);
        }
        if (legendAfter) {
            legendAfter->Draw();
        }

        if (contaminationFunc) {
            contaminationFunc(canvas, binIndex, beforeHists, afterHists);
        }

        canvas->cd(0);
        canvas->Update();

        if (firstPage) {
            std::cout << "Opening PDF for first valid bin: " << binIndex << std::endl;
            canvas->Print(Form("%soverlap_%s(", pdfDir.c_str(), outputFile_.c_str()));
            firstPage = false;
        } else {
            std::cout << "Adding page for bin " << binIndex << std::endl;
        }

        std::string binName = Form("%.2f_%.2f", bins[binIndex]->min, bins[binIndex]->max);
        std::cout << "Bin " << binIndex << ": min=" << bins[binIndex]->min << ", max=" << bins[binIndex]->max << std::endl;
        std::cout << "Printing canvas for bin " << binIndex << " (Page " << ++pageCount << ")" << std::endl;
        canvas->Print(Form("%soverlap_%s", pdfDir.c_str(), outputFile_.c_str()));

        delete legendBefore;
        delete legendAfter;
    }

    if (pageCount > 0) {
        std::cout << "Closing PDF (Page " << pageCount << ")" << std::endl;
        canvas->Print(Form("%soverlap_%s)", pdfDir.c_str(), outputFile_.c_str()));
    }

    return canvas;
}