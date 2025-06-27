#ifndef LEGEND_UTILS_H
#define LEGEND_UTILS_H

#include <TLegend.h>
#include <TH1F.h>
#include <TF1.h>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

TLegend* addLegend(const std::vector<TH1F*>& histograms, const std::vector<std::string>& labels,
                   double x1 = 0.6, double y1 = 0.6, double x2 = 0.95, double y2 = 0.9,
                   const char* header = "") {
    if (histograms.size() != labels.size()) {
        std::cerr << "Error: Mismatch between histograms and labels." << std::endl;
        return nullptr;
    }

    TLegend* legend = new TLegend(x1, y1, x2, y2);
    legend->SetHeader(header);
    legend->SetBorderSize(1);
    legend->SetTextSize(0.015); // Smaller text to fit more content
    legend->SetFillStyle(0); // Transparent fill for better visibility

    for (size_t i = 0; i < histograms.size(); ++i) {
        if (histograms[i]) {
            // Add number of entries at the top
            legend->AddEntry(histograms[i], (std::stringstream() << "Entries: " << static_cast<int>(histograms[i]->GetEntries())).str().c_str(), "");

            // Gaussian fit
            if (TF1* gausFit = histograms[i]->GetFunction("gausFit")) {
                double chi2Ndf = gausFit->GetNDF() > 0 ? gausFit->GetChisquare() / gausFit->GetNDF() : 0.0;
                legend->AddEntry(histograms[i], labels[i].c_str(), "p"); // Marker for histogram
                legend->AddEntry(gausFit, "Gauss:", "l"); // Line for Gauss
                legend->AddEntry(histograms[i], (std::stringstream() << "#mu=" << std::fixed << std::setprecision(2) << gausFit->GetParameter(1)).str().c_str(), ""); // No indent
                legend->AddEntry(histograms[i], (std::stringstream() << "#sigma=" << std::setprecision(2) << gausFit->GetParameter(2)).str().c_str(), ""); // No indent
                legend->AddEntry(histograms[i], (std::stringstream() << "Constant=" << std::fixed << std::setprecision(2) << gausFit->GetParameter(0)).str().c_str(), ""); // No indent
                legend->AddEntry(histograms[i], (std::stringstream() << "#chi^{2}/NDF=" << std::setprecision(2) << gausFit->GetChisquare() << "/"
                                          << gausFit->GetNDF() << "=" << std::setprecision(2) << chi2Ndf).str().c_str(), ""); // No indent
            }

            // Crystal Ball fit
            if (TF1* cbFit = histograms[i]->GetFunction("cbFit")) {
                double chi2Ndf = cbFit->GetNDF() > 0 ? cbFit->GetChisquare() / cbFit->GetNDF() : 0.0;
                legend->AddEntry(cbFit, "CB:", "l"); // Line for CB
                legend->AddEntry(histograms[i], (std::stringstream() << "#mu=" << std::fixed << std::setprecision(2) << cbFit->GetParameter(1)).str().c_str(), ""); // No indent
                legend->AddEntry(histograms[i], (std::stringstream() << "#sigma=" << std::setprecision(2) << cbFit->GetParameter(2)).str().c_str(), ""); // No indent
                legend->AddEntry(histograms[i], (std::stringstream() << "Constant=" << std::fixed << std::setprecision(2) << cbFit->GetParameter(0)).str().c_str(), ""); // No indent
                legend->AddEntry(histograms[i], (std::stringstream() << "#alpha=" << std::setprecision(2) << cbFit->GetParameter(3)).str().c_str(), ""); // No indent
                legend->AddEntry(histograms[i], (std::stringstream() << "n=" << std::setprecision(2) << cbFit->GetParameter(4)).str().c_str(), ""); // No indent
                legend->AddEntry(histograms[i], (std::stringstream() << "#chi^{2}/NDF=" << std::setprecision(2) << cbFit->GetChisquare() << "/"
                                          << cbFit->GetNDF() << "=" << std::setprecision(2) << chi2Ndf).str().c_str(), ""); // No indent
            }
        }
    }

    return legend;
}

#endif // LEGEND_UTILS_H