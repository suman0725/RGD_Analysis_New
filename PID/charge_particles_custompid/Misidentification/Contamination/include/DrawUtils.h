#ifndef DRAW_UTILS_H
#define DRAW_UTILS_H

#include <TH1F.h>
#include <TF1.h>
#include <TLine.h>
#include <vector>

// Reusable function to draw mean, mean - 3σ, and mean + 3σ lines for Gaussian and CB fits
std::vector<TLine*> drawMeanSigmaLines(TH1F* hist, double yMin, double yMax,
                                       int gausColor = kBlue+2, int cbColor = kRed+2,
                                       int lineStyle = kDashed) {
    std::vector<TLine*> lines;
    if (!hist) return lines;

    // Use histogram's actual maximum for yMax if provided yMax is larger
    double effectiveYMax = yMax;
    if (hist->GetMaximum() > 0 && hist->GetMaximum() < yMax) {
        effectiveYMax = hist->GetMaximum() * 1.05; // Slightly above max for visibility
    }

    // Gaussian fit lines
    if (TF1* gausFit = hist->GetFunction("gausFit")) {
        double mean = gausFit->GetParameter(1);
        double sigma = gausFit->GetParameter(2);
        double xLow = mean - 3 * sigma;
        double xHigh = mean + 3 * sigma;

        // Draw line at mean - 3σ
        TLine* lineLowGaus = new TLine(xLow, yMin, xLow, effectiveYMax);
        lineLowGaus->SetLineColor(gausColor);
        lineLowGaus->SetLineStyle(lineStyle);
        lineLowGaus->Draw();
        lines.push_back(lineLowGaus);

        // Draw line at mean
        TLine* lineMeanGaus = new TLine(mean, yMin, mean, effectiveYMax);
        lineMeanGaus->SetLineColor(gausColor);
        lineMeanGaus->SetLineStyle(kSolid); // Solid line for mean to distinguish
        lineMeanGaus->Draw();
        lines.push_back(lineMeanGaus);

        // Draw line at mean + 3σ
        TLine* lineHighGaus = new TLine(xHigh, yMin, xHigh, effectiveYMax);
        lineHighGaus->SetLineColor(gausColor);
        lineHighGaus->SetLineStyle(lineStyle);
        lineHighGaus->Draw();
        lines.push_back(lineHighGaus);
    }

    // Crystal Ball fit lines
    if (TF1* cbFit = hist->GetFunction("cbFit")) {
        double mean = cbFit->GetParameter(1);
        double sigma = cbFit->GetParameter(2);
        double xLow = mean - 3 * sigma;
        double xHigh = mean + 3 * sigma;

        // Draw line at mean - 3σ
        TLine* lineLowCB = new TLine(xLow, yMin, xLow, effectiveYMax);
        lineLowCB->SetLineColor(cbColor);
        lineLowCB->SetLineStyle(lineStyle);
        lineLowCB->Draw();
        lines.push_back(lineLowCB);

        // Draw line at mean
        TLine* lineMeanCB = new TLine(mean, yMin, mean, effectiveYMax);
        lineMeanCB->SetLineColor(cbColor);
        lineMeanCB->SetLineStyle(kSolid); // Solid line for mean to distinguish
        lineMeanCB->Draw();
        lines.push_back(lineMeanCB);

        // Draw line at mean + 3σ
        TLine* lineHighCB = new TLine(xHigh, yMin, xHigh, effectiveYMax);
        lineHighCB->SetLineColor(cbColor);
        lineHighCB->SetLineStyle(lineStyle);
        lineHighCB->Draw();
        lines.push_back(lineHighCB);
    }

    return lines;
}

#endif // DRAW_UTILS_H