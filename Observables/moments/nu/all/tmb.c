#include <TApplication.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <cmath> // For fabs

void tmb() {
    // Define the nu ranges, ratio values, and errors
    double nuRanges[] = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0}; // Midpoints of the nu ranges
    double ratios[] = {0.994712, 1.0083, 0.999728, 1.01373, 
                       0.991925, 1.0025, 0.99498};
    double errors[] = {0.0167646, 0.0113713, 0.010157, 0.00955827, 
                       0.00870167, 0.00968034, 0.00874495};

    // Number of points
    const int nPoints = sizeof(nuRanges) / sizeof(nuRanges[0]);

    // Find Y-axis limits dynamically
    double minY = ratios[0] - errors[0];
    double maxY = ratios[0] + errors[0];
    for (int i = 1; i < nPoints; ++i) {
        if (ratios[i] - errors[i] < minY)
            minY = ratios[i] - errors[i];
        if (ratios[i] + errors[i] > maxY)
            maxY = ratios[i] + errors[i];
    }

    // Add some padding to the range for better visibility
    double padding = 0.1 * fabs(maxY - minY);
    minY -= padding;
    maxY += padding;

    // Create a TGraphErrors object
    TGraphErrors *graph = new TGraphErrors(nPoints, nuRanges, ratios, 0, errors);

    // Create a canvas
    TCanvas *canvas = new TCanvas("canvas", "Ratio vs #nu", 1900, 900);

    // Set graph attributes
    graph->SetTitle("Ratio of cos #phi moment vs #nu; #nu(GeV); #frac{<cos #phi_{h}>_{A}}{<cos #phi_{h}>_{LD2}}");
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);

    // Set Y-axis limits
    graph->GetYaxis()->SetRangeUser(minY, maxY);

    // Draw the graph
    graph->Draw("AP");

    canvas->SaveAs("ratio_vs_nu.png");

    // Display the canvas
    canvas->Update();
}
