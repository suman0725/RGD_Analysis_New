#include <TApplication.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <cmath> // For fabs

void tmb() {
    // Define the nu midpoints, deltaPt2 values, and errors
    double nuRanges[] = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0}; // Midpoints of the nu ranges
    double deltaPt2[] = {0.014605, 0.0114315, 0.00771431, 0.00131515, 
                         -0.00116962, 0.00547917, -0.0000119963};
    double errors[] = {0.00694052, 0.00550872, 0.00501, 0.00475021, 
                       0.00478253, 0.00481494, 0.00470132};

    // Number of points
    const int nPoints = sizeof(nuRanges) / sizeof(nuRanges[0]);

    // Find Y-axis limits dynamically
    double minY = deltaPt2[0] - errors[0];
    double maxY = deltaPt2[0] + errors[0];
    for (int i = 1; i < nPoints; ++i) {
        if (deltaPt2[i] - errors[i] < minY)
            minY = deltaPt2[i] - errors[i];
        if (deltaPt2[i] + errors[i] > maxY)
            maxY = deltaPt2[i] + errors[i];
    }

    // Add some padding to the range for better visibility
    double padding = 0.1 * fabs(maxY - minY);
    minY -= padding;
    maxY += padding;

    // Create a TGraphErrors object
    TGraphErrors *graph = new TGraphErrors(nPoints, nuRanges, deltaPt2, 0, errors);

    // Create a canvas
    TCanvas *canvas = new TCanvas("canvas", "#Delta pT^{2} vs #nu", 800, 600);

    // Set graph attributes
    graph->SetTitle("#DeltaP_{T}^{2} vs #nu; #nu (GeV); #DeltaP_{T}^{2} (GeV^{2})");
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kRed);
    graph->SetLineColor(kRed);

    // Draw the graph
    graph->Draw("AP");

    // Adjust the Y-axis range dynamically
    graph->GetYaxis()->SetRangeUser(minY, maxY);

    /* // Add a legend
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->SetHeader("Legend", "C");
    legend->AddEntry(graph, "Data points with errors", "lep");
    legend->Draw(); */

    // Save the canvas as an image (optional)
    canvas->SaveAs("deltaPt2_vs_nu.png");

    // Display the canvas
    canvas->Update();
}
