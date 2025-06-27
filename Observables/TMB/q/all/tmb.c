#include <TApplication.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <algorithm> // For std::min_element and std::max_element
#include <cmath>     // For fabs

void tmb() {
    // Define the Q2 midpoints, deltaPt2 values, and errors
    double Q2Ranges[] = {1.7, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.5}; // Midpoints of the Q2 ranges
    double deltaPt2[] = {0.00234206, 0.00264174, 0.00573692, 0.00714587, 
                         0.00993663, 0.0101643, 0.00448805, 0.00530052};
    double errors[] = {0.00320626, 0.0037374, 0.00464664, 0.00577216, 
                       0.00730564, 0.00896101, 0.0105102, 0.0105645};

    // Number of points
    const int nPoints = sizeof(Q2Ranges) / sizeof(Q2Ranges[0]);

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
    TGraphErrors *graph = new TGraphErrors(nPoints, Q2Ranges, deltaPt2, 0, errors);

    // Create a canvas
    TCanvas *canvas = new TCanvas("canvas", "#Delta pT^{2} vs Q^{2}", 800, 600);

    // Set graph attributes
    graph->SetTitle("#DeltaP_{T}^{2} vs Q^{2}; Q^{2} (GeV^{2}); #DeltaP_{T}^{2} (GeV^{2})");
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);

    // Draw the graph
    graph->Draw("AP");

    // Adjust the Y-axis range dynamically
    graph->GetYaxis()->SetRangeUser(minY, maxY);

    /* // Add a legend
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->SetHeader("Legend", "C");
    legend->AddEntry(graph, "Data points with errors", "lep");
    legend->Draw();
 */
    // Save the canvas as an image (optional)
    canvas->SaveAs("deltaPt2_vs_Q2.png");

    // Display the canvas
    canvas->Update();
}
