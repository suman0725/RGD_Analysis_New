#include <TApplication.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <cmath> // For fabs

void tmb() {
    // Define the z midpoints, deltaPt2 values, and errors
    double zRanges[] = {0.13, 0.23, 0.33, 0.43, 0.53, 0.63, 0.73, 0.83, 0.94}; // Midpoints of the z ranges
    double deltaPt2[] = {0.00407615, 0.00944956, 0.0101942, 0.0101942, 
                         0.0136208, 0.0107613, 0.0153118, 0.0115875, 0.0190322};
    double errors[] = {0.00359176, 0.0057764, 0.0109552, 0.0109552, 
                       0.0136737, 0.0154988, 0.0180084, 0.0183754, 0.0187732};

    // Number of points
    const int nPoints = sizeof(zRanges) / sizeof(zRanges[0]);

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
    TGraphErrors *graph = new TGraphErrors(nPoints, zRanges, deltaPt2, 0, errors);

    // Create a canvas
    TCanvas *canvas = new TCanvas("canvas", "#Delta pT^{2} vs z", 800, 600);

    // Set graph attributes
    graph->SetTitle("#DeltaP_{T}^{2} vs z; z; #DeltaP_{T}^{2} (GeV^{2})");
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
    legend->Draw(); */

    // Save the canvas as an image (optional)
    canvas->SaveAs("deltaPt2_vs_z.png");

    // Display the canvas
    canvas->Update();
}
