#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLegend.h>

void conta() {
    // Midpoints of momentum bins
    double x[20] = {
        1.15, 1.45, 1.75, 2.05, 2.35, 2.65, 2.95, 3.25, 3.55, 3.85,
        4.15, 4.45, 4.75, 5.05, 5.35, 5.65, 5.95, 6.25, 6.55, 6.85
    };

    // Kaon-to-pion contamination (%)
    double y_kaon[20] = {
        1.70167e-14, 0.0104302, 0.0874723, 0.193868, 0.380622, 0.573549,
        0.618257, 0.946795, 1.60536, 5.58763, 7.5085, 8.46838,
        12.1787, 12.2555, 12.9358, 13.3861, 18.5198, 23.7642, 21.3723, 21.2439
    };

    // Proton-to-pion contamination (%)
    double y_proton[20] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // Below 4.15 GeV/c
        0.0916901, 0.18385, 0.191383, 0.327718, 0.169948, 0.414707,
        1.32798, 1.06248, 1.94156, 4.19467 // From 4.15 to 6.85 GeV/c
    };

    TCanvas *c = new TCanvas("c", "Kaon and Proton Contamination", 1200, 800);
    c->SetGrid(); // Grid on both axes

    // Kaon graph
    TGraph *gr_kaon = new TGraph(20, x, y_kaon);
    gr_kaon->SetMarkerStyle(20);     // Filled circles
    gr_kaon->SetMarkerSize(1.2);     // Same size as before
    gr_kaon->SetMarkerColor(kGreen); // Green for kaon
    gr_kaon->SetLineColor(kGreen);
    gr_kaon->SetName("Kaon");        // Name for legend

    // Proton graph
    TGraph *gr_proton = new TGraph(20, x, y_proton);
    gr_proton->SetMarkerStyle(20);     // Filled circles
    gr_proton->SetMarkerSize(1.2);     // Same size as kaon
    gr_proton->SetMarkerColor(kRed);   // Red for proton
    gr_proton->SetLineColor(kRed);
    gr_proton->SetName("Proton");      // Name for legend

    // Draw kaon graph first
    gr_kaon->Draw("AP");

    // Axis formatting
    gr_kaon->GetXaxis()->SetTitle("Momentum (GeV/c)");
    gr_kaon->GetYaxis()->SetTitle("Contamination (%)");
    gr_kaon->GetXaxis()->SetTitleSize(0.05);
    gr_kaon->GetYaxis()->SetTitleSize(0.05);
    gr_kaon->GetXaxis()->SetLabelSize(0.04);
    gr_kaon->GetYaxis()->SetLabelSize(0.04);
    gr_kaon->GetXaxis()->SetTitleOffset(1.0);
    gr_kaon->GetYaxis()->SetTitleOffset(1.0);

    // Axis ranges and divisions
    gr_kaon->GetYaxis()->SetRangeUser(0, 25);
    gr_kaon->GetYaxis()->SetNdivisions(25, kTRUE); // y-axis ticks every 1 unit
    gr_kaon->GetXaxis()->SetNdivisions(8);

    // Draw proton graph on same canvas
    gr_proton->Draw("P SAME");

    // Create and draw legend without box
    TLegend *legend = new TLegend(0.1, 0.8, 0.3, 0.9); // Position in top-right
    legend->AddEntry(gr_kaon, "Kaon", "p");
    legend->AddEntry(gr_proton, "Proton", "p");
    legend->SetBorderSize(0); // Remove box around legend
    legend->SetFillStyle(0);  // Make legend background transparent
    legend->Draw();

    // Save as PNG
    c->SaveAs("kaon_proton_contamination.png");
}