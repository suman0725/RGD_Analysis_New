#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TApplication.h>

void plot_contamination() {
    // Define momentum bins
    const int nBins = 14; // 14 bin edges
    double pBins[nBins] = {1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 3.4, 3.7, 4.0, 4.3, 4.7, 5.0};
    double pCenters[nBins-1], pErrors[nBins-1];
    for (int i = 0; i < nBins-1; i++) {
        pCenters[i] = (pBins[i] + pBins[i+1]) / 2.0;
        pErrors[i] = (pBins[i+1] - pBins[i]) / 2.0; // Still calculated but not used
    }

    // Data for dt cut
    double kaon_dt[nBins-1] = {0, 0.870574, 0.735735, 0.544588, 1.05471, 1.59476, 2.55969, 5.35737, 8.5803, 11.9623, 11.4651, 11.2707, 8.3503};
    double kaon_dt_err[nBins-1] = {0, 0.005115, 0.00556963, 0.00566099, 0.00930031, 0.0133954, 0.0197241, 0.0330324, 0.0484279, 0.0665714, 0.0750051, 0.0763, 0.0903};
    double proton_dt[nBins-1] = {0, 0, 0, 0, 0, 1.16733, 1.05757, 1.27436, 1.41988, 2.30473, 1.39857, 2.8237, 5.9973};
    double proton_dt_err[nBins-1] = {0, 0, 0, 0, 0, 0.0114364, 0.0125851, 0.0157953, 0.0190395, 0.0279321, 0.0249857, 0.0367, 0.0757};

    // Data for chi2pid cut
    double kaon_chi2pid[nBins-1] = {0, 1.57622, 1.2919, 1.20897, 2.09853, 2.2527, 2.03686, 2.96657, 3.63595, 5.26463, 5.4968, 4.0294, 2.8295};
    double kaon_chi2pid_err[nBins-1] = {0, 0.00673, 0.00722554, 0.00825404, 0.0128403, 0.0156533, 0.0174219, 0.0246679, 0.0320352, 0.0454563, 0.054536, 0.0482, 0.0570};
    double proton_chi2pid[nBins-1] = {0, 0, 0, 0, 0, 1.82545, 1.58939, 1.7215, 1.81276, 2.81852, 2.71435, 2.6580, 5.8287};
    double proton_chi2pid_err[nBins-1] = {0, 0, 0, 0, 0, 0.0140615, 0.0153559, 0.0186773, 0.0224199, 0.0328712, 0.0378143, 0.0389, 0.0830};

    // Create TGraphErrors for each dataset with no x-axis errors
    TGraphErrors *gr_kaon_dt = new TGraphErrors(nBins-1, pCenters, kaon_dt, nullptr, kaon_dt_err);
    TGraphErrors *gr_proton_dt = new TGraphErrors(nBins-1, pCenters, proton_dt, nullptr, proton_dt_err);
    TGraphErrors *gr_kaon_chi2pid = new TGraphErrors(nBins-1, pCenters, kaon_chi2pid, nullptr, kaon_chi2pid_err);
    TGraphErrors *gr_proton_chi2pid = new TGraphErrors(nBins-1, pCenters, proton_chi2pid, nullptr, proton_chi2pid_err);

    // Set styles
    gr_kaon_dt->SetMarkerStyle(20); gr_kaon_dt->SetMarkerColor(kBlue); gr_kaon_dt->SetLineColor(kBlue);
    gr_proton_dt->SetMarkerStyle(21); gr_proton_dt->SetMarkerColor(kRed); gr_proton_dt->SetLineColor(kRed);
    gr_kaon_chi2pid->SetMarkerStyle(24); gr_kaon_chi2pid->SetMarkerColor(kBlue); gr_kaon_chi2pid->SetLineColor(kBlue);
    gr_proton_chi2pid->SetMarkerStyle(25); gr_proton_chi2pid->SetMarkerColor(kRed); gr_proton_chi2pid->SetLineColor(kRed);

    // Create canvases
    TCanvas *c1 = new TCanvas("c1", "Contamination after dt cut", 800, 600);
    TCanvas *c2 = new TCanvas("c2", "Contamination after chi2pid cut", 800, 600);
    TCanvas *c3 = new TCanvas("c3", "Kaon Contamination Comparison", 800, 600);
    TCanvas *c4 = new TCanvas("c4", "Proton Contamination Comparison", 800, 600);

    // Set global style
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // Plot 1: Contamination after dt cut
    c1->cd();
    gr_kaon_dt->SetTitle("Contamination after 3#sigma dt cut;Momentum (GeV/c);Contamination (%)");
    gr_kaon_dt->GetYaxis()->SetRangeUser(0, 15);
    gr_kaon_dt->Draw("AP");
    gr_proton_dt->Draw("P SAME");
    TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg1->AddEntry(gr_kaon_dt, "Kaon Contamination", "lep");
    leg1->AddEntry(gr_proton_dt, "Proton Contamination", "lep");
    leg1->Draw();

    // Plot 2: Contamination after chi2pid cut
    c2->cd();
    gr_kaon_chi2pid->SetTitle("Contamination after p-dependent #chi^{2} PID cut;Momentum (GeV/c);Contamination (%)");
    gr_kaon_chi2pid->GetYaxis()->SetRangeUser(0, 15);
    gr_kaon_chi2pid->Draw("AP");
    gr_proton_chi2pid->Draw("P SAME");
    TLegend *leg2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg2->AddEntry(gr_kaon_chi2pid, "Kaon Contamination", "lep");
    leg2->AddEntry(gr_proton_chi2pid, "Proton Contamination", "lep");
    leg2->Draw();

    // Plot 3: Kaon contamination comparison
    c3->cd();
    gr_kaon_dt->SetTitle("Kaon Contamination Comparison;Momentum (GeV/c);Contamination (%)");
    gr_kaon_dt->GetYaxis()->SetRangeUser(0, 15);
    gr_kaon_dt->Draw("AP");
    gr_kaon_chi2pid->Draw("P SAME");
    TLegend *leg3 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg3->AddEntry(gr_kaon_dt, "dt cut", "lep");
    leg3->AddEntry(gr_kaon_chi2pid, "#chi^{2} PID cut", "lep");
    leg3->Draw();

    // Plot 4: Proton contamination comparison
    c4->cd();
    gr_proton_dt->SetTitle("Proton Contamination Comparison;Momentum (GeV/c);Contamination (%)");
    gr_proton_dt->GetYaxis()->SetRangeUser(0, 15);
    gr_proton_dt->Draw("AP");
    gr_proton_chi2pid->Draw("P SAME");
    TLegend *leg4 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg4->AddEntry(gr_proton_dt, "dt cut", "lep");
    leg4->AddEntry(gr_proton_chi2pid, "#chi^{2} PID cut", "lep");
    leg4->Draw();

    // Save plots
    c1->SaveAs("contamination_dt.pdf");
    c2->SaveAs("contamination_chi2pid.pdf");
    c3->SaveAs("kaon_contamination_comparison.pdf");
    c4->SaveAs("proton_contamination_comparison.pdf");
}