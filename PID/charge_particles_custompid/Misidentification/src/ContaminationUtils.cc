#include "ContaminationUtils.h"
#include <TF1.h>

double calculateContamination(TH1F* histPions, TH1F* histKaons, double min, double max,
                              TCanvas* canvas, const char* binName, TLegend* legend) {
    if (!histPions || !histKaons) {
        throw std::invalid_argument("Invalid histogram pointers");
    }

    // Clone histograms to avoid modifying originals
    TH1F* hPions = (TH1F*)histPions->Clone(Form("hPions_%s", binName));
    TH1F* hKaons = (TH1F*)histKaons->Clone(Form("hKaons_%s", binName));

    // Fit Gaussians
    TF1* gausPion = new TF1(Form("gaus_pion_%s", binName), "gaus", min, max);
    TF1* gausKaon = new TF1(Form("gaus_kaon_%s", binName), "gaus", min, max);

    gausPion->SetParameters(hPions->GetMaximum(), hPions->GetMean(), hPions->GetRMS());
    gausKaon->SetParameters(hKaons->GetMaximum(), hKaons->GetMean(), hKaons->GetRMS());

    hPions->Fit(gausPion, "RQN");
    hKaons->Fit(gausKaon, "RQN");

    // Calculate integrals
    double pionIntegral = gausPion->Integral(min, max);
    double kaonIntegral = gausKaon->Integral(min, max);
    double totalIntegral = pionIntegral + kaonIntegral;

    // Calculate contamination
    double contamination = (totalIntegral > 0) ? kaonIntegral / totalIntegral : 0.0;

    // Update canvas if provided
    if (canvas) {
        hPions->SetLineColor(kBlue);
        hKaons->SetLineColor(kGreen);
        hPions->SetFillColorAlpha(kBlue, 0.3);
        hKaons->SetFillColorAlpha(kGreen, 0.3);

        hPions->Draw("HIST");
        hKaons->Draw("HIST SAME");

        gausPion->SetLineColor(kBlue);
        gausKaon->SetLineColor(kGreen);
        gausPion->Draw("SAME");
        gausKaon->Draw("SAME");

        if (legend) {
            legend->AddEntry(hPions, "Pions", "lf");
            legend->AddEntry(hKaons, "Kaons", "lf");
            legend->AddEntry(gausPion, "Pion Fit", "l");
            legend->AddEntry(gausKaon, "Kaon Fit", "l");
            legend->AddEntry((TObject*)0, Form("Contamination: %.2f%%", contamination * 100), "");
        }
    }

    delete gausPion;
    delete gausKaon;
    delete hPions;
    delete hKaons;

    return contamination;
}