#include "include/BinUtils.h"
#include "include/HistogramUtils.h"
#include "include/FitUtils.h"
#include "include/CanvasUtils.h"
#include "include/LegendUtils.h"
#include "include/DrawUtils.h"
#include "include/ContaminationUtils.h"
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TString.h>
#include <TF1.h>
#include <TLine.h>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "TLatex.h"


// Helper function to create output directories
void createoutputDirs() {
    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/pdf", kTRUE);
    gSystem->mkdir("output/root", kTRUE);
}

int main() {
    // Load ROOT file
    TFile* file = new TFile("../pkptreeCxC_4.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file ../pkptreeCxC_4.root" << std::endl;
        return 1;
    }

    TTree* treePions = (TTree*)file->Get("EB_pid_pions");
    TTree* treeKaons = (TTree*)file->Get("EB_pid_kaons");
    if (!treePions || !treeKaons) {
        std::cerr << "Error: Cannot load trees EB_pid_pions or EB_pid_kaons" << std::endl;
        file->Close();
        delete file;
        return 1;
    }

    // Create output directories
    createoutputDirs();

    system("mkdir -p output/csv");
    std::ofstream csvFile("output/csv/contamination.csv");
    csvFile << "Momentum Bin (GeV/c),c1,c2,Contamination (%)\n";

    // Open CSV file for Task 1 fit parameters (Gaussian and Crystal Ball)
    std::ofstream csvTask1("output/task1_gauscbfitparameter.csv");
    csvTask1 << "Momentum Bin (GeV/c),Pion Gaus Mean,Pion Gaus Sigma,Pion Gaus Mean-3Sigma,Pion Gaus Mean+3Sigma,"
             << "Pion CB Mean,Pion CB Sigma,Pion CB Mean-3Sigma,Pion CB Mean+3Sigma,"
             << "Kaon Gaus Mean,Kaon Gaus Sigma,Kaon Gaus Mean-3Sigma,Kaon Gaus Mean+3Sigma,"
             << "Kaon CB Mean,Kaon CB Sigma,Kaon CB Mean-3Sigma,Kaon CB Mean+3Sigma\n";

    // Open CSV file for Task 1 Gaussian-only fit parameters
    std::ofstream csvTask1Gaus("output/task1_gaus_only_fitparameter.csv");
    csvTask1Gaus << "Momentum Bin (GeV/c),Pion Gaus Mean,Pion Gaus Sigma,Pion Gaus Mean-3Sigma,Pion Gaus Mean+3Sigma,"
                 << "Kaon Gaus Mean,Kaon Gaus Sigma,Kaon Gaus Mean-3Sigma,Kaon Gaus Mean+3Sigma\n";

    // Task 1: 1D chi2pid histograms in p bins with Gaussian and Crystal Ball fits
    // Define momentum bins
    std::vector<double> pBins = createBins(1, 10, 30); // 10 linear bins from 1 to 4 GeV/c
    if (pBins.empty()) {
        std::cerr << "Error: Failed to create p bins." << std::endl;
        file->Close();
        delete file;
        return 1;
    }

    // Open PDF for Task 1 (original)
    TCanvas* pdfCanvas = new TCanvas("pdfCanvas", "pdfCanvas", 1600, 800);
    pdfCanvas->SetCanvasSize(1600, 800);
    pdfCanvas->Print("output/pdf/task1_chi2pid.pdf[");

    // Open PDF for Task 1 (Gaussian only, new)
    TCanvas* pdfCanvasGaus = new TCanvas("pdfCanvasGaus", "pdfCanvasGaus", 1600, 800);
    pdfCanvasGaus->SetCanvasSize(1600, 800);
    pdfCanvasGaus->Print("output/pdf/task1_chi2pid_gaus_only.pdf[");

    // Store Gaussian and CB fit parameters for chi2pid (needed for cuts in Task 2)
    std::map<std::string, std::vector<double>> fitParams;
    fitParams["pions_gaus_mean"] = std::vector<double>(pBins.size() - 1, 0);
    fitParams["pions_gaus_sigma"] = std::vector<double>(pBins.size() - 1, 0);
    fitParams["pions_cb_mean"] = std::vector<double>(pBins.size() - 1, 0);
    fitParams["pions_cb_sigma"] = std::vector<double>(pBins.size() - 1, 0);
    fitParams["kaons_gaus_mean"] = std::vector<double>(pBins.size() - 1, 0);
    fitParams["kaons_gaus_sigma"] = std::vector<double>(pBins.size() - 1, 0);
    fitParams["kaons_cb_mean"] = std::vector<double>(pBins.size() - 1, 0);
    fitParams["kaons_cb_sigma"] = std::vector<double>(pBins.size() - 1, 0);

    // Store Gaussian-only fit parameters for the new PDF
    std::map<std::string, std::vector<double>> fitParamsGaus;
    fitParamsGaus["pions_gaus_mean"] = std::vector<double>(pBins.size() - 1, 0);
    fitParamsGaus["pions_gaus_sigma"] = std::vector<double>(pBins.size() - 1, 0);
    fitParamsGaus["kaons_gaus_mean"] = std::vector<double>(pBins.size() - 1, 0);
    fitParamsGaus["kaons_gaus_sigma"] = std::vector<double>(pBins.size() - 1, 0);

    // Process each p bin for Task 1
    for (size_t i = 0; i < pBins.size() - 1; ++i) {
        // Define p cut
        TString pCut = TString::Format("p >= %f && p < %f", pBins[i], pBins[i+1]);

        // Create chi2pid histograms for original Task 1
        TH1F* histPions = createHistogram(treePions, "chi2pid", {},
                                          TString::Format("p: [%.2f-%.2f) GeV/c", pBins[i], pBins[i+1]),
                                          "chi2pid", "Counts", kBlue, kBlue, 3001, "hist_pions",
                                          -10.0, 10.0);
        TH1F* histKaons = createHistogram(treeKaons, "chi2pid", {},
                                          TString::Format("p: [%.2f-%.2f) GeV/c", pBins[i], pBins[i+1]),
                                          "chi2pid", "Counts", kGreen, kGreen, 3002, "hist_kaons",
                                          -10.0, 10.0);
        if (i > 10) {
            histKaons->GetXaxis()->SetRangeUser(-5, 4);
        }

        // Create chi2pid histograms for Gaussian-only PDF
        TH1F* histPionsGaus = createHistogram(treePions, "chi2pid", {},
                                              TString::Format("p: [%.2f-%.2f) GeV/c", pBins[i], pBins[i+1]),
                                              "chi2pid", "Counts", kBlue, kBlue, 3001, "hist_pions_gaus",
                                              -4.0, 4.0);
        TH1F* histKaonsGaus = createHistogram(treeKaons, "chi2pid", {},
                                              TString::Format("p: [%.2f-%.2f) GeV/c", pBins[i], pBins[i+1]),
                                              "chi2pid", "Counts", kGreen, kGreen, 3002, "hist_kaons_gaus",
                                              -4.0, 4.0);

        // Apply p cut for original histograms
        if (histPions) {
            histPions->Reset();
            treePions->Draw(TString::Format("chi2pid>>%s", histPions->GetName()), pCut, "goff");
            std::cout << "Pions histogram (bin " << i << ") entries: " << histPions->GetEntries() << std::endl;
        }
        if (histKaons) {
            histKaons->Reset();
            treeKaons->Draw(TString::Format("chi2pid>>%s", histKaons->GetName()), pCut, "goff");
            std::cout << "Kaons histogram (bin " << i << ") entries: " << histKaons->GetEntries() << std::endl;
        }

        // Apply p cut for Gaussian-only histograms
        if (histPionsGaus) {
            histPionsGaus->Reset();
            treePions->Draw(TString::Format("chi2pid>>%s", histPionsGaus->GetName()), pCut, "goff");
            std::cout << "Pions Gaussian-only histogram (bin " << i << ") entries: " << histPionsGaus->GetEntries() << std::endl;
        }
        if (histKaonsGaus) {
            histKaonsGaus->Reset();
            treeKaons->Draw(TString::Format("chi2pid>>%s", histKaonsGaus->GetName()), pCut, "goff");
            std::cout << "Kaons Gaussian-only histogram (bin " << i << ") entries: " << histKaonsGaus->GetEntries() << std::endl;
        }

        // Fit histograms with both Gaussian and Crystal Ball for original Task 1
        std::map<std::string, double> fitPions, fitKaons;
        if (histPions && histPions->GetEntries() > 50) {
            fitPions = fitHistogram(histPions, "fit_pions", -15, 15, "both");
            fitParams["pions_gaus_mean"][i] = fitPions["gaus_mean"];
            fitParams["pions_gaus_sigma"][i] = fitPions["gaus_sigma"];
            fitParams["pions_cb_mean"][i] = fitPions["cb_mean"];
            fitParams["pions_cb_sigma"][i] = fitPions["cb_sigma"];
        } else {
            std::cerr << "Skipping fit for pions (bin " << i << ") due to insufficient entries." << std::endl;
        }
        if (histKaons && histKaons->GetEntries() > 50) {
            if (i >= 8) {
                fitKaons = fitHistogram(histKaons, "fit_kaons", -1, 5, "both");
            } else {
                fitKaons = fitHistogram(histKaons, "fit_kaons", -10, 15, "both");
            }
            fitParams["kaons_gaus_mean"][i] = fitKaons["gaus_mean"];
            fitParams["kaons_gaus_sigma"][i] = fitKaons["gaus_sigma"];
            fitParams["kaons_cb_mean"][i] = fitKaons["cb_mean"];
            fitParams["kaons_cb_sigma"][i] = fitKaons["cb_sigma"];
        } else {
            std::cerr << "Skipping fit for kaons (bin " << i << ") due to insufficient entries." << std::endl;
        }

                // Fit histograms with Gaussian only for new PDF with adjusted range
                std::map<std::string, double> fitPionsGaus, fitKaonsGaus;
                if (histPionsGaus && histPionsGaus->GetEntries() > 50) {
                    TF1* gausFitPions = new TF1("gausFit", "gaus", -1.5, 1.5);
                    // Set initial parameters to help convergence (optional)
                    gausFitPions->SetParameters(histPionsGaus->GetMaximum(), 0.0, 1.0);
                    histPionsGaus->Fit(gausFitPions, "R", "", -1.5, 1.5);
                    fitPionsGaus["gaus_mean"] = gausFitPions->GetParameter(1);
                    fitPionsGaus["gaus_sigma"] = gausFitPions->GetParameter(2);
                    fitPionsGaus["gaus_constant"] = gausFitPions->GetParameter(0);
                    fitParamsGaus["pions_gaus_mean"][i] = fitPionsGaus["gaus_mean"];
                    fitParamsGaus["pions_gaus_sigma"][i] = fitPionsGaus["gaus_sigma"];
                } else {
                    std::cerr << "Skipping Gaussian fit for pions (bin " << i << ") due to insufficient entries." << std::endl;
                }
                if (histKaonsGaus && histKaonsGaus->GetEntries() > 50) {
                    TF1* gausFitKaons = new TF1("gausFit", "gaus", -1.5, 1.5);
                    gausFitKaons->SetParameters(histKaonsGaus->GetMaximum(), 0.0, 1.0);
                    histKaonsGaus->Fit(gausFitKaons, "R", "", -1.5, 1.5);
                    fitKaonsGaus["gaus_mean"] = gausFitKaons->GetParameter(1);
                    fitKaonsGaus["gaus_sigma"] = gausFitKaons->GetParameter(2);
                    fitKaonsGaus["gaus_constant"] = gausFitKaons->GetParameter(0);
                    fitParamsGaus["kaons_gaus_mean"][i] = fitKaonsGaus["gaus_mean"];
                    fitParamsGaus["kaons_gaus_sigma"][i] = fitKaonsGaus["gaus_sigma"];
                } else {
                    std::cerr << "Skipping Gaussian fit for kaons (bin " << i << ") due to insufficient entries." << std::endl;
                }

        // Write fit parameters to CSV for original Task 1
        csvTask1 << pBins[i] << "-" << pBins[i+1] << ",";
        csvTask1 << std::fixed << std::setprecision(6) << fitParams["pions_gaus_mean"][i] << ","
                 << fitParams["pions_gaus_sigma"][i] << ","
                 << (fitParams["pions_gaus_mean"][i] - 3 * fitParams["pions_gaus_sigma"][i]) << ","
                 << (fitParams["pions_gaus_mean"][i] + 3 * fitParams["pions_gaus_sigma"][i]) << ","
                 << fitParams["pions_cb_mean"][i] << ","
                 << fitParams["pions_cb_sigma"][i] << ","
                 << (fitParams["pions_cb_mean"][i] - 3 * fitParams["pions_cb_sigma"][i]) << ","
                 << (fitParams["pions_cb_mean"][i] + 3 * fitParams["pions_cb_sigma"][i]) << ","
                 << fitParams["kaons_gaus_mean"][i] << ","
                 << fitParams["kaons_gaus_sigma"][i] << ","
                 << (fitParams["kaons_gaus_mean"][i] - 3 * fitParams["kaons_gaus_sigma"][i]) << ","
                 << (fitParams["kaons_gaus_mean"][i] + 3 * fitParams["kaons_gaus_sigma"][i]) << ","
                 << fitParams["kaons_cb_mean"][i] << ","
                 << fitParams["kaons_cb_sigma"][i] << ","
                 << (fitParams["kaons_cb_mean"][i] - 3 * fitParams["kaons_cb_sigma"][i]) << ","
                 << (fitParams["kaons_cb_mean"][i] + 3 * fitParams["kaons_cb_sigma"][i]) << "\n";

        // Write fit parameters to CSV for Gaussian-only PDF
        csvTask1Gaus << pBins[i] << "-" << pBins[i+1] << ",";
        csvTask1Gaus << std::fixed << std::setprecision(6) << fitParamsGaus["pions_gaus_mean"][i] << ","
                     << fitParamsGaus["pions_gaus_sigma"][i] << ","
                     << (fitParamsGaus["pions_gaus_mean"][i] - 3 * fitParamsGaus["pions_gaus_sigma"][i]) << ","
                     << (fitParamsGaus["pions_gaus_mean"][i] + 3 * fitParamsGaus["pions_gaus_sigma"][i]) << ","
                     << fitParamsGaus["kaons_gaus_mean"][i] << ","
                     << fitParamsGaus["kaons_gaus_sigma"][i] << ","
                     << (fitParamsGaus["kaons_gaus_mean"][i] - 3 * fitParamsGaus["kaons_gaus_sigma"][i]) << ","
                     << (fitParamsGaus["kaons_gaus_mean"][i] + 3 * fitParamsGaus["kaons_gaus_sigma"][i]) << "\n";

        // Create canvas for original Task 1 (1x2: pions left, kaons right)
        std::vector<TH1F*> histograms = {histPions, histKaons};
        TCanvas* canvas = createCanvas(histograms, 1, 2, {1, 2},
                                      TString::Format("canvas_task1_bin%zu", i), 1600, 800);

        // Add legend, draw fits, and mean ± 3σ lines for original Task 1
        if (canvas) {
            canvas->cd(1);
            if (histPions) {
                histPions->Draw("HIST");
                if (histPions->GetFunction("gausFit")) {
                    histPions->GetFunction("gausFit")->SetLineColor(kBlue+2);
                    histPions->GetFunction("gausFit")->SetLineStyle(1);
                    histPions->GetFunction("gausFit")->SetLineWidth(2);
                    histPions->GetFunction("gausFit")->Draw("SAME");
                }
                if (histPions->GetFunction("cbFit")) {
                    histPions->GetFunction("cbFit")->SetLineColor(kRed+2);
                    histPions->GetFunction("cbFit")->SetLineStyle(1);
                    histPions->GetFunction("cbFit")->SetLineWidth(2);
                    histPions->GetFunction("cbFit")->Draw("SAME");
                }
                double yMaxPions = histPions->GetMaximum() * 1.1;
                std::vector<TLine*> linesPions = drawMeanSigmaLines(histPions, 0, yMaxPions);
                std::vector<TObject*> legendObjects = {histPions};
                std::vector<std::string> legendLabels = {"Pions"};
                if (linesPions.size() >= 3) {
                    legendObjects.push_back(linesPions[1]);
                    legendLabels.push_back("Gaus #mu");
                    legendObjects.push_back(linesPions[0]);
                    legendLabels.push_back("Gaus #mu #pm 3#sigma");
                }
                if (linesPions.size() >= 6) {
                    legendObjects.push_back(linesPions[4]);
                    legendLabels.push_back("CB #mu");
                    legendObjects.push_back(linesPions[3]);
                    legendLabels.push_back("CB #mu #pm 3#sigma");
                }
                TLegend* legPions = addLegend({histPions}, {"Pions"}, 0.7, 0.6, 0.9, 0.9, "");
                if (legPions) {
                    legPions->Clear();
                    for (size_t j = 0; j < legendObjects.size(); ++j) {
                        legPions->AddEntry(legendObjects[j], legendLabels[j].c_str(), j == 0 ? "p" : "l");
                    }
                    if (TF1* gausFit = histPions->GetFunction("gausFit")) {
                        double chi2Ndf = gausFit->GetNDF() > 0 ? gausFit->GetChisquare() / gausFit->GetNDF() : 0.0;
                        legPions->AddEntry((TObject*)0, (std::stringstream() << "Entries: " << static_cast<int>(histPions->GetEntries())).str().c_str(), "");
                        legPions->AddEntry(gausFit, "Gauss:", "l");
                        legPions->AddEntry((TObject*)0, (std::stringstream() << "Constant=" << std::fixed << std::setprecision(2) << gausFit->GetParameter(0)).str().c_str(), "");
                        legPions->AddEntry((TObject*)0, (std::stringstream() << "#mu=" << std::fixed << std::setprecision(3) << gausFit->GetParameter(1)).str().c_str(), "");
                        legPions->AddEntry((TObject*)0, (std::stringstream() << "#sigma=" << std::setprecision(3) << gausFit->GetParameter(2)).str().c_str(), "");
                        legPions->AddEntry((TObject*)0, (std::stringstream() << "#chi^{2}/NDF=" << std::setprecision(2) << chi2Ndf).str().c_str(), "");
                    }
                    if (TF1* cbFit = histPions->GetFunction("cbFit")) {
                        double chi2Ndf = cbFit->GetNDF() > 0 ? cbFit->GetChisquare() / cbFit->GetNDF() : 0.0;
                        legPions->AddEntry(cbFit, "CB:", "l");
                        legPions->AddEntry((TObject*)0, (std::stringstream() << "Constant=" << std::fixed << std::setprecision(2) << cbFit->GetParameter(0)).str().c_str(), "");
                        legPions->AddEntry((TObject*)0, (std::stringstream() << "#mu=" << std::fixed << std::setprecision(3) << cbFit->GetParameter(1)).str().c_str(), "");
                        legPions->AddEntry((TObject*)0, (std::stringstream() << "#sigma=" << std::setprecision(3) << cbFit->GetParameter(2)).str().c_str(), "");
                        legPions->AddEntry((TObject*)0, (std::stringstream() << "#alpha=" << std::setprecision(2) << cbFit->GetParameter(3)).str().c_str(), "");
                        legPions->AddEntry((TObject*)0, (std::stringstream() << "n=" << std::setprecision(2) << cbFit->GetParameter(4)).str().c_str(), "");
                        legPions->AddEntry((TObject*)0, (std::stringstream() << "#chi^{2}/NDF=" << std::setprecision(2) << chi2Ndf).str().c_str(), "");
                    }
                    legPions->Draw();
                }
            }
            canvas->cd(2);
            if (histKaons) {
                histKaons->Draw("HIST");
                if (histKaons->GetFunction("gausFit")) {
                    histKaons->GetFunction("gausFit")->SetLineColor(kBlue+2);
                    histKaons->GetFunction("gausFit")->SetLineStyle(1);
                    histKaons->GetFunction("gausFit")->SetLineWidth(2);
                    histKaons->GetFunction("gausFit")->Draw("SAME");
                }
                if (histKaons->GetFunction("cbFit")) {
                    histKaons->GetFunction("cbFit")->SetLineColor(kRed+2);
                    histKaons->GetFunction("cbFit")->SetLineStyle(1);
                    histKaons->GetFunction("cbFit")->SetLineWidth(2);
                    histKaons->GetFunction("cbFit")->Draw("SAME");
                }
                double yMaxKaons = histKaons->GetMaximum() * 1.1;
                std::vector<TLine*> linesKaons = drawMeanSigmaLines(histKaons, 0, yMaxKaons);
                std::vector<TObject*> legendObjects = {histKaons};
                std::vector<std::string> legendLabels = {"Kaons"};
                if (linesKaons.size() >= 3) {
                    legendObjects.push_back(linesKaons[1]);
                    legendLabels.push_back("Gaus #mu");
                    legendObjects.push_back(linesKaons[0]);
                    legendLabels.push_back("Gaus #mu #pm 3#sigma");
                }
                if (linesKaons.size() >= 6) {
                    legendObjects.push_back(linesKaons[4]);
                    legendLabels.push_back("CB #mu");
                    legendObjects.push_back(linesKaons[3]);
                    legendLabels.push_back("CB #mu #pm 3#sigma");
                }
                TLegend* legKaons = addLegend({histKaons}, {"Kaons"}, 0.7, 0.6, 0.9, 0.9, "");
                if (legKaons) {
                    legKaons->Clear();
                    for (size_t j = 0; j < legendObjects.size(); ++j) {
                        legKaons->AddEntry(legendObjects[j], legendLabels[j].c_str(), j == 0 ? "p" : "l");
                    }
                    if (TF1* gausFit = histKaons->GetFunction("gausFit")) {
                        double chi2Ndf = gausFit->GetNDF() > 0 ? gausFit->GetChisquare() / gausFit->GetNDF() : 0.0;
                        legKaons->AddEntry((TObject*)0, (std::stringstream() << "Entries: " << static_cast<int>(histKaons->GetEntries())).str().c_str(), "");
                        legKaons->AddEntry(gausFit, "Gauss:", "l");
                        legKaons->AddEntry((TObject*)0, (std::stringstream() << "Constant=" << std::fixed << std::setprecision(2) << gausFit->GetParameter(0)).str().c_str(), "");
                        legKaons->AddEntry((TObject*)0, (std::stringstream() << "#mu=" << std::fixed << std::setprecision(3) << gausFit->GetParameter(1)).str().c_str(), "");
                        legKaons->AddEntry((TObject*)0, (std::stringstream() << "#sigma=" << std::setprecision(3) << gausFit->GetParameter(2)).str().c_str(), "");
                        legKaons->AddEntry((TObject*)0, (std::stringstream() << "#chi^{2}/NDF=" << std::setprecision(2) << chi2Ndf).str().c_str(), "");
                    }
                    if (TF1* cbFit = histKaons->GetFunction("cbFit")) {
                        double chi2Ndf = cbFit->GetNDF() > 0 ? cbFit->GetChisquare() / cbFit->GetNDF() : 0.0;
                        legKaons->AddEntry(cbFit, "CB:", "l");
                        legKaons->AddEntry((TObject*)0, (std::stringstream() << "Constant=" << std::fixed << std::setprecision(2) << cbFit->GetParameter(0)).str().c_str(), "");
                        legKaons->AddEntry((TObject*)0, (std::stringstream() << "#mu=" << std::fixed << std::setprecision(3) << cbFit->GetParameter(1)).str().c_str(), "");
                        legKaons->AddEntry((TObject*)0, (std::stringstream() << "#sigma=" << std::setprecision(3) << cbFit->GetParameter(2)).str().c_str(), "");
                        legKaons->AddEntry((TObject*)0, (std::stringstream() << "#alpha=" << std::setprecision(2) << cbFit->GetParameter(3)).str().c_str(), "");
                        legKaons->AddEntry((TObject*)0, (std::stringstream() << "n=" << std::setprecision(2) << cbFit->GetParameter(4)).str().c_str(), "");
                        legKaons->AddEntry((TObject*)0, (std::stringstream() << "#chi^{2}/NDF=" << std::setprecision(2) << chi2Ndf).str().c_str(), "");
                    }
                    legKaons->Draw();
                }
            }

            // Save to PDF
            canvas->Print("output/pdf/task1_chi2pid.pdf");

            // Save to ROOT file
            TFile* outRoot = new TFile(TString::Format("output/root/task1_bin%zu.root", i), "RECREATE");
            canvas->Write();
            outRoot->Close();
            delete outRoot;
        }

        // Create canvas for Gaussian-only PDF (1x2: pions left, kaons right)
        std::vector<TH1F*> histogramsGaus = {histPionsGaus, histKaonsGaus};
        TCanvas* canvasGaus = createCanvas(histogramsGaus, 1, 2, {1, 2},
                                          TString::Format("canvas_task1_gaus_bin%zu", i), 1600, 800);

        // Add legend, draw fits, and mean ± 3σ lines for Gaussian-only PDF
        if (canvasGaus) {
            canvasGaus->cd(1);
            if (histPionsGaus) {
                histPionsGaus->Draw("HIST");
                if (histPionsGaus->GetFunction("gausFit")) {
                    histPionsGaus->GetFunction("gausFit")->SetLineColor(kRed+2);
                    histPionsGaus->GetFunction("gausFit")->SetLineStyle(1);
                    histPionsGaus->GetFunction("gausFit")->SetLineWidth(2);
                    histPionsGaus->GetFunction("gausFit")->Draw("SAME");
                }
                double yMaxPions = histPionsGaus->GetMaximum() * 1.1;
                std::vector<TLine*> linesPions = drawMeanSigmaLines(histPionsGaus, 0, yMaxPions);
                std::vector<TObject*> legendObjects = {histPionsGaus};
                std::vector<std::string> legendLabels = {"Pions"};
                if (linesPions.size() >= 3) {
                    legendObjects.push_back(linesPions[1]);
                    legendLabels.push_back("Gaus #mu");
                    legendObjects.push_back(linesPions[0]);
                    legendLabels.push_back("Gaus #mu #pm 3#sigma");
                }
                TLegend* legPions = addLegend({histPionsGaus}, {"Pions"}, 0.7, 0.6, 0.9, 0.9, "");
                if (legPions) {
                    legPions->Clear();
                    for (size_t j = 0; j < legendObjects.size(); ++j) {
                        legPions->AddEntry(legendObjects[j], legendLabels[j].c_str(), j == 0 ? "p" : "l");
                    }
                    if (TF1* gausFit = histPionsGaus->GetFunction("gausFit")) {
                        double chi2Ndf = gausFit->GetNDF() > 0 ? gausFit->GetChisquare() / gausFit->GetNDF() : 0.0;
                        legPions->AddEntry((TObject*)0, (std::stringstream() << "Entries: " << static_cast<int>(histPionsGaus->GetEntries())).str().c_str(), "");
                        legPions->AddEntry(gausFit, "Gauss:", "l");
                        legPions->AddEntry((TObject*)0, (std::stringstream() << "Constant=" << std::fixed << std::setprecision(2) << gausFit->GetParameter(0)).str().c_str(), "");
                        legPions->AddEntry((TObject*)0, (std::stringstream() << "#mu=" << std::fixed << std::setprecision(3) << gausFit->GetParameter(1)).str().c_str(), "");
                        legPions->AddEntry((TObject*)0, (std::stringstream() << "#sigma=" << std::setprecision(3) << gausFit->GetParameter(2)).str().c_str(), "");
                        legPions->AddEntry((TObject*)0, (std::stringstream() << "#chi^{2}/NDF=" << std::setprecision(2) << chi2Ndf).str().c_str(), "");
                    }
                    legPions->Draw();
                }
            }
            canvasGaus->cd(2);
            if (histKaonsGaus) {
                histKaonsGaus->Draw("HIST");
                if (histKaonsGaus->GetFunction("gausFit")) {
                    histKaonsGaus->GetFunction("gausFit")->SetLineColor(kRed+2);
                    histKaonsGaus->GetFunction("gausFit")->SetLineStyle(1);
                    histKaonsGaus->GetFunction("gausFit")->SetLineWidth(2);
                    histKaonsGaus->GetFunction("gausFit")->Draw("SAME");
                }
                double yMaxKaons = histKaonsGaus->GetMaximum() * 1.1;
                std::vector<TLine*> linesKaons = drawMeanSigmaLines(histKaonsGaus, 0, yMaxKaons);
                std::vector<TObject*> legendObjects = {histKaonsGaus};
                std::vector<std::string> legendLabels = {"Kaons"};
                if (linesKaons.size() >= 3) {
                    legendObjects.push_back(linesKaons[1]);
                    legendLabels.push_back("Gaus #mu");
                    legendObjects.push_back(linesKaons[0]);
                    legendLabels.push_back("Gaus #mu #pm 3#sigma");
                }
                TLegend* legKaons = addLegend({histKaonsGaus}, {"Kaons"}, 0.7, 0.6, 0.9, 0.9, "");
                if (legKaons) {
                    legKaons->Clear();
                    for (size_t j = 0; j < legendObjects.size(); ++j) {
                        legKaons->AddEntry(legendObjects[j], legendLabels[j].c_str(), j == 0 ? "p" : "l");
                    }
                    if (TF1* gausFit = histKaonsGaus->GetFunction("gausFit")) {
                        double chi2Ndf = gausFit->GetNDF() > 0 ? gausFit->GetChisquare() / gausFit->GetNDF() : 0.0;
                        legKaons->AddEntry((TObject*)0, (std::stringstream() << "Entries: " << static_cast<int>(histKaonsGaus->GetEntries())).str().c_str(), "");
                        legKaons->AddEntry(gausFit, "Gauss:", "l");
                        legKaons->AddEntry((TObject*)0, (std::stringstream() << "Constant=" << std::fixed << std::setprecision(2) << gausFit->GetParameter(0)).str().c_str(), "");
                        legKaons->AddEntry((TObject*)0, (std::stringstream() << "#mu=" << std::fixed << std::setprecision(3) << gausFit->GetParameter(1)).str().c_str(), "");
                        legKaons->AddEntry((TObject*)0, (std::stringstream() << "#sigma=" << std::setprecision(3) << gausFit->GetParameter(2)).str().c_str(), "");
                        legKaons->AddEntry((TObject*)0, (std::stringstream() << "#chi^{2}/NDF=" << std::setprecision(2) << chi2Ndf).str().c_str(), "");
                    }
                    legKaons->Draw();
                }
            }

            // Save to PDF
            canvasGaus->Print("output/pdf/task1_chi2pid_gaus_only.pdf");

            // Save to ROOT file
            TFile* outRootGaus = new TFile(TString::Format("output/root/task1_gaus_bin%zu.root", i), "RECREATE");
            canvasGaus->Write();
            outRootGaus->Close();
            delete outRootGaus;
        }

        // Clean up
        delete canvas;
        delete canvasGaus;
        delete histPions;
        delete histKaons;
        delete histPionsGaus;
        delete histKaonsGaus;
    }

    // Close PDFs
    pdfCanvas->Print("output/pdf/task1_chi2pid.pdf]");
    pdfCanvasGaus->Print("output/pdf/task1_chi2pid_gaus_only.pdf]");
    delete pdfCanvas;
    delete pdfCanvasGaus;



        // Task 2: 1D beta histograms in p bins, before and after chi2pid cut
    // Open PDF for original Task 2 (using fitParams with CB or Gaussian)
    TCanvas* task2Canvas = new TCanvas("task2Canvas", "task2Canvas", 1600, 800);
    task2Canvas->SetCanvasSize(1600, 800);
    task2Canvas->Print("output/pdf/task2_beta.pdf[");

    // Configuration for cut type (true for CB, false for Gaussian)
    bool useCBcut = true;

    for (size_t i = 0; i < pBins.size() - 1; ++i) {
        // Define cuts for original Task 2
        TString pCut = TString::Format("p >= %f && p < %f && beta >= 0.96 && beta <= 1.03", pBins[i], pBins[i+1]);
        double meanPions = useCBcut ? fitParams["pions_cb_mean"][i] : fitParams["pions_gaus_mean"][i];
        double sigmaPions = useCBcut ? fitParams["pions_cb_sigma"][i] : fitParams["pions_gaus_sigma"][i];
        double meanKaons = useCBcut ? fitParams["kaons_cb_mean"][i] : fitParams["kaons_gaus_mean"][i];
        double sigmaKaons = useCBcut ? fitParams["kaons_cb_sigma"][i] : fitParams["kaons_gaus_sigma"][i];
        TString cutPions = TString::Format("%s && abs(chi2pid - %f) <= %f", pCut.Data(), meanPions, 3 * sigmaPions);
        TString cutKaons = TString::Format("%s && abs(chi2pid - %f) <= %f", pCut.Data(), meanKaons, 3 * sigmaKaons);

        // Create beta histograms (before cut) with fixed range, as points
        TH1F* betaPionsBefore = createHistogram(treePions, "beta", {},
                                                TString::Format("p: [%.2f-%.2f) (GeV/c)", pBins[i], pBins[i+1]),
                                                "beta", "Counts", kBlue+2, kBlue+2, 3001, "beta_pions_before",
                                                0.96, 1.03, "P", 0.8);
        TH1F* betaKaonsBefore = createHistogram(treeKaons, "beta", {},
                                                TString::Format("p: [%.2f-%.2f) (GeV/c)", pBins[i], pBins[i+1]),
                                                "beta", "Counts", kGreen+2, kGreen+2, 3002, "beta_kaons_before",
                                                0.96, 1.03, "P", 0.8);

        // Apply p and beta cut
        if (betaPionsBefore) {
            betaPionsBefore->Reset();
            treePions->Draw(TString::Format("beta>>%s", betaPionsBefore->GetName()), pCut, "goff");
            std::cout << "Beta pions before (bin " << i << ") entries: " << betaPionsBefore->GetEntries() << std::endl;
        }
        if (betaKaonsBefore) {
            betaKaonsBefore->Reset();
            treeKaons->Draw(TString::Format("beta>>%s", betaKaonsBefore->GetName()), pCut, "goff");
            std::cout << "Beta kaons before (bin " << i << ") entries: " << betaKaonsBefore->GetEntries() << std::endl;
        }

        // Create beta histograms (after cut) with fixed range, as points
        TH1F* betaPionsAfter = createHistogram(treePions, "beta", {},
                                               TString::Format("p: [%.2f-%.2f) (GeV/c)", pBins[i], pBins[i+1]),
                                               "beta", "Counts", kBlue+2, kBlue+2, 3001, "beta_pions_after",
                                               0.96, 1.03, "P", 0.8);
        TH1F* betaKaonsAfter = createHistogram(treeKaons, "beta", {},
                                               TString::Format("p:[%.2f-%.2f] (GeV/c)", pBins[i], pBins[i+1]),
                                               "beta", "Counts", kGreen+2, kGreen+2, 3002, "beta_kaons_after",
                                               0.96, 1.03, "P", 0.8);

        // Apply p, beta, and chi2pid cut
        if (betaPionsAfter) {
            betaPionsAfter->Reset();
            treePions->Draw(TString::Format("beta>>%s", betaPionsAfter->GetName()), cutPions, "goff");
            std::cout << "Beta pions after (bin " << i << ") entries: " << betaPionsAfter->GetEntries() << std::endl;
        }
        if (betaKaonsAfter) {
            betaKaonsAfter->Reset();
            treeKaons->Draw(TString::Format("beta>>%s", betaKaonsAfter->GetName()), cutKaons, "goff");
            std::cout << "Beta kaons after (bin " << i << ") entries: " << betaKaonsAfter->GetEntries() << std::endl;
        }

        // Fit before-cut histograms with fitHistogram (Gaussian only)
        if (betaPionsBefore && betaPionsBefore->GetEntries() > 50) {
            fitHistogram(betaPionsBefore, "fit_pions_before", 0.96, 1.03);
        } else {
            std::cerr << "Skipping fit for beta pions before (bin " << i << ") due to insufficient entries." << std::endl;
        }
        if (betaKaonsBefore && betaKaonsBefore->GetEntries() > 50) {
            fitHistogram(betaKaonsBefore, "fit_kaons_before", 0.96, 1.03);
        } else {
            std::cerr << "Skipping fit for beta kaons before (bin " << i << ") due to insufficient entries." << std::endl;
        }

        // Fit after-cut histograms with fitHistogram (Gaussian only)
        std::map<std::string, double> fitPionsAfter, fitKaonsAfter;
        if (betaPionsAfter && betaPionsAfter->GetEntries() > 50) {
            fitPionsAfter = fitHistogram(betaPionsAfter, "fit_pions_after", 0.96, 1.03);
        } else {
            std::cerr << "Skipping fit for beta pions after (bin " << i << ") due to insufficient entries." << std::endl;
        }
        if (betaKaonsAfter && betaKaonsAfter->GetEntries() > 50) {
            fitKaonsAfter = fitHistogram(betaKaonsAfter, "fit_kaons_after", 0.96, 1.03);
        } else {
            std::cerr << "Skipping fit for beta kaons after (bin " << i << ") due to insufficient entries." << std::endl;
        }

        // Calculate contamination using the after-cut beta histograms
        if (fitPionsAfter.count("gaus_constant") && fitPionsAfter.count("gaus_mean") && fitPionsAfter.count("gaus_sigma") &&
            fitKaonsAfter.count("gaus_constant") && fitKaonsAfter.count("gaus_mean") && fitKaonsAfter.count("gaus_sigma") &&
            fitPionsAfter["gaus_sigma"] > 0 && fitKaonsAfter["gaus_sigma"] > 0 &&
            fitPionsAfter["gaus_constant"] > 0 && fitKaonsAfter["gaus_constant"] > 0) {
            std::cout << "Fit parameters for p bin [" << pBins[i] << "-" << pBins[i+1] << "]:" << std::endl;
            std::cout << "Pions - Constant: " << std::fixed << std::setprecision(6) << fitPionsAfter["gaus_constant"]
                      << ", Mean: " << fitPionsAfter["gaus_mean"]
                      << ", Sigma: " << fitPionsAfter["gaus_sigma"] << std::endl;
            std::cout << "Kaons - Constant: " << std::fixed << std::setprecision(6) << fitKaonsAfter["gaus_constant"]
                      << ", Mean: " << fitKaonsAfter["gaus_mean"]
                      << ", Sigma: " << fitKaonsAfter["gaus_sigma"] << std::endl;

            double c1, c2;
            double contamination = calculateOverlapContamination(fitPionsAfter, fitKaonsAfter, c1, c2);
            csvFile << pBins[i] << "-" << pBins[i+1] << ",";
            if (contamination >= 0) {
                csvFile << std::fixed << std::setprecision(6) << c1 << "," << c2 << ","
                        << std::fixed << std::setprecision(4) << contamination << "\n";
                std::cout << "Contamination for p bin [" << pBins[i] << "-" << pBins[i+1] << "] GeV/c (after chi2pid cut): "
                          << std::fixed << std::setprecision(4) << contamination << "%" << std::endl;
            } else {
                csvFile << "N/A,N/A,N/A\n";
                std::cerr << "Contamination calculation failed for p bin [" << pBins[i] << "-" << pBins[i+1] << "]" << std::endl;
            }
        } else {
            std::cerr << "Skipping contamination calculation for p bin [" << pBins[i] << "-" << pBins[i+1]
                      << "] due to missing or invalid fit parameters after cut." << std::endl;
            csvFile << pBins[i] << "-" << pBins[i+1] << ",N/A,N/A,N/A\n";
        }

        // Create canvas (1x2: before cut left, after cut right) for original Task 2
        std::vector<TH1F*> histograms = {betaPionsBefore, betaKaonsBefore, betaPionsAfter, betaKaonsAfter};
        TCanvas* canvas = createCanvas(histograms, 1, 2, {1, 1, 2, 2},
                                      TString::Format("canvas_task2_bin%zu", i), 1600, 800);

        // Add legends, draw fits, and set logarithmic scale
        if (canvas) {
            canvas->cd(1);
            gPad->SetLogy();
            if (betaPionsBefore) {betaPionsBefore->SetMarkerSize(1.5); betaPionsBefore->Draw("P");}
            if (betaKaonsBefore) {betaKaonsBefore->SetMarkerSize(1.5); betaKaonsBefore->Draw("P SAME");}
            if (betaPionsBefore && betaPionsBefore->GetFunction("gausFit")) {
                betaPionsBefore->GetFunction("gausFit")->SetLineColor(kBlue);
                betaPionsBefore->GetFunction("gausFit")->SetLineStyle(kDashed);
                betaPionsBefore->GetFunction("gausFit")->Draw("SAME");
            }
            if (betaKaonsBefore && betaKaonsBefore->GetFunction("gausFit")) {
                betaKaonsBefore->GetFunction("gausFit")->SetLineColor(kGreen);
                betaKaonsBefore->GetFunction("gausFit")->SetLineStyle(kDashed);
                betaKaonsBefore->GetFunction("gausFit")->Draw("SAME");
            }
            TLegend* legBefore = addLegend({betaPionsBefore, betaKaonsBefore}, {"EB Pions", "EB Kaons"},
                                          0.65, 0.6, 0.9, 0.9, "");
            if (legBefore) legBefore->Draw();

            canvas->cd(2);
            gPad->SetLogy();
            if (betaPionsAfter) {betaPionsAfter->SetMarkerSize(1.5); betaPionsAfter->Draw("P");}
            if (betaKaonsAfter) {betaKaonsAfter->SetMarkerSize(1.5); betaKaonsAfter->Draw("P SAME");}
            if (betaPionsAfter && betaPionsAfter->GetFunction("gausFit")) {
                betaPionsAfter->GetFunction("gausFit")->SetLineColor(kBlue);
                betaPionsAfter->GetFunction("gausFit")->SetLineStyle(kDashed);
                betaPionsAfter->GetFunction("gausFit")->Draw("SAME");
            }
            if (betaKaonsAfter && betaKaonsAfter->GetFunction("gausFit")) {
                betaKaonsAfter->GetFunction("gausFit")->SetLineColor(kGreen);
                betaKaonsAfter->GetFunction("gausFit")->SetLineStyle(kDashed);
                betaKaonsAfter->GetFunction("gausFit")->Draw("SAME");
            }
            TLegend* legAfter = addLegend({betaPionsAfter, betaKaonsAfter}, {"EB Pions After Cut", "EB Kaons After Cut"},
                                         0.65, 0.6, 0.9, 0.9, "");
            if (legAfter) legAfter->Draw();

            // Save to PDF for original Task 2
            canvas->Print("output/pdf/task2_beta.pdf");

            // Save to ROOT file
            TFile* outRoot = new TFile(TString::Format("output/root/task2_bin%zu.root", i), "RECREATE");
            canvas->Write();
            outRoot->Close();
            delete outRoot;
        }

        // Clean up for original Task 2 histograms
        delete betaPionsBefore;
        delete betaKaonsBefore;
        delete betaPionsAfter;
        delete betaKaonsAfter;
    }

    // Close PDF for original Task 2
    task2Canvas->Print("output/pdf/task2_beta.pdf]");
    delete task2Canvas;

    std::ofstream csvFileGaus("output/csv/task2_beta_gaus_contamination.csv");
    csvFileGaus << "p-bin,c1,c2,contamination\n";

    // Open PDF for new Task 2 (using fitParamsGaus from Task 1 Gaussian-only)
// Open PDF for new Task 2 (using fitParamsGaus from Task 1 Gaussian-only)
TCanvas* task2GausCanvas = new TCanvas("pdfCanvasGaus2", "pdfCanvasGaus2", 1600, 800);
task2GausCanvas->SetCanvasSize(1600, 800);
task2GausCanvas->Print("output/pdf/task2_beta_gaus_only.pdf[");

for (size_t i = 0; i < pBins.size() - 1; ++i) {
    // Define cuts for new Task 2 using 3-sigma from Task 1 Gaussian-only fits
    TString pCut = TString::Format("p >= %f && p < %f && beta >= 0.96 && beta <= 1.03", pBins[i], pBins[i+1]);
    double meanPionsGaus = fitParamsGaus["pions_gaus_mean"][i];
    double sigmaPionsGaus = fitParamsGaus["pions_gaus_sigma"][i];
    double meanKaonsGaus = fitParamsGaus["kaons_gaus_mean"][i];
    double sigmaKaonsGaus = fitParamsGaus["kaons_gaus_sigma"][i];
    TString cutPionsGaus = TString::Format("%s && abs(chi2pid - %f) <= %f", pCut.Data(), meanPionsGaus, 3 * sigmaPionsGaus);
    TString cutKaonsGaus = TString::Format("%s && abs(chi2pid - %f) <= %f", pCut.Data(), meanKaonsGaus, 3 * sigmaKaonsGaus);

    // Create beta histograms (before cut) with fixed range, as points
    TH1F* betaPionsBeforeGaus = createHistogram(treePions, "beta", {},
                                                TString::Format("p: [%.2f-%.2f) (GeV/c)", pBins[i], pBins[i+1]),
                                                "beta", "Counts", kBlue+2, kBlue+2, 3001, "beta_pions_before_gaus",
                                                0.96, 1.03, "P", 0.8);
    TH1F* betaKaonsBeforeGaus = createHistogram(treeKaons, "beta", {},
                                                TString::Format("p: [%.2f-%.2f) (GeV/c)", pBins[i], pBins[i+1]),
                                                "beta", "Counts", kGreen+2, kGreen+2, 3002, "beta_kaons_before_gaus",
                                                0.96, 1.03, "P", 0.8);

    // Apply p and beta cut
    if (betaPionsBeforeGaus) {
        betaPionsBeforeGaus->Reset();
        treePions->Draw(TString::Format("beta>>%s", betaPionsBeforeGaus->GetName()), pCut, "goff");
        std::cout << "Beta pions before (bin " << i << ") entries (Gaus): " << betaPionsBeforeGaus->GetEntries() << std::endl;
    }
    if (betaKaonsBeforeGaus) {
        betaKaonsBeforeGaus->Reset();
        treeKaons->Draw(TString::Format("beta>>%s", betaKaonsBeforeGaus->GetName()), pCut, "goff");
        std::cout << "Beta kaons before (bin " << i << ") entries (Gaus): " << betaKaonsBeforeGaus->GetEntries() << std::endl;
    }

    // Create beta histograms (after cut) with fixed range, as points
    TH1F* betaPionsAfterGaus = createHistogram(treePions, "beta", {},
                                               TString::Format("p: [%.2f-%.2f) (GeV/c)", pBins[i], pBins[i+1]),
                                               "beta", "Counts", kBlue+2, kBlue+2, 3001, "beta_pions_after_gaus",
                                               0.96, 1.03, "P", 0.8);
    TH1F* betaKaonsAfterGaus = createHistogram(treeKaons, "beta", {},
                                               TString::Format("p:[%.2f-%.2f] (GeV/c)", pBins[i], pBins[i+1]),
                                               "beta", "Counts", kGreen+2, kGreen+2, 3002, "beta_kaons_after_gaus",
                                               0.96, 1.03, "P", 0.8);

    // Apply p, beta, and chi2pid cut using 3-sigma from Gaussian-only fits
    if (betaPionsAfterGaus) {
        betaPionsAfterGaus->Reset();
        treePions->Draw(TString::Format("beta>>%s", betaPionsAfterGaus->GetName()), cutPionsGaus, "goff");
        std::cout << "Beta pions after (bin " << i << ") entries (Gaus): " << betaPionsAfterGaus->GetEntries() << std::endl;
    }
    if (betaKaonsAfterGaus) {
        betaKaonsAfterGaus->Reset();
        treeKaons->Draw(TString::Format("beta>>%s", betaKaonsAfterGaus->GetName()), cutKaonsGaus, "goff");
        std::cout << "Beta kaons after (bin " << i << ") entries (Gaus): " << betaKaonsAfterGaus->GetEntries() << std::endl;
    }

    // Fit before-cut histograms with fitHistogram (Gaussian only) - no threshold
    if (betaPionsBeforeGaus) {
        fitHistogram(betaPionsBeforeGaus, "fit_pions_before_gaus", 0.96, 1.03);
    }
    if (betaKaonsBeforeGaus) {
        fitHistogram(betaKaonsBeforeGaus, "fit_kaons_before_gaus", 0.96, 1.03);
    }

    // Fit after-cut histograms with fitHistogram (Gaussian only) - no threshold
    std::map<std::string, double> fitPionsAfterGaus, fitKaonsAfterGaus;
    if (betaPionsAfterGaus) {
        if (i == 7) {
            fitPionsAfterGaus = fitHistogram(betaPionsAfterGaus, "fit_pions_after_gaus", 0.98, 1.02);
        } else {
            fitPionsAfterGaus = fitHistogram(betaPionsAfterGaus, "fit_pions_after_gaus", 0.96, 1.03);
        }
    }
    if (betaKaonsAfterGaus) {
        if (i == 7) {
            fitKaonsAfterGaus = fitHistogram(betaKaonsAfterGaus, "fit_kaons_after_gaus", 0.972, 1.01);
        } else {
            fitKaonsAfterGaus = fitHistogram(betaKaonsAfterGaus, "fit_kaons_after_gaus", 0.96, 1.03);
        }
    }

    // For i == 7, save the after-cut plot as PNG and calculate contamination for testing
    if (i == 7) {
        // Create canvas for after-cut histograms only
        TCanvas* testCanvas = new TCanvas(TString::Format("test_canvas_gaus_bin%zu", i), "Test Canvas", 800, 600);
        testCanvas->cd();
        gPad->SetLogy();
        if (betaPionsAfterGaus) {
            betaPionsAfterGaus->SetMarkerSize(1.5);
            betaPionsAfterGaus->Draw("P");
        }
        if (betaKaonsAfterGaus) {
            betaKaonsAfterGaus->SetMarkerSize(1.5);
            betaKaonsAfterGaus->Draw("P SAME");
        }
        if (betaPionsAfterGaus && betaPionsAfterGaus->GetFunction("gausFit")) {
            betaPionsAfterGaus->GetFunction("gausFit")->SetLineColor(kBlue);
            betaPionsAfterGaus->GetFunction("gausFit")->SetLineStyle(kDashed);
            betaPionsAfterGaus->GetFunction("gausFit")->Draw("SAME");
        }
        if (betaKaonsAfterGaus && betaKaonsAfterGaus->GetFunction("gausFit")) {
            betaKaonsAfterGaus->GetFunction("gausFit")->SetLineColor(kGreen);
            betaKaonsAfterGaus->GetFunction("gausFit")->SetLineStyle(kDashed);
            betaKaonsAfterGaus->GetFunction("gausFit")->Draw("SAME");
        }
        TLegend* legTest = addLegend({betaPionsAfterGaus, betaKaonsAfterGaus}, 
                                    {"EB Pions After 3-sigma Cut", "EB Kaons After 3-sigma Cut"},
                                    0.65, 0.6, 0.9, 0.9, "");
        if (legTest) legTest->Draw();

        // Save the canvas as PNG
        testCanvas->SaveAs(TString::Format("output/png/task2_gaus_bin%zu_after.png", i));

        // Calculate contamination for i == 7 (test)
        if (fitPionsAfterGaus.count("gaus_constant") && fitPionsAfterGaus.count("gaus_mean") && 
            fitPionsAfterGaus.count("gaus_sigma") && fitKaonsAfterGaus.count("gaus_constant") && 
            fitKaonsAfterGaus.count("gaus_mean") && fitKaonsAfterGaus.count("gaus_sigma") &&
            fitPionsAfterGaus["gaus_sigma"] > 0 && fitKaonsAfterGaus["gaus_sigma"] > 0 &&
            fitPionsAfterGaus["gaus_constant"] > 0 && fitKaonsAfterGaus["gaus_constant"] > 0) {
            std::cout << "Test fit parameters for p bin [" << pBins[i] << "-" << pBins[i+1] << "] (Gaus, i=7):" << std::endl;
            std::cout << "Pions - Constant: " << std::fixed << std::setprecision(6) << fitPionsAfterGaus["gaus_constant"]
                      << ", Mean: " << fitPionsAfterGaus["gaus_mean"]
                      << ", Sigma: " << fitPionsAfterGaus["gaus_sigma"] << std::endl;
            std::cout << "Kaons - Constant: " << std::fixed << std::setprecision(6) << fitKaonsAfterGaus["gaus_constant"]
                      << ", Mean: " << fitKaonsAfterGaus["gaus_mean"]
                      << ", Sigma: " << fitKaonsAfterGaus["gaus_sigma"] << std::endl;

            double c1, c2;
            double contamination = calculateOverlapContamination(fitPionsAfterGaus, fitKaonsAfterGaus, c1, c2);
            if (contamination >= 0) {
                std::cout << "Test contamination for p bin [" << pBins[i] << "-" << pBins[i+1] << "] GeV/c (i=7, after 3-sigma Gaus cut): "
                          << std::fixed << std::setprecision(4) << contamination << "%" << std::endl;
            } else {
                std::cerr << "Test contamination calculation failed for p bin [" << pBins[i] << "-" << pBins[i+1] << "] (Gaus, i=7)" << std::endl;
            }
        } else {
            std::cerr << "Skipping test contamination calculation for p bin [" << pBins[i] << "-" << pBins[i+1] 
                      << "] (i=7) due to missing or invalid fit parameters." << std::endl;
        }

        delete testCanvas;
    }

    // Calculate contamination using the after-cut beta histograms for Gaussian-only case
    if (fitPionsAfterGaus.count("gaus_constant") && fitPionsAfterGaus.count("gaus_mean") && fitPionsAfterGaus.count("gaus_sigma") &&
        fitKaonsAfterGaus.count("gaus_constant") && fitKaonsAfterGaus.count("gaus_mean") && fitKaonsAfterGaus.count("gaus_sigma") &&
        fitPionsAfterGaus["gaus_sigma"] > 0 && fitKaonsAfterGaus["gaus_sigma"] > 0 &&
        fitPionsAfterGaus["gaus_constant"] > 0 && fitKaonsAfterGaus["gaus_constant"] > 0) {
        std::cout << "Fit parameters for p bin [" << pBins[i] << "-" << pBins[i+1] << "] (Gaus):" << std::endl;
        std::cout << "Pions - Constant: " << std::fixed << std::setprecision(6) << fitPionsAfterGaus["gaus_constant"]
                  << ", Mean: " << fitPionsAfterGaus["gaus_mean"]
                  << ", Sigma: " << fitPionsAfterGaus["gaus_sigma"] << std::endl;
        std::cout << "Kaons - Constant: " << std::fixed << std::setprecision(6) << fitKaonsAfterGaus["gaus_constant"]
                  << ", Mean: " << fitKaonsAfterGaus["gaus_mean"]
                  << ", Sigma: " << fitKaonsAfterGaus["gaus_sigma"] << std::endl;

        double c1, c2;
        double contamination = calculateOverlapContamination(fitPionsAfterGaus, fitKaonsAfterGaus, c1, c2);
        csvFileGaus << pBins[i] << "-" << pBins[i+1] << ",";
        if (contamination >= 0) {
            csvFileGaus << std::fixed << std::setprecision(6) << c1 << "," << c2 << ","
                        << std::fixed << std::setprecision(4) << contamination << "\n";
            std::cout << "Contamination for p bin [" << pBins[i] << "-" << pBins[i+1] << "] GeV/c (after 3-sigma Gaus cut): "
                      << std::fixed << std::setprecision(4) << contamination << "%" << std::endl;
        } else {
            csvFileGaus << "N/A,N/A,N/A\n";
            std::cerr << "Contamination calculation failed for p bin [" << pBins[i] << "-" << pBins[i+1] << "] (Gaus)" << std::endl;
        }
    } else {
        std::cerr << "Skipping contamination calculation for p bin [" << pBins[i] << "-" << pBins[i+1] << "] due to missing or invalid fit parameters after 3-sigma Gaus cut." << std::endl;
        csvFileGaus << pBins[i] << "-" << pBins[i+1] << ",N/A,N/A,N/A\n";
    }

    // Create canvas (1x2: before cut left, after cut right) for new Task 2
    std::vector<TH1F*> histogramsGaus = {betaPionsBeforeGaus, betaKaonsBeforeGaus, betaPionsAfterGaus, betaKaonsAfterGaus};
    TCanvas* canvasGaus = createCanvas(histogramsGaus, 1, 2, {1, 1, 2, 2},
                                      TString::Format("canvas_task2_gaus_bin%zu", i), 1600, 800);

    // Add legends, draw fits, and set logarithmic scale
    if (canvasGaus) {
        canvasGaus->cd(1);
        gPad->SetLogy();
        if (betaPionsBeforeGaus) {betaPionsBeforeGaus->SetMarkerSize(1.5); betaPionsBeforeGaus->Draw("P");}
        if (betaKaonsBeforeGaus) {betaKaonsBeforeGaus->SetMarkerSize(1.5); betaKaonsBeforeGaus->Draw("P SAME");}
        if (betaPionsBeforeGaus && betaPionsBeforeGaus->GetFunction("gausFit")) {
            betaPionsBeforeGaus->GetFunction("gausFit")->SetLineColor(kBlue);
            betaPionsBeforeGaus->GetFunction("gausFit")->SetLineStyle(kDashed);
            betaPionsBeforeGaus->GetFunction("gausFit")->Draw("SAME");
        }
        if (betaKaonsBeforeGaus && betaKaonsBeforeGaus->GetFunction("gausFit")) {
            betaKaonsBeforeGaus->GetFunction("gausFit")->SetLineColor(kGreen);
            betaKaonsBeforeGaus->GetFunction("gausFit")->SetLineStyle(kDashed);
            betaKaonsBeforeGaus->GetFunction("gausFit")->Draw("SAME");
        }
        TLegend* legBeforeGaus = addLegend({betaPionsBeforeGaus, betaKaonsBeforeGaus}, {"EB Pions", "EB Kaons"},
                                          0.65, 0.6, 0.9, 0.9, "");
        if (legBeforeGaus) legBeforeGaus->Draw();

        canvasGaus->cd(2);
        gPad->SetLogy();
        if (betaPionsAfterGaus) {betaPionsAfterGaus->SetMarkerSize(1.5); betaPionsAfterGaus->Draw("P");}
        if (betaKaonsAfterGaus) {betaKaonsAfterGaus->SetMarkerSize(1.5); betaKaonsAfterGaus->Draw("P SAME");}
        if (betaPionsAfterGaus && betaPionsAfterGaus->GetFunction("gausFit")) {
            betaPionsAfterGaus->GetFunction("gausFit")->SetLineColor(kBlue);
            betaPionsAfterGaus->GetFunction("gausFit")->SetLineStyle(kDashed);
            betaPionsAfterGaus->GetFunction("gausFit")->Draw("SAME");
        }
        if (betaKaonsAfterGaus && betaKaonsAfterGaus->GetFunction("gausFit")) {
            betaKaonsAfterGaus->GetFunction("gausFit")->SetLineColor(kGreen);
            betaKaonsAfterGaus->GetFunction("gausFit")->SetLineStyle(kDashed);
            betaKaonsAfterGaus->GetFunction("gausFit")->Draw("SAME");
        }
        TLegend* legAfterGaus = addLegend({betaPionsAfterGaus, betaKaonsAfterGaus}, {"EB Pions After 3-sigma Cut", "EB Kaons After 3-sigma Cut"},
                                         0.65, 0.6, 0.9, 0.9, "");
        if (legAfterGaus) legAfterGaus->Draw();

        // Draw and label the chosen intersection point and endpoints on the after-cut plot
        if (fitPionsAfterGaus.count("gaus_mean") && fitPionsAfterGaus.count("gaus_sigma") &&
            fitKaonsAfterGaus.count("gaus_mean") && fitKaonsAfterGaus.count("gaus_sigma")) {
            double c1, c2;
            double contamination = calculateOverlapContamination(fitPionsAfterGaus, fitKaonsAfterGaus, c1, c2);

            // Determine the chosen intersection point
            double mu1 = fitPionsAfterGaus["gaus_mean"];
            double mu2 = fitKaonsAfterGaus["gaus_mean"];
            if (mu1 < mu2) std::swap(mu1, mu2);
            if (c1 > c2) std::swap(c1, c2);
            double c = -1.0;
            if (c1 >= mu2 && c1 <= mu1) {
                c = c1;
            } else if (c2 >= mu2 && c2 <= mu1) {
                c = c2;
            }

            // Calculate endpoints (where Gaussian drops to 10^-4 of peak height or < 1, capped between 0.96 and 1.0)
            double C1 = fitPionsAfterGaus["gaus_constant"];
            double sigma1 = fitPionsAfterGaus["gaus_sigma"];
            double C2 = fitKaonsAfterGaus["gaus_constant"];
            double sigma2 = fitKaonsAfterGaus["gaus_sigma"];

            TF1* pionFit = new TF1("pionFit", "gaus", 0.96, 1.0);
            pionFit->SetParameters(C1, mu1, sigma1);
            TF1* kaonFit = new TF1("kaonFit", "gaus", 0.96, 1.0);
            kaonFit->SetParameters(C2, mu2, sigma2);

            double pionLeftEndpoint = mu1;
            double threshold = std::min(1.0, C1 * 1e-4); // Use 10^-4 of peak or 1, whichever is smaller
            while (pionFit->Eval(pionLeftEndpoint) > threshold && pionLeftEndpoint > 0.96) {
                pionLeftEndpoint -= 0.001 * sigma1;
            }
            pionLeftEndpoint = std::max(0.96, pionLeftEndpoint);

            double pionRightEndpoint = mu1;
            while (pionFit->Eval(pionRightEndpoint) > threshold && pionRightEndpoint < 1.0) {
                pionRightEndpoint += 0.001 * sigma1;
            }
            pionRightEndpoint = std::min(1.0, pionRightEndpoint);

            double kaonRightEndpoint = mu2;
            threshold = std::min(1.0, C2 * 1e-4); // Use 10^-4 of peak or 1, whichever is smaller
            while (kaonFit->Eval(kaonRightEndpoint) > threshold && kaonRightEndpoint < 1.0) {
                kaonRightEndpoint += 0.001 * sigma2;
            }
            kaonRightEndpoint = std::min(1.0, kaonRightEndpoint);

            delete pionFit;
            delete kaonFit;

            // Draw the chosen intersection point if valid
            if (c != -1.0) {
                TLine* lineC = new TLine(c, gPad->GetUymin(), c, gPad->GetUymax());
                lineC->SetLineColor(kRed);
                lineC->SetLineStyle(kDashed);
                lineC->Draw();
                TLatex* labelC = new TLatex(c, gPad->GetUymax() * 0.9, Form("c = %.6f", c));
                labelC->SetTextSize(0.03);
                labelC->SetTextColor(kRed);
                labelC->Draw();
            }
        }

        // Save to PDF for new Task 2
        canvasGaus->Print("output/pdf/task2_beta_gaus_only.pdf");

        // Save to ROOT file
        TFile* outRootGaus = new TFile(TString::Format("output/root/task2_gaus_bin%zu.root", i), "RECREATE");
        canvasGaus->Write();
        outRootGaus->Close();
        delete outRootGaus;
    }

    // Clean up for new Task 2 histograms
    delete betaPionsBeforeGaus;
    delete betaKaonsBeforeGaus;
    delete betaPionsAfterGaus;
    delete betaKaonsAfterGaus;
}

// Close PDF for new Task 2
task2GausCanvas->Print("output/pdf/task2_beta_gaus_only.pdf]");
delete task2GausCanvas;
// Close CSV file
csvFileGaus.close();

// Clean up
file->Close();
delete file;

return 0;
} 


