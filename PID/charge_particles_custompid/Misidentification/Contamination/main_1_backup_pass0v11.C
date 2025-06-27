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

// Helper function to create output_pass0v11 directories
void createoutput_pass0v11Dirs() {
    gSystem->mkdir("output_pass0v11", kTRUE);
    gSystem->mkdir("output_pass0v11/pdf", kTRUE);
    gSystem->mkdir("output_pass0v11/root", kTRUE);
}

int main() {
    // Load ROOT file
    TFile* file = new TFile("../pkptreeCxC_5.root");
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

    // Create output_pass0v11 directories
    createoutput_pass0v11Dirs();

    // Open CSV file for contamination results
    std::ofstream csvFile("output_pass0v11/contamination.csv");
    csvFile << "Momentum Bin (GeV/c),c1,c2,Contamination (%)\n";

    // Task 1: 1D chi2pid histograms in p bins with Gaussian and Crystal Ball fits
    // Define momentum bins
    std::vector<double> pBins = createBins(1, 10, 30); // 10 linear bins from 1 to 4 GeV/c
    if (pBins.empty()) {
        std::cerr << "Error: Failed to create p bins." << std::endl;
        file->Close();
        delete file;
        return 1;
    }

    // Open PDF for Task 1
    TCanvas* pdfCanvas = new TCanvas("pdfCanvas", "pdfCanvas", 1600, 800);
    pdfCanvas->SetCanvasSize(1600, 800);
    pdfCanvas->Print("output_pass0v11/pdf/task1_chi2pid.pdf[");

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

    // Process each p bin for Task 1
    for (size_t i = 0; i < pBins.size() - 1; ++i) {
        // Define p cut
        TString pCut = TString::Format("p >= %f && p < %f", pBins[i], pBins[i+1]);

        // Create chi2pid histograms with fixed range
        TH1F* histPions = createHistogram(treePions, "chi2pid", {},
                                          TString::Format("p: [%.2f-%.2f) GeV/c", pBins[i], pBins[i+1]),
                                          "chi2pid", "Counts", kBlue, kBlue, 3001, "hist_pions",
                                          -10.0, 10.0);
        TH1F* histKaons = createHistogram(treeKaons, "chi2pid", {},
                                          TString::Format("p: [%.2f-%.2f) GeV/c", pBins[i], pBins[i+1]),
                                          "chi2pid", "Counts", kGreen, kGreen, 3002, "hist_kaons",
                                          -10.0, 10.0);
        if (i > 10) {
            //histPions->SetXaxisLimits(-4, 4);
            histKaons->GetXaxis()->SetRangeUser(-5, 4);
        }

        // Apply p cut using TTree::Draw
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

        // Fit histograms with both Gaussian and Crystal Ball
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

        // Create canvas (1x2: pions left, kaons right)
        std::vector<TH1F*> histograms = {histPions, histKaons};
        TCanvas* canvas = createCanvas(histograms, 1, 2, {1, 2},
                                      TString::Format("canvas_task1_bin%zu", i), 1600, 800);

        // Add legend, draw fits, and mean ± 3σ lines
        if (canvas) {
            canvas->cd(1);
            if (histPions) {
                histPions->Draw("HIST");
                if (histPions->GetFunction("gausFit")) {
                    histPions->GetFunction("gausFit")->SetLineColor(kBlue+2); // Darker blue for Gaussian
                    histPions->GetFunction("gausFit")->SetLineStyle(1); // Solid line
                    histPions->GetFunction("gausFit")->SetLineWidth(2); // Thicker line for visibility
                    histPions->GetFunction("gausFit")->Draw("SAME");
                }
                if (histPions->GetFunction("cbFit")) {
                    histPions->GetFunction("cbFit")->SetLineColor(kRed+2); // Darker red for Crystal Ball
                    histPions->GetFunction("cbFit")->SetLineStyle(1); // Solid line
                    histPions->GetFunction("cbFit")->SetLineWidth(2); // Thicker line for visibility
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
                    histKaons->GetFunction("gausFit")->SetLineColor(kBlue+2); // Darker blue for Gaussian
                    histKaons->GetFunction("gausFit")->SetLineStyle(1); // Solid line
                    histKaons->GetFunction("gausFit")->SetLineWidth(2); // Thicker line for visibility
                    histKaons->GetFunction("gausFit")->Draw("SAME");
                }
                if (histKaons->GetFunction("cbFit")) {
                    histKaons->GetFunction("cbFit")->SetLineColor(kRed+2); // Darker red for Crystal Ball
                    histKaons->GetFunction("cbFit")->SetLineStyle(1); // Solid line
                    histKaons->GetFunction("cbFit")->SetLineWidth(2); // Thicker line for visibility
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
            canvas->Print("output_pass0v11/pdf/task1_chi2pid.pdf");

            // Save to ROOT file
            TFile* outRoot = new TFile(TString::Format("output_pass0v11/root/task1_bin%zu.root", i), "RECREATE");
            canvas->Write();
            outRoot->Close();
            delete outRoot;
        }

        // Clean up
        delete canvas;
        delete histPions;
        delete histKaons;
    }

    // Close PDF
    pdfCanvas->Print("output_pass0v11/pdf/task1_chi2pid.pdf]");
    delete pdfCanvas;

    // Task 2: 1D beta histograms in p bins, before and after chi2pid cut
    // Open PDF for Task 2
    pdfCanvas = new TCanvas("pdfCanvas2", "pdfCanvas2", 1600, 800);
    pdfCanvas->SetCanvasSize(1600, 800);
    pdfCanvas->Print("output_pass0v11/pdf/task2_beta.pdf[");

    // Configuration for cut type (true for CB, false for Gaussian)
    bool useCBcut = true;

    for (size_t i = 0; i < pBins.size() - 1; ++i) {
        // Define cuts
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
            // Log fit parameters for debugging
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

        // Create canvas (1x2: before cut left, after cut right)
        std::vector<TH1F*> histograms = {betaPionsBefore, betaKaonsBefore, betaPionsAfter, betaKaonsAfter};
        TCanvas* canvas = createCanvas(histograms, 1, 2, {1, 1, 2, 2}, // Overlay in each pad
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

            // Save to PDF
            canvas->Print("output_pass0v11/pdf/task2_beta.pdf");

            // Save to ROOT file
            TFile* outRoot = new TFile(TString::Format("output_pass0v11/root/task2_bin%zu.root", i), "RECREATE");
            canvas->Write();
            outRoot->Close();
            delete outRoot;
        }

        // Clean up
        delete canvas;
        delete betaPionsBefore;
        delete betaKaonsBefore;
        delete betaPionsAfter;
        delete betaKaonsAfter;
    }

    // Close PDF
    pdfCanvas->Print("output_pass0v11/pdf/task2_beta.pdf]");
    delete pdfCanvas;

    // Close CSV file
    csvFile.close();

    // Clean up
    file->Close();
    delete file;

    return 0;
}