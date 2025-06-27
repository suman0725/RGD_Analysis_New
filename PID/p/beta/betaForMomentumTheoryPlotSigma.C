#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <ROOT/TTreeProcessorMT.hxx>
#include <TStopwatch.h>
#include <cmath>
#include <ROOT/TThreadedObject.hxx>
#include <TLegend.h>
#include <sstream>

void betaForMomentumTheoryPlotSigma() {
    TStopwatch timer;
    timer.Start();

    // Open the ROOT file
    TFile* file = TFile::Open("../charged_particles.root", "READ");
    if (!file || !file->IsOpen()) {
        std::cerr << "Error: Unable to open the ROOT file!" << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("charged_particle");
    if (!tree) {
        std::cerr << "Error: Unable to access the TTree!" << std::endl;
        return;
    }

    ROOT::TThreadedObject<TH2F> threaded_histo(
        "threaded_histo",
        "Beta vs Momentum;Momentum (GeV);Beta",
        200, 0, 8, 200, 0.44, 1.02);

    ROOT::TTreeProcessorMT processor(*tree);
    processor.Process([&](TTreeReader& reader) {
        TTreeReaderValue<int> charge(reader, "charge");
        TTreeReaderValue<short> status(reader, "status");
        TTreeReaderValue<float> px(reader, "px");
        TTreeReaderValue<float> py(reader, "py");
        TTreeReaderValue<float> pz(reader, "pz");
        TTreeReaderValue<float> beta(reader, "beta");

        auto local_histo = threaded_histo.Get();

        while (reader.Next()) {
            if ((abs(*status) / 2000 == 1) && (*charge > 0)) {
                float momentum = sqrt((*px) * (*px) + (*py) * (*py) + (*pz) * (*pz));
                local_histo->Fill(momentum, *beta);
            }
        }
    });

    TH2F* h_beta_vs_momentum = new TH2F(
        "h_beta_vs_momentum",
        "Beta vs Momentum;Momentum (GeV);Beta",
        200, 0, 8, 200, 0.44, 1.02);

    threaded_histo.Merge([&](std::shared_ptr<TH2F> main, std::vector<std::shared_ptr<TH2F>>& histos) {
        for (auto& histo : histos) {
            main->Add(histo.get());
        }
        *h_beta_vs_momentum = *main;
    });

    TCanvas* canvas = new TCanvas("canvas", "Beta vs. Momentum", 800, 600);
    canvas->SetLogz();
    h_beta_vs_momentum->SetMarkerStyle(20);
    h_beta_vs_momentum->Draw("COLZ");
    h_beta_vs_momentum->SetStats(0);

    // Read the fit results for pions from the file
    std::ifstream fitFile("fit_results_two_gauss_added_0.2_0.6.txt");
    if (!fitFile.is_open()) {
        std::cerr << "Error: Could not open fit_results_two_gauss.txt!" << std::endl;
        return;
    }

    // Variables to store the mean and sigma values for pions
    double mean_pion, sigma_pion, momentum_range;
    std::vector<double> pion_means, pion_sigmas, momentum_ranges;

    // Read the file line by line
    std::string line;
    while (std::getline(fitFile, line)) {
        std::istringstream ss(line);
        double value1, value2, value3, value4, value5, value6;
        if (ss >> value1 >> value2 >> value3 >> value4 >> mean_pion >> sigma_pion) {
            // Store momentum range, mean and sigma for pions
            momentum_ranges.push_back(value1); // Assuming value1 is momentum range
            pion_means.push_back(mean_pion);
            pion_sigmas.push_back(sigma_pion);
        }
    }
    fitFile.close();

    // Now print the momentum range along with mean ± 3σ for pions
    std::cout << "Momentum Range (GeV)\tMean\tMean - 3σ\tMean + 3σ" << std::endl;
    for (size_t i = 0; i < pion_means.size(); ++i) {
        double mean = pion_means[i];
        double sigma = pion_sigmas[i];
        double lowerBound = mean - 3 * sigma;
        double upperBound = mean + 3 * sigma;

        // Print the momentum range and the corresponding mean ± 3σ
        std::cout << momentum_ranges[i] << "\t" << mean << "\t" << lowerBound << "\t" << upperBound << std::endl;
    }

    // Now plot the pion data points with 3-sigma range
    TGraph* pionData = new TGraph(2 * pion_means.size()); // 2 points per momentum range

    for (size_t i = 0; i < pion_means.size(); ++i) {
        double mean = pion_means[i];
        double sigma = pion_sigmas[i];
        double lowerBound = mean - 3 * sigma;
        double upperBound = mean + 3 * sigma;

        // Plot data points for mean + 3σ and mean - 3σ
        pionData->SetPoint(2 * i, momentum_ranges[i], lowerBound); // mean - 3σ
        pionData->SetPoint(2 * i + 1, momentum_ranges[i], upperBound); // mean + 3σ

        // Set small marker styles for these data points
        pionData->SetMarkerStyle(21); // Small square markers
        pionData->SetMarkerSize(0.7); // Smaller size for points
        pionData->SetMarkerColor(kBlack); // Color black for the points
    }

    pionData->Draw("P SAME"); // Add the data points to the existing plot

    // Generate theoretical curves
    const int nPoints = 500;
    double kaonMass = 0.4937;   // GeV/c^2
    double pionMass = 0.1396;   // GeV/c^2
    double protonMass = 0.9383; // GeV/c^2
    double momentum[nPoints];
    double betaKaon[nPoints];
    double betaPion[nPoints];
    double betaProton[nPoints];

    for (int i = 0; i < nPoints; ++i) {
        momentum[i] = i * 8.0 / (nPoints - 1); // Linear spacing from 0 to 8 GeV
        betaKaon[i] = momentum[i] / sqrt(kaonMass * kaonMass + momentum[i] * momentum[i]);
        betaPion[i] = momentum[i] / sqrt(pionMass * pionMass + momentum[i] * momentum[i]);
        betaProton[i] = momentum[i] / sqrt(protonMass * protonMass + momentum[i] * momentum[i]);
    }

    TGraph* kaonCurve = new TGraph(nPoints, momentum, betaKaon);
    kaonCurve->SetLineColor(kRed);
    kaonCurve->SetLineWidth(1); // Thinner line for Kaon
    kaonCurve->Draw("L SAME");

    TGraph* pionCurve = new TGraph(nPoints, momentum, betaPion);
    pionCurve->SetLineColor(kBlack);
    pionCurve->SetLineWidth(1); // Thinner line for Pion
    pionCurve->Draw("L SAME");

    TGraph* protonCurve = new TGraph(nPoints, momentum, betaProton);
    protonCurve->SetLineColor(kBlue + 2); // Green for Proton
    protonCurve->SetLineWidth(1);          // Thinner line for Proton
    protonCurve->Draw("L SAME");

    // Add legend for clarity
    auto legend = new TLegend(0.8, 0.8, 0.9, 0.9);
    legend->AddEntry(pionCurve, "Pion", "l");
    legend->AddEntry(kaonCurve, "Kaon", "l");
    legend->AddEntry(protonCurve, "Proton", "l");
    legend->SetBorderSize(0);
    legend->SetTextSize(0.02);
    legend->Draw();

    canvas->SaveAs("output_with_all_curves_and_pion_data_points.pdf");

    timer.Stop();
    std::cout << "Execution completed in " << timer.RealTime() << " seconds (real time), "
              << timer.CpuTime() << " seconds (CPU time)." << std::endl;

    file->Close();
}
