#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TFitResult.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>

using namespace std;

double fitFunction(double *x, double *par) {
    if (x[0] <= 2.15) return 5.0; // Constant 5 up to 2.15 GeV/c
    return par[0] + par[1] * exp(-x[0] / par[2]) + par[3] * exp(-x[0] / par[4]); // Double exponential for p > 2.15
}

void mid() {
    // Load threshold data from file
    vector<double> pValues, thresholds;
    ifstream inFile("chi2pid_thresholds.txt");
    if (!inFile) {
        cout << "Error: Could not open chi2pid_thresholds.txt" << endl;
        return;
    }
    string line;
    getline(inFile, line); // Skip header
    while (getline(inFile, line)) {
        stringstream ss(line);
        string token;
        getline(ss, token, ',');
        pValues.push_back(stod(token));
        getline(ss, token, ',');
        thresholds.push_back(stod(token));
    }
    inFile.close();

    // Create TGraph
    TCanvas *canvas = new TCanvas("canvas", "Chi2pid Threshold vs Momentum", 800, 600);
    TGraph *graph = new TGraph(pValues.size(), pValues.data(), thresholds.data());
    graph->SetTitle("");
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlue);
    graph->GetXaxis()->SetRangeUser(0.2, 8.5);
    graph->GetYaxis()->SetRangeUser(0.0, 5.5);
    graph->GetYaxis()->SetTitle("Cut Position (chi2pid)");
    graph->GetXaxis()->SetTitle("Momentum (GeV/c)");
    graph->Draw("AP");

    // Define fit function with adjusted range
    TF1 *fitFunc = new TF1("fitFunc", fitFunction, 2.15, 8.0, 5); // Fit only for p > 2.15
    fitFunc->SetParameters(-0.21, 3.05, 4.51, 36.23, 0.87); // Use previous fit results
    fitFunc->SetParNames("Offset", "Amp1", "Decay1", "Amp2", "Decay2");
    
    // Set initial step sizes
    fitFunc->SetParLimits(0, -1.0, 1.0); // Offset
    fitFunc->SetParLimits(1, 0.0, 10.0); // Amp1
    fitFunc->SetParLimits(2, 0.1, 10.0); // Decay1
    fitFunc->SetParLimits(3, 0.0, 50.0); // Amp2
    fitFunc->SetParLimits(4, 0.1, 2.0); // Decay2

    // Perform fit and get result
    TFitResultPtr fitResult = graph->Fit(fitFunc, "RMWS"); // 'R' for range, 'M' for improved Migrad, 'W' for weighting, 'S' to store result
    if (fitResult.Get() == nullptr || !fitResult->IsValid()) {
        cout << "Fit failed or is invalid" << endl;
        return;
    }

    // Ensure full range is drawn
    fitFunc->SetRange(0.2, 8.0);
    fitFunc->SetLineColor(kRed);
    fitFunc->Draw("SAME");

    // Manually calculate chi2 (optional, for verification)
    double chi2 = 0.0;
    int nPoints = pValues.size();
    for (int i = 0; i < nPoints; i++) {
        if (pValues[i] > 2.15) {
            double x = pValues[i];
            double y = thresholds[i];
            if (y > 0) {
                double fitValue = fitFunc->Eval(x);
                double weight = 1.0 / (y * y);
                chi2 += (y - fitValue) * (y - fitValue) * weight;
            }
        }
    }

    // Get NDF and parameters from fit result
    int ndf = fitResult->Ndf();
    double chi2FromFit = fitResult->Chi2();

    // Save results to text file
    ofstream outFile("fit_results.txt");
    if (outFile.is_open()) {
        outFile << fixed << setprecision(6);
        outFile << "Fit Results\n";
        outFile << "==========\n";
        outFile << "Parameters:\n";
        for (int i = 0; i < 5; i++) {
            outFile << fitFunc->GetParName(i) << " = " << fitResult->Parameter(i)
                    << " +/- " << fitResult->ParError(i) << "\n";
        }
        outFile << "Chi2 (Manual) = " << chi2 << "\n";
        outFile << "Chi2 (Fit) = " << chi2FromFit << "\n";
        outFile << "NDF = " << ndf << "\n";
        outFile << "Chi2/NDF (Manual) = " << chi2 / ndf << "\n";
        outFile << "Chi2/NDF (Fit) = " << chi2FromFit / ndf << "\n";
        outFile.close();
        cout << "Fit results saved to fit_results.txt" << endl;
    } else {
        cout << "Error: Could not open fit_results.txt for writing" << endl;
    }

    // Add manual text using TLatex
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.02);
    latex->DrawLatex(0.5, 0.79, "Fit: 5 for p <= 2.15 GeV/c");
    latex->DrawLatex(0.5, 0.76, Form("%.2f + %.2f*exp(-p/%.2f) + %.2f*exp(-p/%.2f) for p > 2.15 GeV/c",
                                    fitResult->Parameter(0), fitResult->Parameter(1),
                                    fitResult->Parameter(2), fitResult->Parameter(3),
                                    fitResult->Parameter(4)));

    // Save as PNG
    canvas->SaveAs("chi2pid_threshold_vs_p_fit.png");
    cout << "PNG saved as: chi2pid_threshold_vs_p_fit.png" << endl;

    // Clean up
    delete graph;
    delete fitFunc;
    delete latex;
    delete canvas;
}
