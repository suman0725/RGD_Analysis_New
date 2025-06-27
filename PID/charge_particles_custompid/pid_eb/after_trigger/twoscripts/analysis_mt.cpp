#include <iostream>
#include <vector>
#include <atomic> // For thread-safe counter
#include <thread> // For hardware_concurrency
#include <filesystem> // For directory creation
#include <iomanip>
#include "TFile.h"
#include "TH1F.h" // For 1D histograms
#include "TH2F.h"
#include "TGraph.h" // For expected beta curves
#include "TCanvas.h"
#include "TPad.h"
#include "TStopwatch.h" // For timing
#include "TLegend.h" // For legends in combined plots
#include "TF1.h" // For Gaussian fitting
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TStyle.h"
#include "TLatex.h"
#include <fstream>

using namespace std;
using namespace ROOT;
std::string formatDouble(double value) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(3) << value;
    return stream.str();
}



// Function to compute momentum from px, py, pz
float computeMomentum(const RVec<float>& px, const RVec<float>& py, const RVec<float>& pz, size_t index) {
    return sqrt(px[index] * px[index] + py[index] * py[index] + pz[index] * pz[index]);
}

// Function to compute expected beta for a given mass and momentum
float computeExpectedBeta(float p, float mass) {
    return p / sqrt(p * p + mass * mass);
}

void AddFitLegend(TCanvas* c, TH1* h, const char* label, double mean, double sigma) {
    c->cd();
    TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(h, label, "l");
    leg->AddEntry((TObject*)0, Form("Mean = %.3f", mean), "");
    leg->AddEntry((TObject*)0, Form("Sigma = %.3f", sigma), "");
    leg->Draw();
}


int main() {
    // Determine the number of threads supported by the CPU
    unsigned int nThreads = std::thread::hardware_concurrency();
    if (nThreads == 0) {
        cout << "Could not determine the number of threads; defaulting to 4." << endl;
        nThreads = 4; // Fallback value if hardware_concurrency fails
    }
    cout << "Number of threads supported by CPU: " << nThreads << endl;

    // Enable multi-threading in ROOT using the detected number of threads
    ROOT::EnableImplicitMT(nThreads);
    cout << "Multi-threading enabled for RDataFrame processing with " << nThreads << " threads." << endl;

    // Open the ROOT file and create an RDataFrame from the TTree
    TFile* inputFile = new TFile("output/particle_data.root", "READ");
    if (!inputFile->IsOpen()) {
        cerr << "Error: Could not open output/particle_data.root!" << endl;
        return 1;
    }

    // Create RDataFrame from the TTree
    ROOT::RDataFrame df("ParticleTree", inputFile);

   
    Long64_t totalEvents = *df.Count();
    cout << "Total number of events to process: " << totalEvents << endl;

    
    //TH2F hChi2PIDvsPPID211("hChi2PIDvsPPID211", "Chi2PID vs. Momentum for Positive Pions (PID 211, Forward);Chi2PID;Momentum p (GeV)", 200, -10, 10, 200, 0, 10);
    TH2F hChi2PIDvsPPID211("hChi2PIDvsPPID211", "Chi2PID vs. Momentum for Positive Pions (PID 211, Forward);Momentum p (GeV);Chi2PID", 200, 0, 10, 200, -5, 5);
    TH2F hChi2PIDvsPPID211_MeanBetaCut("hChi2PIDvsPPID211_MeanBetaCut", "Chi2PID vs. Momentum for Positive Pions (PID 211, Forward, Beta > Mean(Pion,Kaon));Momentum p (GeV);Chi2PID", 200, 0, 10, 200, -5, 5);
    TH2F hChi2PIDvsPPID321("hChi2PIDvsPPID321", "Chi2PID vs. Momentum for Positive Kaons (PID 321, Forward);Chi2PID;Momentum p (GeV)", 200, -5, 5, 200, 0, 10);
    TH2F hChi2PIDvsPPID2212("hChi2PIDvsPPID2212", "Chi2PID vs. Momentum for Protons (PID 2212, Forward);Chi2PID;Momentum p (GeV)", 200, -5, 5, 200, 0, 10);
    TH2F hBetaVsPPositive("hBetaVsPPositive", "Beta vs. Momentum for Positive Charge Particles (Excl. PID -11, Forward);Momentum p (GeV);Beta", 200, 0, 4, 200, 0.8, 1.05);
    TH2F hBetaVsPPKP("hBetaVsPPKP", "Beta vs. Momentum for EB (211, 321, 2212) Forward);Momentum p (GeV);Beta", 200, 0, 4, 200, 0.8, 1.05);
    TH2F hBetaVsPPion("hBetaVsPPion", "Beta vs. Momentum for EB PID 211: Forward; Momentum p (GeV); Beta", 200,0, 10, 200, 0.8, 1.05); 
    TH2F hBetaVsPPion_MeanBetaCut("hBetaVsPPion_MeanBetaCut", "Beta vs. Momentum for EB PID 211: Forward, Beta > Mean(Pion,Kaon); Momentum p (GeV); Beta", 200,0, 10, 200, 0.8, 1.05); 


    // Define momentum slices (1 to 8 GeV, width 0.35 GeV)
    const double pMin = 1.0;
    const double pMax = 8;
    const double pWidth = 0.35;
    const int nSlices = static_cast<int>((pMax - pMin) / pWidth); // 20 slices
    const int sliceThreshold = 5; // Slices 0-4 (up to 2.05-2.40 GeV) use beta 0.8-1.05, slices 5-19 use beta 0.95-1.05

    const double massPion = 0.13957;  // GeV/c²
    const double massKaon = 0.49367;  // GeV/c²
    const double massProton = 0.93827; // GeV/c²

    // Define Gaussian fit ranges for each slice (kaon, pion)
    const double kaonFitRanges[20][2] = {
        {0.9, 0.94},    // Slice 1
        {0.935, 0.97},  // Slice 2
        {0.955, 0.98},  // Slice 3
        {0.965, 0.99},  // Slice 4
        {0.975, 0.99},  // Slice 5
        {0.98, 0.994},  // Slice 6
        {0.982, 0.996}, // Slice 7
        {0.984, 0.998}, // Slice 8
        {0.986, 0.998}, // Slice 9
        {0.988, 0.998}, // Slice 10
        {0.99, 0.998},  // Slice 11
        {0.99, 0.998},  // Slice 12
        {0.992, 0.999}, // Slice 13
        {0.992, 0.999}, // Slice 14
        {0.993, 0.999}, // Slice 15
        {0.994, 0.999}, // Slice 16
        {0.994, 0.999}, // Slice 17
        {0.995, 0.999}, // Slice 18
        {0.995, 0.999}, // Slice 19
        {0.995, 0.999}  // Slice 20 (assumed)
    };

    const double pionFitRanges[20][2] = {
        {0.98, 1.05},   // Slice 1
        {0.985, 1.05},  // Slice 2
        {0.99, 1.05},   // Slice 3
        {0.99, 1.05},   // Slice 4
        {0.99, 1.05},   // Slice 5
        {0.994, 1.008}, // Slice 6
        {0.994, 1.006}, // Slice 7
        {0.994, 1.006}, // Slice 8
        {0.994, 1.006}, // Slice 9
        {0.996, 1.006}, // Slice 10
        {0.996, 1.006}, // Slice 11
        {0.996, 1.006}, // Slice 12
        {0.996, 1.006}, // Slice 13
        {0.996, 1.006}, // Slice 14
        {0.996, 1.006}, // Slice 15
        {0.996, 1.005}, // Slice 16
        {0.996, 1.005}, // Slice 17
        {0.996, 1.005}, // Slice 18
        {0.996, 1.005}, // Slice 19
        {0.996, 1.005}  // Slice 20 (assumed)
    };

    // Create 1D histograms for beta in each momentum slice for pions and kaons
    vector<TH1F*> hBetaPionSlices(nSlices);
    vector<TH1F*> hBetaKaonSlices(nSlices);
    for (int i = 0; i < nSlices; ++i) {
        double pLow = pMin + i * pWidth;
        double pHigh = pLow + pWidth;
        string namePion = "hBetaPion_p_" + to_string(pLow) + "_" + to_string(pHigh);
        string nameKaon = "hBetaKaon_p_" + to_string(pLow) + "_" + to_string(pHigh);
        string titlePion = "Beta for Pions (p: " + to_string(pLow) + " to " + to_string(pHigh) + " GeV);Beta;Counts";
        string titleKaon = "Beta for Kaons (p: " + to_string(pLow) + " to " + to_string(pHigh) + " GeV);Beta;Counts";
        // Use different beta ranges based on the slice index
        if (i < sliceThreshold) {
            // Slices 0-4: beta range 0.8 to 1.05
            hBetaPionSlices[i] = new TH1F(namePion.c_str(), titlePion.c_str(), 100, 0.8, 1.05);
            hBetaKaonSlices[i] = new TH1F(nameKaon.c_str(), titleKaon.c_str(), 100, 0.8, 1.05);
        } else {
            // Slices 5-19: beta range 0.95 to 1.05
            hBetaPionSlices[i] = new TH1F(namePion.c_str(), titlePion.c_str(), 100, 0.95, 1.05);
            hBetaKaonSlices[i] = new TH1F(nameKaon.c_str(), titleKaon.c_str(), 100, 0.95, 1.05);
        }
    }

    // Create 1D histograms for chi2pid in each momentum slice for pions (PID 211) and kaons (PID 321)
    vector<TH1F*> hChi2PIDSlices(nSlices);
    vector<TH1F*> hChi2PidKaons(nSlices); 
    vector<TH1F*> hChi2PIDSlicesMeanBetaCut(nSlices);

    for (int i = 0; i < nSlices; ++i) {
        double pLow = pMin + i * pWidth;
        double pHigh = pLow + pWidth;
        // Chi2PID histograms for pions
        string nameChi2PID = "hChi2PID_p_" + to_string(pLow) + "_" + to_string(pHigh);
        string titleChi2PID = "Chi2PID for Pions (p: " + formatDouble(pLow) + " to " + formatDouble(pHigh) + " GeV);Chi2PID;Counts";
        hChi2PIDSlices[i] = new TH1F(nameChi2PID.c_str(), titleChi2PID.c_str(), 100, -3.2, 3.2);

        // Kaons (321)
        string nameKaon = "hChi2PidKaon_p_" + to_string(pLow) + "_" + to_string(pHigh);
        string titleKaon = "Chi2PID for Kaons (p: " + to_string(pLow) + " to " + to_string(pHigh) + " GeV);Chi2PID;Counts";
        hChi2PidKaons[i] = new TH1F(nameKaon.c_str(), titleKaon.c_str(), 100, -5, 5);

        // Chi2PID histograms for pions with mean beta cut
        string nameChi2PIDMeanBetaCut = "hChi2PID_MeanBetaCut_p_" + to_string(pLow) + "_" + to_string(pHigh);
        string titleChi2PIDMeanBetaCut = "Chi2PID for Pions (p: " + formatDouble(pLow) + " to " + formatDouble(pHigh) + " GeV, Beta > Mean(Pion,Kaon));Chi2PID;Counts";
        hChi2PIDSlicesMeanBetaCut[i] = new TH1F(nameChi2PIDMeanBetaCut.c_str(), titleChi2PIDMeanBetaCut.c_str(), 100, -3.2, 3.2);
    }


    // Thread-safe counter for processed events
    std::atomic<Long64_t> processedEvents(0);
    Long64_t progressInterval = totalEvents / 20; // Update progress every 5%
    if (progressInterval == 0) progressInterval = 1; // Avoid division by zero for small datasets

    // Set up a timer to measure the processing time
    TStopwatch timer;
    timer.Start();

    // Process the data in parallel using RDataFrame with Foreach
    df.Foreach([&](const RVec<int>& pid, const RVec<float>& chi2pid, const RVec<float>& beta, const RVec<float>& px, const RVec<float>& py, const RVec<float>& pz, const RVec<int>& charge, const RVec<int>& status) {
        // Increment the counter for each event (thread-safe)
        Long64_t currentEvent = processedEvents.fetch_add(1, std::memory_order_relaxed);

        // Print progress at regular intervals
        if (currentEvent % progressInterval == 0 || currentEvent == totalEvents - 1) {
            double percent = (currentEvent + 1) * 100.0 / totalEvents;
            cout << "\rProgress: " << fixed << setprecision(1) << percent << "%" << flush;
        }

        // Process the event
        for (size_t i = 0; i < pid.size(); ++i) {
            // Check if the particle is in the forward region: abs(status)/2000 == 1
            if (abs(status[i]) / 2000 != 1) continue;

            float p = computeMomentum(px, py, pz, i);

            // Fill chi2pid vs. momentum for pions (PID 211)
            if (pid[i] == 211 && chi2pid[i] != 9999.0) {
               // hChi2PIDvsPPID211.Fill(chi2pid[i], p); 
                hChi2PIDvsPPID211.Fill(p,chi2pid[i]); 
               
                // Fill 1D chi2pid histograms for pions in momentum slices
                if (p >= pMin && p < pMax) {
                    int sliceIndex = static_cast<int>((p - pMin) / pWidth);
                    hChi2PIDSlices[sliceIndex]->Fill(chi2pid[i]);
                }
                // Fill histograms with mean beta cut
                if (beta[i] != 9999.0) {
                    float expectedBetaPion = computeExpectedBeta(p, massPion);
                    float expectedBetaKaon = computeExpectedBeta(p, massKaon);
                    float meanBeta = 0.5 * (expectedBetaPion + expectedBetaKaon);
                    if (beta[i] > meanBeta) {
                        hChi2PIDvsPPID211_MeanBetaCut.Fill(p, chi2pid[i]);
                        // Fill 1D chi2pid histograms for pions in momentum slices (with beta cut)
                        if (p >= pMin && p < pMax) {
                            int sliceIndex = static_cast<int>((p - pMin) / pWidth);
                            hChi2PIDSlicesMeanBetaCut[sliceIndex]->Fill(chi2pid[i]);
                        }
                    }
                }
            }

            // Fill chi2pid vs. momentum for kaons (PID 321)
            if (pid[i] == 321 && chi2pid[i] != 9999.0) {
                hChi2PIDvsPPID321.Fill(chi2pid[i], p); // chi2pid on x-axis, p on y-axis
                if (p >= pMin && p < pMax) {
                    int sliceIndex = static_cast<int>((p - pMin) / pWidth);
                    hChi2PidKaons[sliceIndex]->Fill(chi2pid[i]);
                }
            }

            // Fill chi2pid vs. momentum for protons (PID 2212)
            if (pid[i] == 2212 && chi2pid[i] != 9999.0) {
                hChi2PIDvsPPID2212.Fill(chi2pid[i], p); // chi2pid on x-axis, p on y-axis
                
            }
            

            // Fill beta vs. momentum for positive charge particles
            if (charge[i] > 0 && pid[i]!= -11 && beta[i]!=9999.0) { // Use the charge branch directly
                hBetaVsPPositive.Fill(p, beta[i]); // p on x-axis, beta on y-axis
            }

            if (pid[i] == 211 || pid[i] == 321 || pid[i] == 2212 && beta[i]!=9999.0){
                if (pid[i] == 211){ 
                    hBetaVsPPion.Fill(p, beta[i]);
                    float expectedBetaPion = computeExpectedBeta(p, massPion);
                    float expectedBetaKaon = computeExpectedBeta(p, massKaon);
                    float meanBeta = 0.5 * (expectedBetaPion + expectedBetaKaon);
                    if (beta[i] > meanBeta) {
                        hBetaVsPPion_MeanBetaCut.Fill(p, beta[i]);
                    }
                }
                hBetaVsPPKP.Fill(p, beta[i]);

            }
           
            

            // Fill 1D beta histograms for pions and kaons in momentum slices
            if (charge[i] > 0 && (pid[i] == 211 || pid[i] == 321) && p >= pMin && p < pMax) {
                int sliceIndex = static_cast<int>((p - pMin) / pWidth);
                // Only fill if beta is within the histogram's range
                double betaMin = (sliceIndex < sliceThreshold) ? 0.8 : 0.95;
                double betaMax = 1.05;
                if (beta[i] >= betaMin && beta[i] <= betaMax) {
                    if (pid[i] == 211) {
                        hBetaPionSlices[sliceIndex]->Fill(beta[i]);
                    } else if (pid[i] == 321) {
                        hBetaKaonSlices[sliceIndex]->Fill(beta[i]);
                    }
                }
            }
        }
    }, {"pid", "chi2pid", "beta", "px", "py", "pz", "charge", "status"});

    // Ensure the progress line ends with a newline
    cout << endl;

    // Stop the timer and print the elapsed time in minutes
    timer.Stop();
    double timeInMinutes = timer.RealTime() / 60.0; // Convert seconds to minutes
    cout << "Time taken for RDataFrame processing: " << timeInMinutes << " minutes" << endl;

    // Debug: Check entries in the 2D histograms
    cout << "Entries in hChi2PIDvsPPID211 (Pions, Forward): " << hChi2PIDvsPPID211.GetEntries() << endl;
    cout << "Entries in hChi2PIDvsPPID321 (Kaons, Forward): " << hChi2PIDvsPPID321.GetEntries() << endl;
    cout << "Entries in hChi2PIDvsPPID2212 (Protons, Forward): " << hChi2PIDvsPPID2212.GetEntries() << endl;
    cout << "Entries in hBetaVsPPositive (Positive Charge Particles, Forward): " << hBetaVsPPositive.GetEntries() << endl;

    // Debug: Check entries in the 1D beta histograms
    for (int i = 0; i < nSlices; ++i) {
        double pLow = pMin + i * pWidth;
        double pHigh = pLow + pWidth;
        cout << "Momentum slice " << pLow << " to " << pHigh << " GeV: "
             << "Pions: " << hBetaPionSlices[i]->GetEntries() << " entries, "
             << "Kaons: " << hBetaKaonSlices[i]->GetEntries() << " entries" 
             << "Pions (Chi2PID): " << hChi2PIDSlices[i]->GetEntries() << " entries" << endl;
    }

    // Create the output directory if it doesn't exist
    std::filesystem::create_directories("output_mt");
    cout << "Output directory 'output_mt' created (if it didn't exist)." << endl;

    // Create canvases
    TCanvas* canvasBetaVsP = new TCanvas("canvasBetaVsP", "Beta vs. Momentum for Positive Charge Particles (Forward)", 800, 600);
    TCanvas* canvasBetaVsPPKP = new TCanvas("canvasBetaVsPPKP", "Beta vs. Momentum for EB PID (211, 321, 2212) (Forward)", 800, 600);
    TCanvas* canvasBetaVsPPion = new TCanvas("canvasBetaVsPPion", "Beta vs. Momentum for EB PID 211 (Forward)", 800, 600);
    TCanvas* canvasBetaVsPPion_MeanBetaCut = new TCanvas("canvasBetaVsPPion_MeanBetaCut", "Beta vs. Momentum for EB PID 211 (Forward),Beta > Mean(Pion,Kaon)", 800, 600);
    if (hBetaVsPPositive.GetEntries()>0 ||hBetaVsPPKP.GetEntries()> 0 || hBetaVsPPion.GetEntries()>0 ||hBetaVsPPion_MeanBetaCut.GetEntries()>0) {
        // Draw first histogram

        canvasBetaVsP->cd();  
        hBetaVsPPositive.SetStats(0);
        hBetaVsPPositive.Draw("COLZ");
        gPad->SetLogz();  // Log scale

        // Draw second histogram
        canvasBetaVsPPKP->cd();  
        hBetaVsPPKP.SetStats(0);
        hBetaVsPPKP.Draw("COLZ");
        gPad->SetLogz();  // Log scale

        canvasBetaVsPPion->cd(); 
        hBetaVsPPion.SetStats(0); 
        hBetaVsPPion.Draw("COLZ");
        gPad->SetLogz(); 

        canvasBetaVsPPion_MeanBetaCut->cd();
        hBetaVsPPion_MeanBetaCut.SetStats(0); 
        hBetaVsPPion_MeanBetaCut.Draw("COLZ");
        gPad->SetLogz(); 

        // Add expected beta curves for pion, kaon, and proton
        const int nPoints = 1000;
        double pValues[nPoints], betaPion[nPoints], betaKaon[nPoints], betaProton[nPoints], meanbetaPK[nPoints ];
        const double pMinCurve = 0.0, pMaxCurve = 10.0;
        

        // Compute beta values
        for (int i = 0; i < nPoints; ++i) {
            pValues[i] = pMinCurve + (pMaxCurve - pMinCurve) * i / (nPoints - 1);
            betaPion[i] = computeExpectedBeta(pValues[i], massPion);
            betaKaon[i] = computeExpectedBeta(pValues[i], massKaon);
            betaProton[i] = computeExpectedBeta(pValues[i], massProton);
            meanbetaPK[i] = 0.5 * (betaPion[i] + betaKaon[i]);
        }

        // Create graphs
        TGraph* graphPion = new TGraph(nPoints, pValues, betaPion);
        TGraph* graphKaon = new TGraph(nPoints, pValues, betaKaon);
        TGraph* graphProton = new TGraph(nPoints, pValues, betaProton);
        TGraph* graphmeanPK = new TGraph(nPoints, pValues, meanbetaPK);

        // Style settings
        graphPion->SetLineColor(kRed);
        graphPion->SetLineWidth(1);
        graphPion->SetLineStyle(1);  // Solid line

        graphKaon->SetLineColor(kBlue);
        graphKaon->SetLineWidth(1);
        graphKaon->SetLineStyle(1);  // Solid line

        graphProton->SetLineColor(kGreen);
        graphProton->SetLineWidth(1);
        graphProton->SetLineStyle(1);  // Solid line

        graphmeanPK->SetLineColor(kBlack);
        graphmeanPK->SetLineWidth(1); 
        graphmeanPK->SetLineStyle(2); //dash line

        
        // Draw graphs on both canvases
        canvasBetaVsP->cd();
        graphPion->Draw("L SAME");
        graphKaon->Draw("L SAME");
        graphProton->Draw("L SAME");
        graphmeanPK->Draw("L SAME");

        canvasBetaVsPPKP->cd();
        graphPion->Draw("L SAME");
        graphKaon->Draw("L SAME");
        graphProton->Draw("L SAME");
        graphmeanPK->Draw("L SAME");

        canvasBetaVsPPion->cd(); 
        graphPion->Draw("L SAME");
        graphKaon->Draw("L SAME");
        graphProton->Draw("L SAME");
        graphmeanPK->Draw("L SAME");

        canvasBetaVsPPion_MeanBetaCut->cd();
        graphPion->Draw("L SAME");
        graphKaon->Draw("L SAME");
        graphProton->Draw("L SAME");
        graphmeanPK->Draw("L SAME");
        

        // Save to PDF
        canvasBetaVsP->Print("output_mt/beta_vs_p_positive_forward.pdf");
        canvasBetaVsPPKP->Print("output_mt/beta_vs_p_EBPKP.pdf");
        canvasBetaVsPPion->Print("output_mt/beta_vs_p_EBpion.pdf");
        canvasBetaVsPPion_MeanBetaCut->Print("output_mt/beta_vs_p_EBpion_meanbetacut.pdf");
    }




     // Plot: 1D Chi2PID Histograms for Pions in Momentum Slices (one canvas per slice, saved in a multi-page PDF)
    vector<double> chi2pidMeanPion(nSlices), chi2pidSigmaPion(nSlices);
    vector<double> chi2pidMeanPlus3Sigma(nSlices), chi2pidMeanMinus3Sigma(nSlices);
    vector<double> pValues(nSlices); 
    vector<double> chi2pidMeanKaon(nSlices), chi2pidSigmaKaon(nSlices);
    vector<double> chi2NdfPion(nSlices), chi2NdfKaon(nSlices);
    vector<double> distanceInSigma(nSlices);

    TCanvas* canvasChi2PIDSlices = new TCanvas("canvasChi2PIDSlices", "Chi2PID for Pions in Momentum Slice (Forward)", 800, 600);
    TCanvas* canvasChi2PIDSlicesKaon = new TCanvas("canvasChi2PIDSlicesKaon", "Chi2PID for Kaons in Momentum Slices (Forward)", 800, 600);
    TCanvas* canvasComparison = new TCanvas("canvasComparison", "Pion vs Kaon Chi2PID Comparison", 1200, 800);
    TCanvas* canvasChi2PIDSlicesAll = new TCanvas("canvasChi2PIDSlicesAll", "Chi2PID for All Momentum Slices", 1200, 800);
    canvasChi2PIDSlicesAll->Divide(5, 4); // 5 columns, 4 rows layout
    
    
    canvasChi2PIDSlices->Print("output_mt/chi2pid_slices_pion_forward.pdf["); 
    canvasChi2PIDSlicesKaon->Print("output_mt/chi2pid_slices_kaon_forward.pdf[");
    canvasComparison->Print("output_mt/chi2pid_comparison.pdf[");

    // Open a CSV file to store the results
    ofstream csvFile("fit_results.csv");
    csvFile << "Momentum (GeV),Pion Mean,Pion Sigma,Pion Chi2/ndf,Kaon Mean,Kaon Sigma,Kaon Chi2/ndf\n";
    for (int i = 0; i < nSlices; ++i) {
        canvasComparison->Clear();
        canvasChi2PIDSlices->Clear();
        canvasChi2PIDSlicesKaon->Clear();

        canvasChi2PIDSlices->cd();
        canvasChi2PIDSlices->SetGrid(); 
        hChi2PIDSlices[i]->SetStats(0);
        hChi2PIDSlices[i]->SetLineColor(kBlack);
        hChi2PIDSlices[i]->SetMarkerStyle(20);
        hChi2PIDSlices[i]->SetLineWidth(3);
        //hChi2PIDSlices[i]->Scale(1.0 / hChi2PIDSlices[i]->Integral());
        hChi2PIDSlices[i]->Draw("P");

       
        double fitMin = -1.6;
        double fitMax = 1.6;
        double sliceMin = pMin + i * pWidth;
        double sliceMax = sliceMin + pWidth;
        if (sliceMin >= 2.75 && sliceMin <= 3.45) {
            fitMin = -1;
            fitMax = 1;
        }
        if (sliceMin >3.50 && sliceMin<4.85) {  
            fitMin = -0.8;
            fitMax = 0.8;
        }
        
         if (sliceMin >= 4.85&& sliceMin <= 5.55) {
            fitMin = -0.4;
            fitMax = 0.5;
        }
        if (sliceMin >5.8 && sliceMin < 5.95) {
            fitMin = -0.4;
            fitMax = 0.4;
        }
        if (sliceMin >= 6.25) {
            fitMin = -0.4;
            fitMax = 0.35;
        }
       

        // Define the function using the updated fit range
        TF1 *gaussFitPion = new TF1("gaussFitPion", "[0] * exp(-0.5 * ((x-[1])/[2])^2)", fitMin, fitMax);
        gaussFitPion->SetParameters(hChi2PIDSlices[i]->GetMaximum(), hChi2PIDSlices[i]->GetMean(), hChi2PIDSlices[i]->GetRMS());
        hChi2PIDSlices[i]->Fit(gaussFitPion, "R");

        // Extract mean and sigma from the fit
        double fitMean = gaussFitPion->GetParameter(1);
        double fitSigma = gaussFitPion->GetParameter(2);
        chi2pidMeanPion[i] = gaussFitPion->GetParameter(1);
        chi2pidSigmaPion[i] = gaussFitPion->GetParameter(2);
        chi2pidMeanPlus3Sigma[i] = chi2pidMeanPion[i] + 3 * chi2pidSigmaPion[i];
        chi2pidMeanMinus3Sigma[i] = chi2pidMeanPion[i] - 3 * chi2pidSigmaPion[i];

        AddFitLegend(canvasChi2PIDSlices, hChi2PIDSlices[i], "Pions", chi2pidMeanPion[i], chi2pidSigmaPion[i]);

        // Draw the fit
        gaussFitPion->SetLineColor(kBlack);
        gaussFitPion->Draw("SAME");


        // ============ KAON ANALYSIS (new code) ============
        canvasChi2PIDSlicesKaon->cd();
        hChi2PidKaons[i]->SetStats(0);
        hChi2PidKaons[i]->SetLineColor(kBlack);  // Different color for kaons
        //hChi2PidKaons[i]->Scale(1.0/hChi2PidKaons[i]->Integral());
        hChi2PidKaons[i]->Draw();
        TF1* gaussFitKaon = new TF1("gaussFitKaon", "gaus", -1.6, 1.6);
        gaussFitKaon->SetParameters(hChi2PidKaons[i]->GetMaximum(), 2, 1);  // Different initial guess
        hChi2PidKaons[i]->Fit(gaussFitKaon, "RQN");
        
        chi2pidMeanKaon[i] = gaussFitKaon->GetParameter(1);
        chi2pidSigmaKaon[i] = gaussFitKaon->GetParameter(2);

        AddFitLegend(canvasChi2PIDSlicesKaon, hChi2PidKaons[i], "Kaons", chi2pidMeanKaon[i], chi2pidSigmaKaon[i]);
        // Kaon plot (new)
        
        gaussFitKaon->SetLineColor(kBlack);
        gaussFitKaon->Draw("SAME");

        
        canvasComparison->cd();
        // Configure and draw pions
        hChi2PIDSlices[i]->SetStats(0);
        hChi2PIDSlices[i]->SetLineColor(kRed);
        //hChi2PIDSlices[i]->SetLineWidth(1);
        hChi2PIDSlices[i]->Scale(1.0/hChi2PIDSlices[i]->Integral());
        hChi2PIDSlices[i]->GetXaxis()->SetTitle("chi2pid");
        hChi2PIDSlices[i]->GetYaxis()->SetTitle("Normalized Counts");
        hChi2PIDSlices[i]->Draw("HIST");

        // Configure and draw kaons
        hChi2PidKaons[i]->SetStats(0);
        hChi2PidKaons[i]->SetLineColor(kBlue);
        hChi2PidKaons[i]->SetLineWidth(1);
        hChi2PidKaons[i]->Scale(1.0/hChi2PidKaons[i]->Integral());
        hChi2PidKaons[i]->Draw("HIST SAME");

        // Set y-axis range
        double ymax = max(hChi2PIDSlices[i]->GetMaximum(), hChi2PidKaons[i]->GetMaximum());
        hChi2PIDSlices[i]->GetYaxis()->SetRangeUser(0, ymax*1.2);

        // Fit pions
       
        gaussFitPion->SetParameters(hChi2PIDSlices[i]->GetMaximum(), 0.0, 1.0);
        hChi2PIDSlices[i]->Fit(gaussFitPion, "RQN");
        gaussFitPion->SetLineColor(kRed);
        gaussFitPion->Draw("SAME");

        // Fit kaons
        
        gaussFitKaon->SetParameters(hChi2PidKaons[i]->GetMaximum(), 1.0, 1.5);
        hChi2PidKaons[i]->Fit(gaussFitKaon, "RQN");
        gaussFitKaon->SetLineColor(kBlue);
        gaussFitKaon->Draw("SAME");

        // Store fit results
        chi2pidMeanPion[i] = gaussFitPion->GetParameter(1);
        chi2pidSigmaPion[i] = gaussFitPion->GetParameter(2);
        chi2pidMeanKaon[i] = gaussFitKaon->GetParameter(1);
        chi2pidSigmaKaon[i] = gaussFitKaon->GetParameter(2);

        // Compute chi2/ndf for both fits
        chi2NdfPion[i] = gaussFitPion->GetChisquare() / gaussFitPion->GetNDF();
        chi2NdfKaon[i] = gaussFitKaon->GetChisquare() / gaussFitKaon->GetNDF();
        pValues[i] =  pMin + (i + 0.5) * pWidth;

        if (pValues[i] <= 4.0 && chi2pidSigmaPion[i] > 0) {
            distanceInSigma[i] = fabs(chi2pidMeanKaon[i] - chi2pidMeanPion[i]) / chi2pidSigmaPion[i];
        } 
        // Write to CSV
        csvFile << pValues[i] << ","
                << chi2pidMeanPion[i] << ","
                << chi2pidSigmaPion[i] << ","
                << chi2NdfPion[i] << ","
                << chi2pidMeanKaon[i] << ","
                << chi2pidSigmaKaon[i] << ","
                << chi2NdfKaon[i] << ","
                << distanceInSigma[i] << "\n";

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.SetTextColor(kRed);
        latex.DrawLatex(0.15, 0.85, Form("#mu_{#pi} = %.2f, #sigma_{#pi} = %.2f, #chi^{2}/ndf_{#pi} = %.2f", chi2pidMeanPion[i], chi2pidSigmaPion[i], chi2NdfPion[i]));
        latex.SetTextColor(kBlue);
        latex.DrawLatex(0.15, 0.80, Form("#mu_{K} = %.2f, #sigma_{K} = %.2f, #chi^{2}/ndf_{K} = %.2f", chi2pidMeanKaon[i], chi2pidSigmaKaon[i], chi2NdfKaon[i]));

        
        canvasChi2PIDSlicesAll->cd(i + 1);
        //gPad->SetGrid(0);
        hChi2PIDSlices[i]->SetLineColor(kBlack);
        hChi2PIDSlices[i]->Draw("HIST");
        TF1 *gaussFitPionClone = (TF1*)gaussFitPion->Clone(Form("gaussFitPion_%d", i));
        gaussFitPionClone->SetLineColor(kRed);
        gaussFitPionClone->Draw("SAME");
        // Add a legend with mean and sigma
        TLegend *legend = new TLegend(0.55, 0.65, 0.85, 0.85); // (x1, y1, x2, y2) in NDC coordinates
        legend->SetBorderSize(0); // No border
        legend->SetFillStyle(0);  // Transparent background
       /*  legend->AddEntry(hChi2PIDSlices[i], "Pions", "l");
        legend->AddEntry(gaussFitPionClone, "Gaussian Fit", "l"); */
        legend->AddEntry((TObject*)0, Form("#mu = %.2f", chi2pidMeanPion[i]), "");
        legend->AddEntry((TObject*)0, Form("#sigma = %.2f", chi2pidSigmaPion[i]), "");
        legend->Draw();
        
       
        // Save comparison plot
        canvasComparison->Print("output_mt/chi2pid_comparison.pdf");
        canvasChi2PIDSlices->Print("output_mt/chi2pid_slices_pion_forward.pdf");
        canvasChi2PIDSlicesKaon->Print("output_mt/chi2pid_slices_kaon_forward.pdf");
        canvasChi2PIDSlicesAll->Print("output_mt/chi2pid_slices_pion_forward_same.pdf");
        
        
        // Cleanup
        delete gaussFitPion;
        delete gaussFitKaon;
       
        

    }
    // Close the CSV file
    csvFile.close();

    TCanvas* canvasGapVsP = new TCanvas("canvasGapVsP", "Distance in Sigma vs Momentum", 800, 600);
    TGraph* graphGap = new TGraph(nSlices, pValues.data(), distanceInSigma.data());

    graphGap->SetTitle("Distance in Sigma Between Pions and Kaons vs Momentum");
    graphGap->GetXaxis()->SetTitle("p (GeV/c)");
    graphGap->GetYaxis()->SetTitle("Distance in Sigma (#sigma_{#pi})");
    graphGap->SetMarkerStyle(20);  // Dots
    graphGap->SetMarkerColor(kBlue);
    graphGap->SetLineColor(kBlue);
    graphGap->Draw("APL");  // A for axis, P for points, L for line connecting points

    // Save the plot
    canvasGapVsP->Print("output_mt/distanceInSigmaVsP.pdf");

    TCanvas* canvasSigmaVsP = new TCanvas("canvasSigmaVsP", "Sigma vs Momentum", 800, 600);
    TGraph* graph = new TGraph(nSlices, pValues.data(), chi2pidSigmaPion.data());

    // Style it
    graph->SetTitle("Pion Sigma vs Momentum");
    graph->GetXaxis()->SetTitle("p (GeV/c)");
    graph->GetYaxis()->SetTitle("#sigma_{#pi}");
    graph->SetMarkerStyle(20);  // Dots
    graph->SetMarkerColor(kRed);

    // Draw
    graph->Draw("AP");

    canvasSigmaVsP->Print("output_mt/pionsigmaVsP.pdf");
    canvasChi2PIDSlices->Print("output_mt/chi2pid_slices_pion_forward.pdf]");
    canvasChi2PIDSlicesKaon->Print("output_mt/chi2pid_slices_kaon_forward.pdf]");
    canvasComparison->Print("output_mt/chi2pid_comparison.pdf]");



    // Plot: 1D Chi2PID Histograms for Pions with Mean Beta Cut in Momentum Slices
    TCanvas* canvasChi2PIDSlicesMeanBetaCut = new TCanvas("canvasChi2PIDSlicesMeanBetaCut", "Chi2PID for Pions in Momentum Slice (Forward, Beta > Mean(Pion,Kaon))", 800, 600);
    canvasChi2PIDSlicesMeanBetaCut->Print("output_mt/chi2pid_slices_pion_mean_beta_cut_forward.pdf["); // Open the PDF
    for (int i = 0; i < nSlices; ++i) {
        canvasChi2PIDSlicesMeanBetaCut->Clear();
        hChi2PIDSlicesMeanBetaCut[i]->SetStats(0);
        hChi2PIDSlicesMeanBetaCut[i]->SetLineColor(kRed);
        hChi2PIDSlicesMeanBetaCut[i]->Draw();

        // Fit a Gaussian to the histogram using the same range as for hChi2PIDSlices
        double fitMin = -1.6;
        double fitMax = 1.6;
        double sliceMin = pMin + i * pWidth;
        if (sliceMin >= 3.45) {
            fitMin = -0.8;
            fitMax = 0.8;
        }
       

        TF1* gaussFit = new TF1("gaussFit", "gaus", fitMin, fitMax);
        gaussFit->SetParameters(hChi2PIDSlicesMeanBetaCut[i]->GetMaximum(), hChi2PIDSlicesMeanBetaCut[i]->GetMean(), hChi2PIDSlicesMeanBetaCut[i]->GetRMS());
        hChi2PIDSlicesMeanBetaCut[i]->Fit(gaussFit, "RQN"); // R: range, Q: quiet, N: no draw

        // Extract mean and sigma from the fit
        double fitMean = gaussFit->GetParameter(1);
        double fitSigma = gaussFit->GetParameter(2);

        // Draw the fit
        gaussFit->SetLineColor(kBlack);
        gaussFit->Draw("SAME");

        // Add a legend with the fit parameters
        TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg->AddEntry(hChi2PIDSlicesMeanBetaCut[i], "Pions", "l");
        char fitText[100];
        sprintf(fitText, "Mean = %.4f", fitMean);
        leg->AddEntry((TObject*)0, fitText, "");
        sprintf(fitText, "Sigma = %.4f", fitSigma);
        leg->AddEntry((TObject*)0, fitText, "");
        leg->Draw();

        if (hChi2PIDSlicesMeanBetaCut[i]->GetEntries() == 0) {
            cout << "Warning: Chi2PID histogram for pions with mean beta cut in momentum slice " << (pMin + i * pWidth) << " to " << (pMin + (i + 1) * pWidth) << " GeV is empty." << endl;
        }

        canvasChi2PIDSlicesMeanBetaCut->Print("output_mt/chi2pid_slices_pion_mean_beta_cut_forward.pdf"); // Add to the PDF
        delete gaussFit; // Clean up the fit function
    }
    canvasChi2PIDSlicesMeanBetaCut->Print("output_mt/chi2pid_slices_pion_mean_beta_cut_forward.pdf]"); // Close the PDF




    // Plot: Chi2PID vs. Momentum for Positive Pions (PID 211) with logarithmic z-axis
    TCanvas* canvasChi2PIDvsPPion = new TCanvas("canvasChi2PIDvsPPion", "Chi2PID vs. Momentum for Positive Pions (Forward)", 800, 600);
    TCanvas* canvasChi2PIDvsPPionMeanBetaCut = new TCanvas("canvasChi2PIDvsPPionMeanBetaCut", "Chi2PID vs. Momentum for Positive Pions (Forward, Beta > Mean(Pion,Kaon))", 800, 600);
    if (hChi2PIDvsPPID211.GetEntries() > 0 || hChi2PIDvsPPID211_MeanBetaCut.GetEntries() > 0) {

        canvasChi2PIDvsPPion->cd(); 
        hChi2PIDvsPPID211.SetStats(0);
        hChi2PIDvsPPID211.Draw("COLZ");
        gPad->SetLogz(); // Set logarithmic z-axis

        canvasChi2PIDvsPPionMeanBetaCut->cd(); 
        hChi2PIDvsPPID211_MeanBetaCut.SetStats(0);
        hChi2PIDvsPPID211_MeanBetaCut.Draw("COLZ");
        gPad->SetLogz();



        // Overlay the mean, mean + 3σ, and mean - 3σ as TGraphs
        TGraph* graphChi2PIDMean = new TGraph();
        TGraph* graphChi2PIDMeanPlus3Sigma = new TGraph();
        TGraph* graphChi2PIDMeanMinus3Sigma = new TGraph();
        int pointIndex = 0;
        for (int i = 0; i < nSlices; ++i) {
            if (chi2pidMeanPion[i] != 0.0 && chi2pidSigmaPion[i] != 0.0) { // Only plot slices with valid fits
                double pCenter = pMin + (i + 0.5) * pWidth; // Central momentum of the slice
                /* graphChi2PIDMean->SetPoint(pointIndex, chi2pidMeanPion[i], pCenter); 
                graphChi2PIDMeanPlus3Sigma->SetPoint(pointIndex, chi2pidMeanPlus3Sigma[i], pCenter);
                graphChi2PIDMeanMinus3Sigma->SetPoint(pointIndex, chi2pidMeanMinus3Sigma[i], pCenter);  */
                graphChi2PIDMean->SetPoint(pointIndex, pCenter, chi2pidMeanPion[i]); 
                graphChi2PIDMeanPlus3Sigma->SetPoint(pointIndex, pCenter, chi2pidMeanPlus3Sigma[i]);
                graphChi2PIDMeanMinus3Sigma->SetPoint(pointIndex,pCenter, chi2pidMeanMinus3Sigma[i]); 
                pointIndex++;
            }
        }

        // Set styles for the graphs
        graphChi2PIDMean->SetMarkerColor(kRed);
        graphChi2PIDMean->SetMarkerStyle(20); // Filled circle
        graphChi2PIDMean->SetMarkerSize(1.0);
    

        graphChi2PIDMeanPlus3Sigma->SetMarkerColor(kBlack);
        graphChi2PIDMeanPlus3Sigma->SetMarkerStyle(22); // Open triangle up
        graphChi2PIDMeanPlus3Sigma->SetMarkerSize(1.0);
        

        graphChi2PIDMeanMinus3Sigma->SetMarkerColor(kBlack);
        graphChi2PIDMeanMinus3Sigma->SetMarkerStyle(23); // Open triangle down
        graphChi2PIDMeanMinus3Sigma->SetMarkerSize(1.0);
        

        canvasChi2PIDvsPPion->cd();
        graphChi2PIDMean->Draw("P SAME");
        graphChi2PIDMeanPlus3Sigma->Draw("P SAME");
        graphChi2PIDMeanMinus3Sigma->Draw("P SAME"); 

        // Add a legend
        TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg->AddEntry(graphChi2PIDMean, "Mean Chi2PID", "p");
        leg->AddEntry(graphChi2PIDMeanPlus3Sigma, "Mean #pm 3#sigma", "p");
        leg->Draw();
    } else {
        cout << "Warning: hChi2PIDvsPPID211 (Pions, Forward) is empty, plot will be empty." << endl;
    }
    canvasChi2PIDvsPPion->Print("output_mt/chi2pid_vs_p_pion_forward.pdf");
    canvasChi2PIDvsPPionMeanBetaCut->Print("output_mt/chi2pid_vs_p_pion_forward_mean_beta_cut.pdf");




    // Plot: Chi2PID vs. Momentum for Positive Kaons (PID 321) with logarithmic z-axis
    TCanvas* canvasChi2PIDvsPKaon = new TCanvas("canvasChi2PIDvsPKaon", "Chi2PID vs. Momentum for Positive Kaons (Forward)", 800, 600);
    if (hChi2PIDvsPPID321.GetEntries() > 0) {
        hChi2PIDvsPPID321.SetStats(0);
        hChi2PIDvsPPID321.Draw("COLZ");
        gPad->SetLogz(); // Set logarithmic z-axis
    } else {
        cout << "Warning: hChi2PIDvsPPID321 (Kaons, Forward) is empty, plot will be empty." << endl;
    }
    canvasChi2PIDvsPKaon->Print("output_mt/chi2pid_vs_p_kaon_forward.pdf");




    // Plot: Chi2PID vs. Momentum for Protons (PID 2212) with logarithmic z-axis
    TCanvas* canvasChi2PIDvsPProton = new TCanvas("canvasChi2PIDvsPProton", "Chi2PID vs. Momentum for Protons (Forward)", 800, 600);
    if (hChi2PIDvsPPID2212.GetEntries() > 0) {
        hChi2PIDvsPPID2212.SetStats(0);
        hChi2PIDvsPPID2212.Draw("COLZ");
        gPad->SetLogz(); // Set logarithmic z-axis
    } else {
        cout << "Warning: hChi2PIDvsPPID2212 (Protons, Forward) is empty, plot will be empty." << endl;
    }
    canvasChi2PIDvsPProton->Print("output_mt/chi2pid_vs_p_proton_forward.pdf");

   


    // Plot: 1D Beta Histograms for Pions in Momentum Slices (one canvas per slice, saved in a multi-page PDF)
    TCanvas* canvasBetaPionSlice = new TCanvas("canvasBetaPionSlice", "Beta for Pions in Momentum Slice (Forward)", 800, 600);
    canvasBetaPionSlice->Print("output_mt/beta_slices_pion_forward.pdf["); // Open the PDF
    for (int i = 0; i < nSlices; ++i) {
        canvasBetaPionSlice->Clear();
        hBetaPionSlices[i]->SetStats(0);
        hBetaPionSlices[i]->SetLineColor(kRed);
        hBetaPionSlices[i]->Draw();

        // Fit a Gaussian to the histogram using the specified range
        double mean = 1.0; // Initial guess for the mean (beta is typically around 1)
        double sigma = 0.02; // Initial guess for sigma
        double betaMin = pionFitRanges[i][0];
        double betaMax = pionFitRanges[i][1];
        TF1* gaussFit = new TF1("gaussFit", "gaus", betaMin, betaMax);
        gaussFit->SetParameters(hBetaPionSlices[i]->GetMaximum(), mean, sigma);
        hBetaPionSlices[i]->Fit(gaussFit, "RQN"); // R: range, Q: quiet, N: no draw

        // Extract mean and sigma from the fit
        double fitMean = gaussFit->GetParameter(1);
        double fitSigma = gaussFit->GetParameter(2);

        // Draw the fit
        gaussFit->SetLineColor(kBlack);
        gaussFit->Draw("SAME");

        // Add a legend with the fit parameters
        TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg->AddEntry(hBetaPionSlices[i], "Pions", "l");
        char fitText[100];
        sprintf(fitText, "Mean = %.4f", fitMean);
        leg->AddEntry((TObject*)0, fitText, "");
        sprintf(fitText, "Sigma = %.4f", fitSigma);
        leg->AddEntry((TObject*)0, fitText, "");
        leg->Draw();

        if (hBetaPionSlices[i]->GetEntries() == 0) {
            cout << "Warning: Beta histogram for pions in momentum slice " << (pMin + i * pWidth) << " to " << (pMin + (i + 1) * pWidth) << " GeV is empty." << endl;
        }

        canvasBetaPionSlice->Print("output_mt/beta_slices_pion_forward.pdf"); // Add to the PDF
        delete gaussFit; // Clean up the fit function
    }
    canvasBetaPionSlice->Print("output_mt/beta_slices_pion_forward.pdf]"); // Close the PDF

       

    // Plot: 1D Beta Histograms for Kaons in Momentum Slices (one canvas per slice, saved in a multi-page PDF)
    TCanvas* canvasBetaKaonSlice = new TCanvas("canvasBetaKaonSlice", "Beta for Kaons in Momentum Slice (Forward)", 800, 600);
    canvasBetaKaonSlice->Print("output_mt/beta_slices_kaon_forward.pdf["); // Open the PDF
    for (int i = 0; i < nSlices; ++i) {
        canvasBetaKaonSlice->Clear();
        hBetaKaonSlices[i]->SetStats(0);
        hBetaKaonSlices[i]->SetLineColor(kBlue);
        hBetaKaonSlices[i]->Draw();

        // Fit a Gaussian to the histogram using the specified range
        double mean = 1.0; // Initial guess for the mean
        double sigma = 0.02; // Initial guess for sigma
        double betaMin = kaonFitRanges[i][0];
        double betaMax = kaonFitRanges[i][1];
        TF1* gaussFit = new TF1("gaussFit", "gaus", betaMin, betaMax);
        gaussFit->SetParameters(hBetaKaonSlices[i]->GetMaximum(), mean, sigma);
        hBetaKaonSlices[i]->Fit(gaussFit, "RQN"); // R: range, Q: quiet, N: no draw

        // Extract mean and sigma from the fit
        double fitMean = gaussFit->GetParameter(1);
        double fitSigma = gaussFit->GetParameter(2);

        // Draw the fit
        gaussFit->SetLineColor(kBlack);
        gaussFit->Draw("SAME");

        // Add a legend with the fit parameters
        TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg->AddEntry(hBetaKaonSlices[i], "Kaons", "l");
        char fitText[100];
        sprintf(fitText, "Mean = %.4f", fitMean);
        leg->AddEntry((TObject*)0, fitText, "");
        sprintf(fitText, "Sigma = %.4f", fitSigma);
        leg->AddEntry((TObject*)0, fitText, "");
        leg->Draw();

        if (hBetaKaonSlices[i]->GetEntries() == 0) {
            cout << "Warning: Beta histogram for kaons in momentum slice " << (pMin + i * pWidth) << " to " << (pMin + (i + 1) * pWidth) << " GeV is empty." << endl;
        }

        canvasBetaKaonSlice->Print("output_mt/beta_slices_kaon_forward.pdf"); // Add to the PDF
        delete gaussFit; // Clean up the fit function
    }
    canvasBetaKaonSlice->Print("output_mt/beta_slices_kaon_forward.pdf]"); // Close the PDF

    // Vectors to store mean, sigma, and momentum for the new plots
    vector<double> momentumValues(nSlices);
    vector<double> meanPion(nSlices), sigmaPion(nSlices);
    vector<double> meanKaon(nSlices), sigmaKaon(nSlices);
    vector<double> diffMeans(nSlices); // Difference of means (pion - kaon)
    vector<double> resolutionPion(nSlices), resolutionKaon(nSlices);

    // Plot: 1D Beta Histograms for Pions and Kaons (Combined) in Momentum Slices (one canvas per slice, saved in a multi-page PDF)
    TCanvas* canvasBetaPionKaonSlice = new TCanvas("canvasBetaPionKaonSlice", "Beta for Pions and Kaons in Momentum Slice (Forward)", 800, 600);
    canvasBetaPionKaonSlice->Print("output_mt/beta_slices_pion_kaon_forward.pdf["); // Open the PDF
    for (int i = 0; i < nSlices; ++i) {
        // Store the central momentum of the slice
        momentumValues[i] = pMin + (i + 0.5) * pWidth;

        canvasBetaPionKaonSlice->Clear();
        // Normalize histograms to maximum for better comparison
        if (hBetaPionSlices[i]->GetMaximum() > 0) hBetaPionSlices[i]->Scale(1.0 / hBetaPionSlices[i]->GetMaximum());
        if (hBetaKaonSlices[i]->GetMaximum() > 0) hBetaKaonSlices[i]->Scale(1.0 / hBetaKaonSlices[i]->GetMaximum());
        // Set styles
        hBetaPionSlices[i]->SetStats(0);
        hBetaPionSlices[i]->SetLineColor(kRed);
        hBetaPionSlices[i]->SetFillColorAlpha(kRed, 0.3); // Semi-transparent fill
        hBetaKaonSlices[i]->SetStats(0);
        hBetaKaonSlices[i]->SetLineColor(kBlue);
        hBetaKaonSlices[i]->SetFillColorAlpha(kBlue, 0.3); // Semi-transparent fill
        // Draw the histograms
        hBetaPionSlices[i]->Draw("HIST");
        hBetaKaonSlices[i]->Draw("HIST SAME");

        // Fit a Gaussian to the pion histogram using the specified range
        double meanPionVal = 1.0;
        double sigmaPionVal = 0.02;
        double betaMinPion = pionFitRanges[i][0];
        double betaMaxPion = pionFitRanges[i][1];
        TF1* gaussFitPion = new TF1("gaussFitPion", "gaus", betaMinPion, betaMaxPion);
        gaussFitPion->SetParameters(hBetaPionSlices[i]->GetMaximum(), meanPionVal, sigmaPionVal);
        if (hBetaPionSlices[i]->GetEntries() > 10) { // Ensure enough entries for a meaningful fit
            hBetaPionSlices[i]->Fit(gaussFitPion, "RQN"); // R: range, Q: quiet, N: no draw
            meanPion[i] = gaussFitPion->GetParameter(1);
            sigmaPion[i] = gaussFitPion->GetParameter(2);
        } else {
            meanPion[i] = 0.0;
            sigmaPion[i] = 0.0;
        }

        // Fit a Gaussian to the kaon histogram using the specified range
        double meanKaonVal = 1.0;
        double sigmaKaonVal = 0.02;
        double betaMinKaon = kaonFitRanges[i][0];
        double betaMaxKaon = kaonFitRanges[i][1];
        TF1* gaussFitKaon = new TF1("gaussFitKaon", "gaus", betaMinKaon, betaMaxKaon);
        gaussFitKaon->SetParameters(hBetaKaonSlices[i]->GetMaximum(), meanKaonVal, sigmaKaonVal);
        if (hBetaKaonSlices[i]->GetEntries() > 10) { // Ensure enough entries for a meaningful fit
            hBetaKaonSlices[i]->Fit(gaussFitKaon, "RQN"); // R: range, Q: quiet, N: no draw
            meanKaon[i] = gaussFitKaon->GetParameter(1);
            sigmaKaon[i] = gaussFitKaon->GetParameter(2);
        } else {
            meanKaon[i] = 0.0;
            sigmaKaon[i] = 0.0;
        }

        // Draw the fits
        gaussFitPion->SetLineColor(kRed);
        gaussFitPion->SetLineStyle(2); // Dashed line for the fit
        gaussFitPion->Draw("SAME");
        gaussFitKaon->SetLineColor(kBlue);
        gaussFitKaon->SetLineStyle(2); // Dashed line for the fit
        gaussFitKaon->Draw("SAME");

        // Add a legend with the fit parameters
        TLegend* leg = new TLegend(0.7, 0.5, 0.9, 0.9); // Adjusted position to fit more entries
        leg->AddEntry(hBetaPionSlices[i], "Pions", "f");
        char fitText[100];
        sprintf(fitText, "Mean = %.4f", meanPion[i]);
        leg->AddEntry((TObject*)0, fitText, "");
        sprintf(fitText, "Sigma = %.4f", sigmaPion[i]);
        leg->AddEntry((TObject*)0, fitText, "");
        leg->AddEntry(hBetaKaonSlices[i], "Kaons", "f");
        sprintf(fitText, "Mean = %.4f", meanKaon[i]);
        leg->AddEntry((TObject*)0, fitText, "");
        sprintf(fitText, "Sigma = %.4f", sigmaKaon[i]);
        leg->AddEntry((TObject*)0, fitText, "");
        leg->Draw();

        // Calculate difference of means and resolution
        diffMeans[i] = meanPion[i] - meanKaon[i];
        resolutionPion[i] = sigmaPion[i]; // Resolution is just σ
        resolutionKaon[i] = sigmaKaon[i]; // Resolution is just σ

        canvasBetaPionKaonSlice->Print("output_mt/beta_slices_pion_kaon_forward.pdf"); // Add to the PDF
        delete gaussFitPion;
        delete gaussFitKaon; // Clean up the fit functions
    }
    canvasBetaPionKaonSlice->Print("output_mt/beta_slices_pion_kaon_forward.pdf]"); // Close the PDF

    // Plot: Difference of Means (Pion - Kaon) vs. Momentum
    TCanvas* canvasDiffMeans = new TCanvas("canvasDiffMeans", "Difference of Means (Pion - Kaon) vs. Momentum (Forward)", 800, 600);
    TGraph* graphDiffMeans = new TGraph();
    int pointIndex = 0;
    for (int i = 0; i < nSlices; ++i) {
        if (meanPion[i] != 0.0 && meanKaon[i] != 0.0) { // Only plot slices with valid data
            graphDiffMeans->SetPoint(pointIndex, momentumValues[i], diffMeans[i]);
            pointIndex++;
        }
    }
    graphDiffMeans->SetTitle("Difference of Means (Pion - Kaon) vs. Momentum (Forward);Momentum p (GeV);Mean Beta (Pion) - Mean Beta (Kaon)");
    graphDiffMeans->SetMarkerStyle(20);
    graphDiffMeans->SetMarkerColor(kBlue);
    graphDiffMeans->Draw("AP");
    canvasDiffMeans->Print("output_mt/diff_means_vs_p_forward.pdf");

    // Plot: Resolution (σ) of Pions vs. Momentum (Separate Canvas)
    TCanvas* canvasResolutionPion = new TCanvas("canvasResolutionPion", "Resolution of Pions vs. Momentum (Forward)", 800, 600);
    TGraph* graphResolutionPion = new TGraph();
    pointIndex = 0;
    for (int i = 0; i < nSlices; ++i) {
        if (resolutionPion[i] != 0.0) { // Only plot slices with valid data
            graphResolutionPion->SetPoint(pointIndex, momentumValues[i], resolutionPion[i]);
            pointIndex++;
        }
    }
    graphResolutionPion->SetTitle("Resolution of Pions vs. Momentum (Forward);Momentum p (GeV);Resolution (#sigma)");
    graphResolutionPion->SetMarkerStyle(20);
    graphResolutionPion->SetMarkerColor(kRed);
    graphResolutionPion->SetLineColor(kRed);
    graphResolutionPion->Draw("APL");
    canvasResolutionPion->Print("output_mt/resolution_vs_p_pion_forward.pdf");

    // Plot: Resolution (σ) of Kaons vs. Momentum (Separate Canvas)
    TCanvas* canvasResolutionKaon = new TCanvas("canvasResolutionKaon", "Resolution of Kaons vs. Momentum (Forward)", 800, 600);
    TGraph* graphResolutionKaon = new TGraph();
    pointIndex = 0;
    for (int i = 0; i < nSlices; ++i) {
        if (resolutionKaon[i] != 0.0) {
            graphResolutionKaon->SetPoint(pointIndex, momentumValues[i], resolutionKaon[i]);
            pointIndex++;
        }
    }
    graphResolutionKaon->SetTitle("Resolution of Kaons vs. Momentum (Forward);Momentum p (GeV);Resolution (#sigma)");
    graphResolutionKaon->SetMarkerStyle(21);
    graphResolutionKaon->SetMarkerColor(kBlue);
    graphResolutionKaon->SetLineColor(kBlue);
    graphResolutionKaon->Draw("APL");
    canvasResolutionKaon->Print("output_mt/resolution_vs_p_kaon_forward.pdf");

    // Plot: Combined Resolution (σ) of Pions and Kaons vs. Momentum with Difference of Means (All on Left y-axis, No Scaling)
    TCanvas* canvasResolutionCombined = new TCanvas("canvasResolutionCombined", "Resolution and Difference of Means vs. Momentum (Forward)", 800, 600);

    // Left y-axis: Resolution (σ) for Pions with logarithmic scale
    TGraph* graphResPion = new TGraph();
    pointIndex = 0;
    for (int i = 0; i < nSlices; ++i) {
        if (resolutionPion[i] != 0.0) {
            graphResPion->SetPoint(pointIndex, momentumValues[i], resolutionPion[i]);
            pointIndex++;
        }
    }
    graphResPion->SetTitle("Resolution and Difference of Means vs. Momentum (Forward);Momentum p (GeV);Resolution (#sigma) or Gap ");
    graphResPion->SetMarkerStyle(20);
    graphResPion->SetMarkerColor(kRed);
    graphResPion->SetLineColor(kRed);
    graphResPion->GetYaxis()->SetRangeUser(1e-3, 1); // Adjusted range to accommodate shifted gap
    graphResPion->Draw("APL");
    gPad->SetLogy(); // Set logarithmic scale for the left y-axis

    // Plot Kaon Resolution on the same left y-axis
    TGraph* graphResKaon = new TGraph();
    pointIndex = 0;
    for (int i = 0; i < nSlices; ++i) {
        if (resolutionKaon[i] != 0.0) {
            graphResKaon->SetPoint(pointIndex, momentumValues[i], resolutionKaon[i]);
            pointIndex++;
        }
    }
    graphResKaon->SetMarkerStyle(21);
    graphResKaon->SetMarkerColor(kBlue);
    graphResKaon->SetLineColor(kBlue);
    graphResKaon->Draw("PL SAME");

    // Plot the difference of means (shifted to be positive for logarithmic scale)
    TGraph* graphDiffMeansShifted = new TGraph();
    pointIndex = 0;
    for (int i = 0; i < nSlices; ++i) {
        if (meanPion[i] != 0.0 && meanKaon[i] != 0.0) {
            if (diffMeans[i] < 0) continue; 
            graphDiffMeansShifted->SetPoint(pointIndex, momentumValues[i], diffMeans[i]);
            pointIndex++;
        }
    }
    graphDiffMeansShifted->SetMarkerStyle(22);
    graphDiffMeansShifted->SetMarkerColor(kGreen+2);
    graphDiffMeansShifted->SetLineColor(kGreen+2);
    graphDiffMeansShifted->Draw("PL SAME");

    // Add a legend
    TLegend* legend = new TLegend(0.7, 0.6, 0.85, 0.85);
    legend->AddEntry(graphResPion, "Pion Resolution", "lp");
    legend->AddEntry(graphResKaon, "Kaon Resolution", "lp");
    legend->AddEntry(graphDiffMeansShifted, "Pion - Kaon Gap", "lp");
    legend->Draw();

    canvasResolutionCombined->Print("output_mt/resolution_and_gap_vs_p_forward.pdf");

    // Clean up
    delete canvasChi2PIDvsPPion;
    delete canvasChi2PIDvsPKaon;
    delete canvasChi2PIDvsPProton;
    delete canvasBetaVsP;
    delete canvasBetaPionSlice;
    delete canvasBetaKaonSlice;
    delete canvasBetaPionKaonSlice;
    delete canvasDiffMeans;
    delete canvasResolutionPion;
    delete canvasResolutionKaon;
    delete canvasResolutionCombined;
    delete canvasChi2PIDSlices;
    delete canvasChi2PIDSlicesKaon;
    delete canvasGapVsP; 
    delete canvasSigmaVsP ; 
    for (int i = 0; i < nSlices; ++i) {
        delete hBetaPionSlices[i];
        delete hBetaKaonSlices[i];
        delete hChi2PIDSlices[i];
    }
    inputFile->Close();
    delete inputFile;

    // Disable multi-threading (good practice)
    ROOT::DisableImplicitMT();

    cout << "Program completed successfully." << endl;
    return 0;
}