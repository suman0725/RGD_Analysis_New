import ROOT
import numpy as np
from ROOT import TFile, TH1F, TCanvas, TLegend, TF1, TLine
import time
import os

# Start timer
start_time = time.time()

# Enable multithreading
ROOT.EnableImplicitMT()
print("Multi-threading enabled for this fit.")

# Open the ROOT file and create RDataFrame
print("Opening file...")
file = TFile("/w/hallb-scshelf2102/clas12/suman/RGD_Analysis/PID/charge_particles_custompid/Misidentification/Skim/pkptreeCxC_9_test_modified.root", "READ")
if not file or file.IsZombie():
    print("Error: Cannot open file pkptreeCxC_9_test_modified.root")
    exit(1)
df_all = ROOT.RDataFrame("EB_all_pion_assumed", file)

# Focus on 4.0-4.3 GeV/c bin
pLow, pHigh = 4.0, 4.3
print("Filtering data...")
filtered_df = df_all.Filter(f"p >= {pLow} && p < {pHigh}")

# Create histogram for beta
print("Creating histogram...")
modelBeta = ROOT.RDF.TH1DModel("beta_fit", "p: [4.0-4.3) GeV/c; #beta; Counts", 40, 0.95, 1.03)
histo_beta = filtered_df.Histo1D(modelBeta, "beta")

# Check histogram content
print(f"Histogram entries: {histo_beta.GetEntries()}")
if histo_beta.GetEntries() < 100:
    print(f"Error: Insufficient statistics ({histo_beta.GetEntries()}) for fitting.")
    exit(1)

# Initialize canvas
canvas = TCanvas("canvas", "Beta Fit with Signals", 1200, 800)
canvas.SetGridx(1)
canvas.SetGridy(1)

# Draw histogram as black points
histo_beta.SetMarkerStyle(20)
histo_beta.SetMarkerSize(1.2)
histo_beta.SetMarkerColor(ROOT.kBlack)
histo_beta.Draw("P")
canvas.Update()

# Define fit function with three Gaussians + quadratic + power background
fitFunc = TF1("fitFunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[4])/[5])^2) + [6]*exp(-0.5*((x-[7])/[8])^2) + [9] + [10]*x + [11]*x*x + [12]/(x^[13])", 0.95, 1.03)

# Set parameter names individually
fitFunc.SetParName(0, "Pion_Amp")
fitFunc.SetParName(1, "Pion_Mean")
fitFunc.SetParName(2, "Pion_Sigma")
fitFunc.SetParName(3, "Kaon_Amp")
fitFunc.SetParName(4, "Kaon_Mean")
fitFunc.SetParName(5, "Kaon_Sigma")
fitFunc.SetParName(6, "Proton_Amp")
fitFunc.SetParName(7, "Proton_Mean")
fitFunc.SetParName(8, "Proton_Sigma")
fitFunc.SetParName(9, "Bg_Const")
fitFunc.SetParName(10, "Bg_Linear")
fitFunc.SetParName(11, "Bg_Quad")
fitFunc.SetParName(12, "Bg_Pow_Amp")
fitFunc.SetParName(13, "Bg_Pow_Exp")

# Theoretical beta values
p_mid = (pLow + pHigh) / 2
beta_pi = p_mid / np.sqrt(p_mid * p_mid + 0.1396 * 0.1396)  # ~0.997
beta_K = p_mid / np.sqrt(p_mid * p_mid + 0.4937 * 0.4937)   # ~0.989
beta_p = p_mid / np.sqrt(p_mid * p_mid + 0.9383 * 0.9383)   # ~0.980

# Find maximum amplitudes
bin_pi = histo_beta.FindBin(beta_pi - 0.002)
bin_pi_max = bin_pi
max_amp_pi = 0
for i in range(bin_pi, histo_beta.FindBin(beta_pi + 0.002) + 1):
    if histo_beta.GetBinContent(i) > max_amp_pi:
        max_amp_pi = histo_beta.GetBinContent(i)
        bin_pi_max = i
amp_pi = max_amp_pi
max_beta_pi = histo_beta.GetBinCenter(bin_pi_max)

bin_K = histo_beta.FindBin(beta_K - 0.002)
bin_K_max = bin_K
max_amp_K = 0
for i in range(bin_K, histo_beta.FindBin(beta_K + 0.002) + 1):
    if histo_beta.GetBinContent(i) > max_amp_K:
        max_amp_K = histo_beta.GetBinContent(i)
        bin_K_max = i
amp_K = max_amp_K
max_beta_K = histo_beta.GetBinCenter(bin_K_max)

bin_p = histo_beta.FindBin(beta_p - 0.002)
bin_p_max = bin_p
max_amp_p = 0
for i in range(bin_p, histo_beta.FindBin(beta_p + 0.002) + 1):
    if histo_beta.GetBinContent(i) > max_amp_p:
        max_amp_p = histo_beta.GetBinContent(i)
        bin_p_max = i
amp_p = max_amp_p
max_beta_p = histo_beta.GetBinCenter(bin_p_max)

print(f"Pion amplitude at beta = {max_beta_pi}: {amp_pi}")
print(f"Kaon amplitude at beta = {max_beta_K}: {amp_K}")
print(f"Proton amplitude at beta = {max_beta_p}: {amp_p}")

# Set initial parameters and limits
fitFunc.SetParameter(0, amp_pi)
fitFunc.SetParLimits(0, amp_pi * 0.5, amp_pi * 1.5)
fitFunc.SetParameter(1, beta_pi)
fitFunc.SetParLimits(1, beta_pi - 0.002, beta_pi + 0.002)
fitFunc.SetParameter(2, 0.007)
fitFunc.SetParLimits(2, 0.001, 0.015)

fitFunc.SetParameter(3, amp_K)
fitFunc.SetParLimits(3, amp_K * 0.5, amp_K * 1.5)
fitFunc.SetParameter(4, beta_K)
fitFunc.SetParLimits(4, beta_K - 0.002, beta_K + 0.002)
fitFunc.SetParameter(5, 0.007)
fitFunc.SetParLimits(5, 0.001, 0.015)

fitFunc.SetParameter(6, amp_p)
fitFunc.SetParLimits(6, amp_p * 0.5, amp_p * 1.5)
fitFunc.SetParameter(7, beta_p)
fitFunc.SetParLimits(7, beta_p - 0.002, beta_p + 0.002)
fitFunc.SetParameter(8, 0.007)
fitFunc.SetParLimits(8, 0.001, 0.015)

fitFunc.SetParameter(9, 1000.0)
fitFunc.SetParLimits(9, 0.0, 5000.0)
fitFunc.SetParameter(10, 0.0)
fitFunc.SetParLimits(10, -500.0, 500.0)
fitFunc.SetParameter(11, 0.0)
fitFunc.SetParLimits(11, -200.0, 200.0)
fitFunc.SetParameter(12, 3000.0)
fitFunc.SetParLimits(12, 1000.0, 5000.0)
fitFunc.SetParameter(13, 1.5)
fitFunc.SetParLimits(13, 1.0, 3.0)

# Multi-step fit
print("Fitting background...")
bgFitTemp = TF1("bgFitTemp", "[0] + [1]*x + [2]*x*x + [3]/(x^[4])", 0.95, 1.03)
bgFitTemp.SetParameters(1000.0, 0.0, 0.0, 3000.0, 1.5)
bgFitTemp.SetParLimits(0, 0.0, 5000.0)
bgFitTemp.SetParLimits(1, -500.0, 500.0)
bgFitTemp.SetParLimits(2, -200.0, 200.0)
bgFitTemp.SetParLimits(3, 1000.0, 5000.0)
bgFitTemp.SetParLimits(4, 1.0, 3.0)
histo_beta.Fit(bgFitTemp, "R0", "", 0.95, 1.03)  # Full range for background
fitFunc.SetParameter(9, bgFitTemp.GetParameter(0))
fitFunc.SetParameter(10, bgFitTemp.GetParameter(1))
fitFunc.SetParameter(11, bgFitTemp.GetParameter(2))
fitFunc.SetParameter(12, bgFitTemp.GetParameter(3))
fitFunc.SetParameter(13, bgFitTemp.GetParameter(4))

print("Fitting full model...")
histo_beta.Fit(fitFunc, "RME")  # Simplified to 'RME' for speed

# Draw components
fitFunc.SetLineColor(ROOT.kBlue)
fitFunc.SetLineStyle(1)
fitFunc.Draw("SAME")

pionFit = TF1("pionFit", "[0]*exp(-0.5*((x-[1])/[2])^2)", 0.95, 1.03)
pionFit.SetParameters(fitFunc.GetParameter(0), fitFunc.GetParameter(1), fitFunc.GetParameter(2))
pionFit.SetLineColor(ROOT.kBlue)
pionFit.SetLineStyle(2)
pionFit.Draw("SAME")

kaonFit = TF1("kaonFit", "[0]*exp(-0.5*((x-[1])/[2])^2)", 0.95, 1.03)
kaonFit.SetParameters(fitFunc.GetParameter(3), fitFunc.GetParameter(4), fitFunc.GetParameter(5))
kaonFit.SetLineColor(ROOT.kGreen)
kaonFit.SetLineStyle(2)
kaonFit.Draw("SAME")

protonFit = TF1("protonFit", "[0]*exp(-0.5*((x-[1])/[2])^2)", 0.95, 1.03)
protonFit.SetParameters(fitFunc.GetParameter(6), fitFunc.GetParameter(7), fitFunc.GetParameter(8))
protonFit.SetLineColor(ROOT.kRed)
protonFit.SetLineStyle(2)
protonFit.Draw("SAME")

bgFit = TF1("bgFit", "[0] + [1]*x + [2]*x*x + [3]/(x^[4])", 0.95, 1.03)
bgFit.SetParameters(fitFunc.GetParameter(9), fitFunc.GetParameter(10), fitFunc.GetParameter(11), fitFunc.GetParameter(12), fitFunc.GetParameter(13))
bgFit.SetLineColor(ROOT.kGray)
bgFit.SetLineStyle(2)
bgFit.Draw("SAME")

# Add theoretical beta lines
line_pi = TLine(beta_pi, 0, beta_pi, histo_beta.GetMaximum())
line_pi.SetLineColor(ROOT.kMagenta)
line_pi.SetLineStyle(2)
line_pi.Draw()

line_K = TLine(beta_K, 0, beta_K, histo_beta.GetMaximum())
line_K.SetLineColor(ROOT.kMagenta)
line_K.SetLineStyle(2)
line_K.Draw()

line_p = TLine(beta_p, 0, beta_p, histo_beta.GetMaximum())
line_p.SetLineColor(ROOT.kMagenta)
line_p.SetLineStyle(2)
line_p.Draw()

# Add legend
leg = TLegend(0.7, 0.6, 0.9, 0.8)
leg.AddEntry(histo_beta.GetPtr(), f"Data ({int(histo_beta.GetEntries())})", "p")
leg.AddEntry(fitFunc, "Total Fit", "l")
leg.AddEntry(pionFit, "Pion Signal", "l")
leg.AddEntry(kaonFit, "Kaon Signal", "l")
leg.AddEntry(protonFit, "Proton Signal", "l")
leg.AddEntry(bgFit, "Background", "l")
leg.AddEntry(line_pi, "Theoretical #beta", "l")
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.03)
leg.Draw()

# Ensure output directory exists
os.makedirs("output/fit_plots", exist_ok=True)

# Save to PDF and PNG
canvas.Update()
canvas.Print("output/fit_plots/beta_fit_4.0-4.3.pdf")
canvas.Print("output/fit_plots/beta_fit_4.0-4.3.png")

# Clean up
del leg, line_pi, line_K, line_p, pionFit, kaonFit, protonFit, bgFit, fitFunc, canvas
file.Close()

# End timer
end_time = time.time()
print(f"Total program time: {(end_time - start_time):.2f} seconds")
print("Program completed successfully.")