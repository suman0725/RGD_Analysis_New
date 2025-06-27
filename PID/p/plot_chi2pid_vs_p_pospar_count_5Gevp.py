import ROOT
import numpy as np

# Open the ROOT file
file = ROOT.TFile("charged_particles.root", "READ")

# Access the TTree
tree = file.Get("charged_particle")

# Create a 2D histogram for chi2pid vs p
h_chi2pid_vs_p = ROOT.TH2F("h_chi2pid_vs_p", "Chi2Pid vs Momentum;Momentum (GeV/c);Chi2Pid", 
                           100, 0, 5, 100, -5, 5)

# Counters for positive pions and kaons
positive_pion_count = 0
kaon_count = 0

# Loop over the tree entries
for entry in tree:
    charge = entry.charge
    if charge > 0:  # Include only positive particles

        px = entry.px
        py = entry.py
        pz = entry.pz
        chi2pid = entry.chi2pid
        pid = entry.pid  # Access the particle ID

        # Count positive pions and kaons
        if pid == 211:
            positive_pion_count += 1
        elif pid == 321:
            kaon_count += 1

        # Calculate momentum magnitude (p)
        p = np.sqrt(px**2 + py**2 + pz**2)

        # Fill the histogram
        h_chi2pid_vs_p.Fill(p, chi2pid)

# Create a canvas and draw the histogram
canvas = ROOT.TCanvas("canvas", "Chi2Pid vs Momentum", 800, 600)
h_chi2pid_vs_p.Draw("COLZ")

# Add the counts to the plot
latex = ROOT.TLatex()
latex.SetTextSize(0.03)
latex.SetNDC()  # Use normalized device coordinates

# Save the plot
canvas.SaveAs("chi2pid_vs_p_with_counts_5GeVp.png")

# Close the ROOT file
file.Close()
