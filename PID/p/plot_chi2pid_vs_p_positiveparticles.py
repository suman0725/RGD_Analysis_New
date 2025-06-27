import ROOT
import numpy as np

# Open the ROOT file
file = ROOT.TFile("charged_particles.root", "READ")

# Access the TTree
tree = file.Get("charged_particle")

# Create a 2D histogram for chi2pid vs p
h_chi2pid_vs_p = ROOT.TH2F("h_chi2pid_vs_p", "Chi2Pid vs Momentum;Momentum (GeV/c);Chi2Pid", 
                           100, 0, 10, 100, -10, 10)

# Loop over the tree entries
for entry in tree:
    charge = entry.charge
    if charge > 0:  # Include only positive particles
        px = entry.px
        py = entry.py
        pz = entry.pz
        chi2pid = entry.chi2pid

        # Calculate momentum magnitude (p)
        p = np.sqrt(px**2 + py**2 + pz**2)

        # Fill the histogram
        h_chi2pid_vs_p.Fill(p, chi2pid)

# Draw the histogram
canvas = ROOT.TCanvas("canvas", "Chi2Pid vs Momentum", 800, 600)
h_chi2pid_vs_p.Draw("COLZ")

# Save the plot
canvas.SaveAs("chi2pid_vs_p_posi_parti.png")

# Close the ROOT file
file.Close()
