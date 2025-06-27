import ROOT
import numpy as np

# Open the ROOT file
file = ROOT.TFile("charged_particles.root", "READ")

# Access the TTree
tree = file.Get("charged_particle")

# Create a 2D histogram for beta vs p
h_beta_vs_p = ROOT.TH2F("h_beta_vs_p", "Beta vs Momentum;Momentum (GeV/c);Beta", 
                        100, 0, 5, 100, 0, 1.2)  # Adjust the y-axis range if needed

# Loop over the tree entries
for entry in tree:
    charge = entry.charge
    if charge > 0:  # Include only positive particles
        px = entry.px
        py = entry.py
        pz = entry.pz
        beta = entry.beta  # Access beta directly

        # Calculate momentum magnitude (p)
        p = np.sqrt(px**2 + py**2 + pz**2)

        # Fill the histogram
        h_beta_vs_p.Fill(p, beta)

# Draw the histogram
canvas = ROOT.TCanvas("canvas", "Beta vs Momentum", 800, 600)
h_beta_vs_p.Draw("COLZ")

# Save the plot
canvas.SaveAs("beta_vs_p_posi_parti.png")

# Close the ROOT file
file.Close()
