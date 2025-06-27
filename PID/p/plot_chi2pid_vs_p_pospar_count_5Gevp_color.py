import ROOT
import numpy as np

# Open the ROOT file
file = ROOT.TFile("charged_particles.root", "READ")

# Access the TTree
tree = file.Get("charged_particle")

# Create a 2D histogram for chi2pid vs p (combined)
h_chi2pid_vs_p_all_positive = ROOT.TH2F("h_chi2pid_vs_p_all_positive", 
                                        "Chi2Pid vs Momentum (All Positive Particles);Momentum (GeV/c);Chi2Pid", 
                                        100, 0, 10, 100, -10, 10)

# Create a 2D histogram for chi2pid vs p for positive pions and positive kaons
h_chi2pid_vs_p_pion_kaon = ROOT.TH2F("h_chi2pid_vs_p_pion_kaon", 
                                     "Chi2Pid vs Momentum (Positive Pions and Kaons);Momentum (GeV/c);Chi2Pid", 
                                     100, 0, 10, 100, -10, 10)

# Counters for positive particles, positive pions, and positive kaons
positive_particles_count = 0
positive_pion_count = 0
positive_kaon_count = 0

# Loop over the tree entries
for entry in tree:
    px = entry.px
    py = entry.py
    pz = entry.pz
    chi2pid = entry.chi2pid
    pid = entry.pid  # Access the particle ID
    charge = entry.charge  # Access the particle charge

    # Calculate momentum magnitude (p)
    p = np.sqrt(px**2 + py**2 + pz**2)

    # Fill histogram for all positive particles based on charge
    if charge > 0:  # Only for positive particles
        positive_particles_count += 1
        h_chi2pid_vs_p_all_positive.Fill(p, chi2pid)
    
    # Fill histogram for positive pions and positive kaons based on pid
    if charge > 0 and (pid == 211 or pid == 321):  # Positive pions (pid=211) and positive kaons (pid=321)
        if pid == 211:
            positive_pion_count += 1
        elif pid == 321:
            positive_kaon_count += 1
        h_chi2pid_vs_p_pion_kaon.Fill(p, chi2pid)

# Create a canvas
canvas = ROOT.TCanvas("canvas", "Chi2Pid vs Momentum", 800, 600)

# Draw the first histogram for all positive particles
h_chi2pid_vs_p_all_positive.SetLineColor(ROOT.kBlue)  # Color for all positive particles
h_chi2pid_vs_p_all_positive.Draw("COLZ")

# Draw the second histogram for positive pions and positive kaons
h_chi2pid_vs_p_pion_kaon.SetLineColor(ROOT.kRed)  # Color for positive pions and positive kaons
h_chi2pid_vs_p_pion_kaon.Draw("COLZ SAME")  # Overlay on top of the first histogram

# Add a legend to distinguish between the histograms
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(h_chi2pid_vs_p_all_positive, "All Positive Particles (Charge > 0)", "l")
legend.AddEntry(h_chi2pid_vs_p_pion_kaon, "Positive Pions (pid=211) and Kaons (pid=321)", "l")
legend.Draw()

# Add the counts to the plot
latex = ROOT.TLatex()
latex.SetTextSize(0.03)
latex.SetNDC()  # Use normalized device coordinates
latex.DrawLatex(0.15, 0.85, f"Positive Particles: {positive_particles_count}")
latex.DrawLatex(0.15, 0.80, f"Positive Pions (pid=211): {positive_pion_count}")
latex.DrawLatex(0.15, 0.75, f"Positive Kaons (pid=321): {positive_kaon_count}")

# Save the plot
canvas.SaveAs("chi2pid_vs_p_combined.png")

# Close the ROOT file
file.Close()
