import matplotlib.pyplot as plt

# Data
momentum_bins = ["1-1.3", "1.3-1.6", "1.6-1.9", "1.9-2.2", "2.2-2.5", "2.5-2.8", "2.8-3.1", "3.1-3.4", "3.4-3.7", "3.7-4"]
intersections = [0.955807, 0.971015, 0.979634, 0.984947, 0.988438, 0.99085, 0.992583, 0.993868, 0.994848, 0.995611]
contaminations = [0, 0.00999753, 0.0611711, 0.0701307, 0.159335, 0.264269, 0.36334, 1.19483, 1.77175, 7.08274]

# Plot
plt.figure(figsize=(8, 6))
plt.plot(momentum_bins, contaminations, 'b.', markersize=10)  # Blue points
plt.xlabel('Momentum (GeV/c)', fontsize=14, labelpad=10)  # Larger font size and padding
plt.ylabel('Contamination (%)', fontsize=14, labelpad=10)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('contamination_plot.png', dpi=300)
plt.close()