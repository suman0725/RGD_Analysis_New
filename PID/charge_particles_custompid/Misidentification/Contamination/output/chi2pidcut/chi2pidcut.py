import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit

# Define the file path (relative to chi2pidcut directory)
file_path = '../csv/proton/chipid_cuts.csv'

# Debug: Print the absolute path to verify
absolute_path = os.path.abspath(file_path)
print(f"Attempting to read file at: {absolute_path}")

# Check if the file exists
if not os.path.exists(absolute_path):
    print(f"Error: File not found at {absolute_path}. Using provided data instead.")
    # Fallback to provided proton data
    data = """Momentum Bin (GeV/c),Pion Mean,Pion Sigma,Chi2 Min (Mean-3σ),Chi2 Max (Mean+3σ)
1-1.3,-0.025385,0.947292,-4.76185,4.71108
1.3-1.6,-0.0655028,0.922497,-4.67799,4.54698
1.6-1.9,-0.0701425,0.918273,-4.66151,4.52122
1.9-2.2,-0.0801434,0.896183,-4.56106,4.40077
2.2-2.5,-0.0889431,0.918253,-4.68021,4.50232
2.5-2.8,-0.0916486,0.91381,-4.6607,4.4774
2.8-3.1,-0.0946585,0.894425,-4.56678,4.37747
3.1-3.4,-0.114579,0.890817,-4.56866,4.3395
3.4-3.7,-0.0818719,0.84869,-4.32532,4.16158
3.7-4,-0.112288,0.877792,-4.50125,4.27667
4-4.3,-0.053756,0.864494,-4.37622,4.26871
4.3-4.6,-0.136295,0.782112,-4.04686,3.77427
4.6-4.9,-0.127355,0.622426,-3.23948,2.98477
4.9-5.2,-0.170755,0.576434,-3.05292,2.71141
5.2-5.5,-0.163832,0.554404,-2.93585,2.60819
5.5-5.8,-0.180414,0.57407,-3.05076,2.68994
5.8-6.1,-0.156212,0.504954,-2.68098,2.36856
6.1-6.4,-0.241216,0.586669,-3.17456,2.69213
6.4-6.7,-0.273223,0.472779,-2.63712,2.09067
6.7-7,-0.248791,0.520706,-2.85232,2.35474"""
    from io import StringIO
    df = pd.read_csv(StringIO(data))
else:
    # Read the CSV file
    df = pd.read_csv(file_path)

# Extract momentum bins and calculate midpoints
momentum_bins = df['Momentum Bin (GeV/c)'].str.split('-', expand=True).astype(float)
momentum_midpoints = (momentum_bins[0] + momentum_bins[1]) / 2

# Extract Chi2 Min and Chi2 Max
chi2_min = df['Chi2 Min (Mean-3σ)']
chi2_max = df['Chi2 Max (Mean+3σ)']

# Convert to numpy arrays
momentum_midpoints = momentum_midpoints.to_numpy()
chi2_min = chi2_min.to_numpy()
chi2_max = chi2_max.to_numpy()

# Define the double exponential + constant function
def double_exp_func(p, a, b, c, d, e):
    return a + b * np.exp(-p/c) + d * np.exp(-p/e)

# Fit the function to Chi2 Min
popt_min, _ = curve_fit(double_exp_func, momentum_midpoints, chi2_min, p0=[-4, 1, 1, 1, 5], maxfev=10000)
# Fit the function to Chi2 Max
popt_max, _ = curve_fit(double_exp_func, momentum_midpoints, chi2_max, p0=[4, -1, 1, -1, 5], maxfev=10000)

# Generate smooth x-values for plotting
x_smooth = np.linspace(min(momentum_midpoints), max(momentum_midpoints), 100)

# Calculate fitted curves
fit_min = double_exp_func(x_smooth, *popt_min)
fit_max = double_exp_func(x_smooth, *popt_max)

# Create the plot
plt.figure(figsize=(10, 6))

# Plot the data points as scatter
plt.scatter(momentum_midpoints, chi2_min, label='Chi2 Min (Mean - 3σ)', marker='o', color='blue')
plt.scatter(momentum_midpoints, chi2_max, label='Chi2 Max (Mean + 3σ)', marker='s', color='red')

# Plot the fitted curves
plt.plot(x_smooth, fit_min, label=f'Chi2 Min Fit: y = {popt_min[0]:.2f} + {popt_min[1]:.2f} * exp(-p/{popt_min[2]:.2f}) + {popt_min[3]:.2f} * exp(-p/{popt_min[4]:.2f})', color='blue')
plt.plot(x_smooth, fit_max, label=f'Chi2 Max Fit: y = {popt_max[0]:.2f} + {popt_max[1]:.2f} * exp(-p/{popt_max[2]:.2f}) + {popt_max[3]:.2f} * exp(-p/{popt_max[4]:.2f})', color='red')

# Add labels and title
plt.xlabel('Momentum (GeV/c)')
plt.ylabel('Chi2 PID')
plt.title('Chi2 PID vs Momentum (Proton) with Fitted Curves')
plt.legend()
plt.grid(True)

# Save the plot to a file
plt.savefig('chi2pid_double_exp_fit_plot.png')
print("Plot saved as chi2pid_double_exp_fit_plot.png")

# Show the plot (if running in an environment with a display)
plt.show()