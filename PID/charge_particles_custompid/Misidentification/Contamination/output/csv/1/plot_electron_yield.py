import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Data: momentum bins (lower edges) and yields
bin_edges = [1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 3.4, 3.7, 4.0, 4.3, 4.6, 4.9, 5.2]
yields = [0.000, 0.0000, 0.0013, 0.0239, 0.1700, 0.5987, 0.7892, 0.6832, 1.7735, 3.1190, 1.8090, 9.6823, 3.2725, 11.4314]

# Calculate bin centers
bin_centers = [(bin_edges[i] + bin_edges[i+1]) / 2 for i in range(len(bin_edges) - 1)]

# Define the power law with offset function: y = a * x^b + c
def power_law(x, a, b, c):
    return a * np.power(x, b) + c

# Fit the power law to the data
popt, pcov = curve_fit(power_law, bin_centers, yields, p0=[1.0, 2.0, 0.0])  # Initial guess: a=1, b=2, c=0
a, b, c = popt

# Generate points for the smooth fit curve
fit_x = np.linspace(min(bin_centers), max(bin_centers), 100)
fit_y = power_law(fit_x, a, b, c)

# Plot as data points and fit
plt.figure(figsize=(8, 6))
plt.plot(bin_centers, yields, 'bo', label='Contamination')  # Blue circles, no line
plt.plot(fit_x, fit_y, 'r-', label=f'Fit: {a:.2f}x^{b:.2f} + {c:.2f}')  # Red line for the fit
plt.xlabel('Momentum (GeV)')
plt.ylabel('Contamination(+k -> +pion)')
plt.title('Contamination vs Momentum (CuSn)')
plt.grid(True)
plt.legend()

# Save the plot
plt.savefig('contamination_yield.png')
plt.close()