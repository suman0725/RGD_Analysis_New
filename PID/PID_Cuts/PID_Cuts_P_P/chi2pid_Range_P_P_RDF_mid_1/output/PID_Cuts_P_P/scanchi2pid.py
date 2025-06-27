import numpy as np
import matplotlib.pyplot as plt

# Fit function parameters from the plot
offset = -0.20
amp1 = 3.00
decay1 = 4.51
amp2 = 37.26
decay2 = 0.87
threshold = 3.0  # Constant value for p <= 2.75

# Define the fit function
def chi2pid(p):
    if p <= 2.75:
        return threshold
    else:
        return offset + amp1 * np.exp(-p / decay1) + amp2 * np.exp(-p / decay2)

# Momentum range to scan (same as your data)
p_values = np.arange(0.2, 8.1, 0.1)  # From 0.2 to 8.0 with 0.1 steps

# Calculate chi2pid for each momentum
chi2pid_values = [chi2pid(p) for p in p_values]

# Print specific value at p = 4 GeV/c
p_target = 4.0
chi2pid_at_4 = chi2pid(p_target)
print(f"Chi2pid at {p_target} GeV/c = {chi2pid_at_4:.2f}")

# Scan and print all values
for p, chi2 in zip(p_values, chi2pid_values):
    print(f"p = {p:.1f} GeV/c: Chi2pid = {chi2:.2f}")

# Optional: Plot the result
plt.plot(p_values, chi2pid_values, 'r-', label='Fit Function')
plt.scatter(p_values, chi2pid_values, color='blue', s=10, label='Calculated Points')
plt.axvline(x=2.75, color='gray', linestyle='--', label='Transition at 2.75 GeV/c')
plt.xlabel('Momentum (GeV/c)')
plt.ylabel('Cut Position (chi2pid)')
plt.title('Chi2pid vs Momentum')
plt.legend()
plt.grid(True)
plt.show()