import numpy as np
import matplotlib.pyplot as plt

# Bin data
bins = [
    (1.0, 1.3, -0.00625907, -0.00625907),
    (1.3, 1.6, -0.0364921, -0.0364921),
    (1.6, 1.9, -2.79599, 2.69247),
    (1.9, 2.2, -2.80432, 2.68414),
    (2.2, 2.5, -2.81214, 2.67632),
    (2.5, 2.8, -2.81467, 2.67379),
    (2.8, 3.1, -2.82121, 1.7525),
    (3.1, 3.4, -2.82711, 1.74661),
    (3.4, 3.7, -2.82313, 1.75059),
    (3.7, 4.0, -2.81815, 1.75557),
    (4.0, 4.3, -2.82645, 1.74726),
    (4.3, 4.6, -2.85151, 1.7222),
    (4.6, 4.9, -2.87243, 1.70129)
]

# Extract bin centers and chi2 values (from 3rd bin onward)
p = [(low + high) / 2 for (low, high, _, _) in bins[2:]]
chi2_neg = [neg for (_, _, neg, _) in bins[2:]]
chi2_pos = [pos for (_, _, _, pos) in bins[2:]]

# Convert to numpy arrays
p = np.array(p)
chi2_neg = np.array(chi2_neg)
chi2_pos = np.array(chi2_pos)

# Fit chi2_neg with a quadratic polynomial
coeffs_neg = np.polyfit(p, chi2_neg, 2)
poly_neg = np.poly1d(coeffs_neg)

# Fitted range
x_fit = np.linspace(1.5, 5, 500)

# Positive chi2pid exponential fit
def chi2_pos_fit(p_vals, C=1.0):
    p_vals = np.array(p_vals)
    y = np.zeros_like(p_vals)

    # Region 1: constant for p < 2.8
    mask1 = p_vals < 2.8
    y[mask1] = C * 2.68  # Based on data average

    # Region 2: exponential decay for p >= 2.8
    mask2 = p_vals >= 2.8
    y[mask2] = C * (1.66 + 1.10 * np.exp(-p_vals[mask2] / 0.8) +
                    0.14 * np.exp(-p_vals[mask2] / 3.2))

    return y

# Generate fits
y_neg_fit = poly_neg(x_fit)
y_pos_fit = chi2_pos_fit(x_fit, C=1.0)

# Plotting
plt.figure(figsize=(8, 5))
plt.scatter(p, chi2_neg, color='blue', marker='o')  # negative
plt.scatter(p, chi2_pos, color='red', marker='x')   # positive
plt.plot(x_fit, y_neg_fit, 'b--')                   # negative fit
plt.plot(x_fit, y_pos_fit, 'r--')                   # positive fit

plt.xlabel('Momentum (GeV/c)')
plt.ylabel('chi²pid')
plt.title('chi²pid vs. Momentum')
plt.grid(True)
plt.xlim(1.5, 5)
plt.tight_layout()
plt.savefig("chi2pid_piecewise_fit_nolegend.png", dpi=300)
plt.show()
