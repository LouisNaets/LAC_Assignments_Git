import numpy as np
from numpy.polynomial.polynomial import Polynomial
import matplotlib.pyplot as plt

plt.rcParams.update({'axes.labelsize': 14, 'xtick.labelsize': 12, 'ytick.labelsize': 12, 'legend.fontsize': 10, 'axes.titlesize': 15})

theta = [0.69990, 1.33032, 1.87337, 2.34564, 2.78268, 3.18625, 3.55478, 6.41687, 8.55115,10.35828,11.97465,13.46319,14.85867,16.18002,17.44155,18.65245,19.81880,20.94841,22.04174,23.10374]

dQ_dtheta = [-1250.21385,-1294.80527,-1336.27132,-1391.04065,-1441.85613,-1486.51244,-1537.42141,-1817.22625,-2099.38349,-2355.84289,-2573.84846,-2803.55880,-3034.19836,-3248.08553,-3467.85105,-3684.03512,-3906.81836,-4126.44262,-4349.84614,-4582.75853]

# Perform quadratic fit
coefficients = np.polyfit(theta, dQ_dtheta, 2)  # 2 means quadratic

# Extract the coefficients
KK2, KK1, intercept = coefficients  # From highest degree to lowest (x^2, x, constant)

print('KK1:' + str(round(KK1,2)))
print('KK2:' + str(round(KK2,2)))
print('dQ_dtheta_0:' + str(round(intercept, 2)))

# Generate a range of theta values for the fit line
theta_fit = np.linspace(min(theta), max(theta), 500)

# Calculate the fitted dq_dtheta values based on the quadratic fit
dQ_dtheta_fit = np.polyval(coefficients, theta_fit)

# Create the plot
plt.figure(0, figsize=(10,6), dpi=500)
plt.scatter(theta, dQ_dtheta, color='tab:orange', label='HAWC2S Data', zorder=5, marker='x')
plt.plot(theta_fit, dQ_dtheta_fit, color='tab:blue', label='Quadratic Fit', zorder=4)

# Add labels and title
plt.xlabel('θ (x)', fontsize=12)
plt.ylabel('dQ/dθ (y)', fontsize=12)
plt.title('Quadratic Fit of dQ/dθ vs θ', fontsize=14)
plt.legend()

# Formula text for the plot
formula_text = f"dQ/dθ = {round(KK2, 2)}θ² + {round(KK1, 2)}θ + {round(intercept, 2)}"
# Add the formula to the plot
plt.text(0.05, 0.25, formula_text, fontsize=12, transform=plt.gca().transAxes, verticalalignment='top')

# Show the plot
plt.grid(True)

plt.savefig('A3 Controller design/Assignment/Figures/Part1_fit.svg', format='svg')
plt.savefig('A3 Controller design/Assignment/Figures/Part1_fit.png', format='png')