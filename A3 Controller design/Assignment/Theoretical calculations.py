import numpy as np
from numpy.polynomial.polynomial import Polynomial
import matplotlib.pyplot as plt
from lacbox.io import load_ctrl_txt

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
plt.figure(0, figsize=(6,4), dpi=500)
plt.scatter(theta, dQ_dtheta, color='tab:orange', label='HAWC2S Data', zorder=5, marker='x')
plt.plot(theta_fit, dQ_dtheta_fit, color='tab:blue', label='Quadratic Fit', zorder=4)

# Add labels and title
plt.xlabel('θ [deg]', fontsize=12)
plt.ylabel('dQ/dθ [kNm/deg]', fontsize=12)
plt.title('Quadratic Fit of dQ/dθ vs θ', fontsize=14)
plt.legend()

# Formula text for the plot
formula_text = f"dQ/dθ = {round(KK2, 2)}θ² {round(KK1, 2)}θ {round(intercept, 2)}"
# Add the formula to the plot
plt.text(0.05, 0.25, formula_text, fontsize=12, transform=plt.gca().transAxes, verticalalignment='top')

# Show the plot
plt.grid(True)

plt.savefig('A3 Controller design/Assignment/Figures/Part1_fit.svg', format='svg')
plt.savefig('A3 Controller design/Assignment/Figures/Part1_fit.png', format='png')

print(intercept/KK1)
print(intercept/KK2)


'''Theoretical estimations'''
txt_dict = load_ctrl_txt('hawc_files/our_design/res_hawc2s/group7_3B_design_controller_tuning_ctrl_tuning.txt')

eta = 1
rho = 1.225
zeta = 0.7
omega = np.array([0.05, 0.06])*(2*np.pi)
I = txt_dict["Irotor_kg*m^2"]
print(I)
dQ_dOmega_0 = -1967.33604*1000
P_r = 10000000
omega_r = 8.627*(2*np.pi)/60 #from RPM to rad/s
dQ_dOmega_op = -P_r/(omega_r**2)
dQ_dtheta_0 = np.rad2deg(-1185.59951)*1000
a_0 = KK2
a_1 = KK1
R = 92.4
C_p_opt = 0.441
lambda_opt = 7.5

K_opt = (eta*rho*np.pi*R**5*C_p_opt)/(2*lambda_opt**3)
print('{:.3e}'.format(K_opt))

k_P_t = 2*eta*zeta*omega[0]*I
print('{:.3e}'.format(k_P_t))

k_I_t = eta*I*omega[0]**2
print('{:.3e}'.format(k_I_t))

k_P_p = (2*zeta*omega[1]*I - 1/eta * dQ_dOmega_op)/-dQ_dtheta_0
print('{:.3e}'.format(k_P_p))

k_I_p = (omega[1]**2 * I) / -dQ_dtheta_0
print('{:.3e}'.format(k_I_p))

print(dQ_dOmega_op)