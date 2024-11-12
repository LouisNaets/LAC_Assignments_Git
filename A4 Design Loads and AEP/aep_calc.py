import numpy as np
import h5py
import matplotlib.pyplot as plt

# File paths
hdf5_file_path = './A4 Design Loads and AEP/stats_files/dtu_10mw_turb_stats.hdf5'
debug_file_path = './A4 Design Loads and AEP/In class material/debug_dtu10mw_designloads_aep/debug_aep_1A.txt'

# Turbine and wind class parameters
GENEFF = 0.94  # Generator efficiency
HOURS_IN_YEAR = 8760  # Total hours in a year

# Weibull parameters for wind class 1A
ubar = 10  # Mean wind speed [m/s]
TI = 0.18  # Turbulence intensity
k = 2      # Shape parameter
c = 1.13 * ubar  # Scale parameter
v_bins = np.arange(5, 25)  # Wind speed bins

def extract_floats(line):
    """Extract float numbers from a line of text."""
    # Remove non-numeric characters like brackets and split the line
    clean_line = line.replace('[', '').replace(']', '').replace(',', '').strip()
    return np.array([float(x) for x in clean_line.split() if x.replace('.', '', 1).isdigit()])

with open(debug_file_path, 'r') as f:
    lines = f.readlines()
    bin_edges, bin_probabilities, power_output_bins = [],[],[]
    # Extract data using the updated function
    for i in range(9, 12):
        bin_edges.extend(extract_floats(lines[i]))
      # Extract bin probabilities spanning multiple lines (lines 11 to 14)
    for i in range(12, 17):
        bin_probabilities.extend(extract_floats(lines[i]))
    for i in range(18, 23):
        power_output_bins.extend(extract_floats(lines[i]))
     
# Convert lists to numpy arrays
bin_edges = np.array(bin_edges, dtype=float)
bin_probabilities = np.array(bin_probabilities, dtype=float)
power_output_bins = np.array(power_output_bins, dtype=float)

# Check if the data was extracted correctly
print("Bin edges:", bin_edges)
print("Bin probabilities:", bin_probabilities)
print('Weibull sum: ', sum(bin_probabilities))
print("Power output bins:", power_output_bins)

# Verify that the sum of probabilities is close to 1
# assert np.isclose(sum(bin_probabilities), 1), "Bin probabilities do not sum up to 1."

# Load power data from HDF5 file
def load_power_data(hdf5_file):
    with h5py.File(hdf5_file, 'r') as f:
        wsp = f['wsp'][:]  # Wind speeds
        el_power = f['pelec'][:]  # Electrical power output [W]
        return wsp, el_power

wsp, el_power = load_power_data(hdf5_file_path)

# Interpolate power output for the wind speed bins
power_curve = np.interp(v_bins, wsp, el_power) / 1e6  # Convert W to MW

# Calculate AEP
AEP = sum(power_curve * bin_probabilities) * HOURS_IN_YEAR  # [MWh]
print(f"AEP calculated: {AEP:.2f} GWh")

# Verify against the debug AEP value
reference_AEP = 49.267999785414275  # GWh from debug_aep_1A.txt
print(f"Reference AEP: {reference_AEP:.2f} GWh")
print(f"Difference: {abs(AEP - reference_AEP):.4f} GWh")

# Plotting the power curve and bin probabilities
plt.figure(figsize=(10, 6))
plt.bar(v_bins, bin_probabilities, width=0.8, alpha=0.7, label='Bin Probabilities')
plt.plot(v_bins, power_curve, 'r-o', label='Power Curve [MW]')
plt.xlabel('Wind Speed [m/s]')
plt.ylabel('Power [MW] / Probability')
plt.title('Power Curve and Wind Speed Distribution')
plt.legend()
plt.grid(True)
plt.show()
