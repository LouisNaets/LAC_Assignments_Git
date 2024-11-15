import numpy as np
from scipy.stats import weibull_min
import h5py


def SafetyFactorsTest():
# Calculate and print partial safety factors for each load channel
    load_channels = [
        {"name": "TbFA", "extreme": 215879.984375, "characteristic": 291437.97890625003, "design": 364297.47363281256},
        {"name": "TbSS", "extreme": 69760.18684895833, "characteristic": 94176.25224609375, "design": 117720.31530761719},
        {"name": "YbTilt", "extreme": 32585.416666666668, "characteristic": 43990.31250000001, "design": 54987.89062500001},
        {"name": "YbRoll", "extreme": 14167.382486979166, "characteristic": 19125.966357421876, "design": 23907.457946777344},
        {"name": "ShftTrs", "extreme": -12068.249674479166, "characteristic": 16292.137060546875, "design": 20365.171325683594},
        {"name": "OoPBRM", "extreme": -42829.711588541664, "characteristic": 57820.11064453125, "design": 72275.13830566406},
        {"name": "IPBRM", "extreme": 24160.418619791668, "characteristic": 32616.565136718753, "design": 40770.70642089844},
    ]

    for channel in load_channels:
        SF1 = channel["design"] / channel["characteristic"]
        SF2 = channel["characteristic"] / channel["extreme"]
        print(f"Load channel: {channel['name']} [kNm]")
        print(f"Extreme value: {channel['extreme']}")
        print(f"Characteristic value: {channel['characteristic']}")
        print(f"Design value: {channel['design']}")
        print(f"Partial safety factor 1: {SF1}\n")
        print(f"Partial safety factor 2: {SF2}\n")
    


def ProbabilityBinsTest(c, k):
    U_ave = 10
    bin_probabilities_verification = [0.06442809, 0.0709227, 0.07472189, 0.07591763, 0.0747455, 0.07155162,
                                      0.06675382, 0.06080169, 0.05413969, 0.04717648, 0.04026251, 0.03367663,
                                      0.02762128, 0.02222493, 0.01755019, 0.01360521, 0.01035688, 0.00774379,
                                      0.00568809, 0.0041053]
    # Define bin edges (as given)
    bin_edges = np.array([4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 
                      14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5])

    bins = bin_edges + 0.5
    print("Bins:", bins)
    
    bin_probabilities = np.zeros(len(bin_edges) - 1)

    for i in range(len(bin_probabilities)):
        # Probability of being in the current bin
        bin_probabilities[i] = weibull_min.cdf(bin_edges[i + 1], k, scale=c) - weibull_min.cdf(bin_edges[i], k, scale=c)

    bin_CDF_IEC = 1 - np.exp(-np.pi * (bin_edges / (2 * U_ave)) ** 2 )
    bin_probabilities_IEC = np.diff(bin_CDF_IEC)

    # Verifiy the results with the given bin probabilities
    assert np.allclose(bin_probabilities, bin_probabilities_verification, rtol=1e-5), "Bin probabilities do not match the verification values"
    print("Bin probabilities:", bin_probabilities)
    print('Bin probabilities IEC:', bin_probabilities_IEC)

    return bin_probabilities


def CombinedDelTest(bin_probability):
    m = 4
    DEL1 = np.array([43546.39930582, 55454.496635, 49645.98229272, 50795.02020986, 36507.24233848, 34222.73837354])
    combined_DEL_10min_verification = 46860.78173262

    # Step 1: Combine the DELs within the bin for a 10-minute equivalent DEL
    combined_DEL_10min = (np.mean(DEL1 ** m)) ** (1 / m)

    # Verifiy the results with the given equivalent DEL and lifetime fatigue load
    print(f"R_eq for the first wind speed bin: {combined_DEL_10min:.8f} kNm")
    print(f"Verification value of the 10-min combined DEL {combined_DEL_10min_verification } kNm")


def LifetimeFatigueLoad(bin_probabilities):
    # Given values
    m = 4  # Wöhler exponent
    n_eq = 10000000  # equivalent cycle count for 10-minute DEL
    N_T = 630720000  # total number of cycles in the lifetime

    # 10-minute DELs for each wind-speed bin
    DELs_10min_per_bin = np.array([
        46860.78173262, 60533.72498661, 59165.75222735, 49669.98028487,
        48249.13096179, 41188.70573419, 40341.05873765, 39077.21951565,
        37340.02415617, 36937.69987989, 39649.01740049, 43179.18330751,
        41795.27733207, 42207.91706777, 46612.72861644, 48856.4751121,
        49077.04325018, 54465.85842442, 53021.0998996, 54357.18578793
    ])

    # Calculate the lifetime cycles in each bin based on probabilities
    n_i_per_bin = N_T * bin_probabilities

    # Step 1: Adjust each DEL by the cycle ratio (n_i / n_eq) raised to the Wöhler exponent
    adjusted_DELs_per_bin = (n_i_per_bin / n_eq) * (DELs_10min_per_bin ** m)

    # Step 2: Sum the adjusted DELs
    sum_adjusted_DELs = np.sum(adjusted_DELs_per_bin)

    # Step 3: Take the 4th root (m-th root) to get the lifetime equivalent DEL
    lifetime_fatigue_load = sum_adjusted_DELs ** (1 / m)

    # Debug information for comparison
    Lifetime_fatigue_load_verification = 129826.22850291798 

    print(f"Calculated lifetime fatigue load: {lifetime_fatigue_load:.8f} kNm")
    print(f"Verification value for lifetime fatigue load: {Lifetime_fatigue_load_verification} kNm")

    return lifetime_fatigue_load


U_ave = 10
U_std = 2   
c = 2/np.sqrt(np.pi) * U_ave
k = 2

# SafetyFactorsTest()

bin_probabilities = ProbabilityBinsTest(c, k)

# CombinedDelTest(bin_probabilities[0])

# LifetimeFatigueLoad(bin_probabilities)



# Specify the path to your HDF5 file
file_path = '.\A4 Design Loads and AEP\Assignment\dtu_10mw_turb_stats.hdf5'
search_value = '5.946875487657857'
found_values = []

with h5py.File(file_path, 'r') as hdf_file:
    # Function to recursively search through the file for the target value
    def recursive_search(group):
        for key in group.keys():
            item = group[key]
            if isinstance(item, h5py.Dataset):
                # Check if the dataset contains the search value
                if search_value in str(item[:]):
                    found_values.append((key, item[:]))
            elif isinstance(item, h5py.Group):
                recursive_search(item)

    # Start recursive search from the root
    recursive_search(hdf_file)

print("Found values:", found_values)


"""
Load channel: TbFA [kNm]
Wöhler exponent: 4
DELs in first wind-speed bin:
 [43546.39930582 55454.496635   49645.98229272 50795.02020986
 36507.24233848 34222.73837354]
10-min DELs combined within each wind-speed bin: [46860.78173262 60533.72498661 59165.75222735 49669.98028487
 48249.13096179 41188.70573419 40341.05873765 39077.21951565
 37340.02415617 36937.69987989 39649.01740049 43179.18330751
 41795.27733207 42207.91706777 46612.72861644 48856.4751121
 49077.04325018 54465.85842442 53021.0998996  54357.18578793]
Lifetime fatigue load: 129826.22850291798
""" 
