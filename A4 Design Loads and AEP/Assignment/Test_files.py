import numpy as np
from scipy.stats import weibull_min

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



# Example DELs and cycle counts for a wind speed bin
DELs = np.array([5579.56307297, 7045.70884184, 5387.71133704, 5410.98404634, 4523.11438317, 4538.78189934])

U_ave = 10
U_std = 2   

c = 2/np.sqrt(np.pi) * U_ave
k = 2
# k = (U_std/U_ave) ** -1.086

# Define bin edges (as given)
bin_edges = np.array([4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 
                      14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5])

# Calculate bin probabilities using the Weibull CDF
bin_probabilities = np.zeros(len(bin_edges) - 1)

for i in range(len(bin_probabilities)):
    # Probability of being in the current bin
    bin_probabilities[i] = weibull_min.cdf(bin_edges[i + 1], k, scale=c) - weibull_min.cdf(bin_edges[i], k, scale=c)

# Normalize probabilities to sum to 1 (just in case)
bin_probabilities /= np.sum(bin_probabilities)

# Print the probabilities for each bin
print("Bin probabilities:", bin_probabilities)
print("Sum of bin probabilities:", np.sum(bin_probabilities))


Wohler_exponent = 4
n_eq = 10000000  # equivalent cycle count
N_T = 630720000  # total number of cycles in the lifetime

ref = 5617.84563021

# Step 1: Adjust each DEL by the cycle ratio (N_T / n_eq) raised to the Wöhler exponent
adjusted_DELs = (N_T / n_eq) * (DELs ** Wohler_exponent)

# Step 2: Sum the adjusted DELs
sum_adjusted_DELs = np.sum(adjusted_DELs)

# Step 3: Take the 4th root (m-th root) to get the equivalent DEL
R_eq = sum_adjusted_DELs ** (1 / Wohler_exponent)

print(f"Equivalent DEL (R_eq) for the wind speed bin: {R_eq:.8f} kNm")
print(f"Reference value: {ref}")
print(f"Difference: {abs(R_eq - ref)} kNm")