# -*- coding: utf-8 -*-
from lacbox.io import load_stats, load_oper
import matplotlib.pyplot as plt
import numpy as np

'''AEP Q1 - Weibull parameters'''
U_mean = 8
sigma_U = 5

C = (2/np.sqrt(np.pi))*U_mean
k = (sigma_U/U_mean)**(-1.086)
print("Weibull scale parameter (C):", round(C, 3))
print("Weibull shape parameter (k):", round(k, 3))
print('Note: These expressions are valid only when 1.6 < k < 3')

'''AEP Q2 - Simple power curve AEP'''
V_bins = [0,5,8,12,14,25,100] #100 is just a high random number
P_bins = [0,20,35,40,45,0]

V_prob = []
for i in range(0,len(V_bins)-1):
    V_prob.append(np.exp(-(V_bins[i]/C)**k) - np.exp(-(V_bins[i+1]/C)**k))

P = [p * v for p, v in zip(P_bins, V_prob)]
P_tot = sum(P)

print("Power generation before reliability:", round(P_tot,3), 'kW')

R = 0.95
P_tot_reliability = P_tot * R
AEP = P_tot_reliability * 365.25 * 24/1000

print('AEP:', str(round(AEP,2)), 'MWh')