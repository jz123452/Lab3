# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 20:21:42 2025

@author: Lenovo
"""

import numpy as np
import matplotlib.pyplot as pyplot
import scipy.optimize 
import glob
import scipy.stats
import csv
import scipy.signal
import pint
import pandas as pd


u = pint.UnitRegistry()

def open_file(file_name):
    with open(file_name,'r') as f:
        raw_data = list(csv.reader(f))
        f.close()
        return raw_data
table6 = open_file('Table 6.csv')
table6 = table6[:-3]
table6 = np.array(table6, dtype=np.float64)
table1 = open_file('Table 1.csv')
table1 = table1[5:]
table1 = pd.DataFrame(table1)
table1.replace('', 0, inplace=True)




def correctsg(d,T):##calculate what density a glass hydrometer would have indicated were the glass calibrated at 60◦F but used at T
    alpha = 25 * 10**-6   ##degC^-1
    return d*(1+alpha*(T-15.5555))/0.99904


def interpolate6(d, T):
    target_value = correctsg(d, T)  # Single calculation for the target value
    nearest_lower_value = None
    nearest_upper_value = None
    min_diff_lower = float('inf')
    min_diff_upper = float('inf')
    index_lower = None
    index_upper = None

    # Iterate through table6 to find nearest lower and upper values
    for i in range(len(table6)):
        current_value = table6[i, 1]  # Access the second column
        
        diff = abs(current_value - target_value)
        
        # Find the nearest lower value
        if current_value < target_value and diff < min_diff_lower:
            nearest_lower_value = current_value
            min_diff_lower = diff
            index_lower = i
        
        # Find the nearest upper value
        elif current_value > target_value and diff < min_diff_upper:
            nearest_upper_value = current_value
            min_diff_upper = diff
            index_upper = i

    # Ensure that both lower and upper values are found before calculating apparent_proof
    if index_lower is not None and index_upper is not None:
        apparent_proof = table6[index_lower, 0] + (
            (target_value - nearest_upper_value) / (nearest_lower_value - nearest_upper_value) * 
            (table6[index_upper, 0] - table6[index_lower, 0])
        )
        print("Apparent Proof:", apparent_proof)  # Debugging output
        return apparent_proof
    else:
        print("Error: Could not find appropriate lower and upper values for interpolation.")
        return None  # Or raise an exception or return a specified value
densities = np.array([0.921866667,
0.8534,
0.969866667,
0.958766667,
0.9743,
0.827233333,
0.874366667,
0.904866667,
0.978333333,

])

densuncertanties = np.array([0.003607468,
0.005066667,
0.002031668,
0.00433757,
0.000405658,
0.000309826,
0.000309826,
0.00223419,
0.000309826,
0.006791996,
])
temps = np.array([23.66666667,
22.73333333,
26.2,
25.9,
25.23333333,
23.66666667,
24.4,
25.93333333,
25.46666667,
])
proofs = []
for i in range(len(densities)):
    proofs.append(interpolate6(densities[i],temps[i]))
proofs = np.array(proofs,dtype=float)
tempuncertainties = np.array([0.117103375,
0.117103375,
0.202828995,
0.202828995,
0.117103375,
0.117103375,
8.82543E-15,
0.309826407,
0.117103375,
0.202828995,
])
def trueproof(AP,T):
    Tf = T*(9/5) + 32 ## Convert to Farenheit
    Tl = int(np.floor(Tf) + 1)
    Th = int(np.ceil(Tf) + 1)
    APl = int(np.floor(AP) + 1)
    APh = int(np.ceil(AP) + 1)
    APlTl = float(table1.iloc[APl][Tl])
    APlTh = float(table1.iloc[APl][Th])
    APhTl = float(table1.iloc[APh][Tl])
    APhTh = float(table1.iloc[APh][Th])
    dfdc = (APhTl - APhTh) / (APh - APl)
    dfdt = (APlTh - APlTl) / (Th - Tl)
    trueproof = APlTl + (AP - APl)*dfdc + (Tf-Tl)*dfdt
    
    return trueproof
    
trueproofs = []
for i in range(len(proofs)):
    trueproofs.append(trueproof(proofs[i], temps[i]))

