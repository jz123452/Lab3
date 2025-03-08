# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 20:21:42 2025

@author: Lenovo
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize 
import glob
import scipy.stats
import csv
import scipy.signal

import pandas as pd




def open_file(file_name):
    with open(file_name,'r') as f:
        raw_data = list(csv.reader(f))
        f.close()
        return raw_data
table6 = open_file('TTB_6.csv')
table6 = table6[1:]
table6 = np.array(table6, dtype=np.float64)
table1 = open_file('Table 1.csv')
table1 = table1[5:]
table1 = pd.DataFrame(table1)
table1.replace('', 0, inplace=True)




def correctsg(d,T):##calculate what density a glass hydrometer would have indicated were the glass calibrated at 60◦F but used at T
    alpha = 25 * 10**-6   ##degC^-1
    return d*(1+alpha*(T-15.5555))/0.99904


def interpolate6(d, T):
    target_value = correctsg(d, T)
    nearest_lower = None
    nearest_upper = None
    min_diff_lower = float('inf')
    min_diff_upper = float('inf')
    index_lower = None
    index_upper = None

    # Find nearest bounds in table6
    for i, row in enumerate(table6):
        current_value = row[3]  # Assume column 3 holds the S.G. values
        if current_value < target_value and (target_value - current_value) < min_diff_lower:
            nearest_lower = current_value
            index_lower = i
            min_diff_lower = target_value - current_value
        elif current_value > target_value and (current_value - target_value) < min_diff_upper:
            nearest_upper = current_value
            index_upper = i
            min_diff_upper = current_value - target_value

    # Edge case handling
    if nearest_upper is None:  # Target above max in table
        nearest_upper = table6[-1, 3]
        index_upper = len(table6) - 1
        nearest_lower = nearest_upper
        index_lower = index_upper
    elif nearest_lower is None:  # Target below min in table
        nearest_lower = table6[0, 3]
        index_lower = 0
        nearest_upper = nearest_lower
        index_upper = index_lower

    # Prevent division by zero
    if nearest_upper == nearest_lower:
        apparent_proof = table6[index_lower, 0]
    else:
        fraction = (target_value - nearest_upper) / (nearest_lower - nearest_upper)
        apparent_proof = table6[index_lower, 0] + fraction * (table6[index_upper, 0] - table6[index_lower, 0])

    return apparent_proof

def trueproof(AP,T):
    
    Tf = T*(9/5) + 32 ## Convert to Farenheit
    Tl = int(np.floor(Tf) )
    Th = int(np.ceil(Tf) )
    APl = int(np.floor(AP) )
    APh = int(np.ceil(AP) )
    APlTl = float(table1.iloc[APl][Tl])
    APlTh = float(table1.iloc[APl][Th])
    APhTl = float(table1.iloc[APh][Tl])
    APhTh = float(table1.iloc[APh][Th])
    if APh == APl:
        dfdc = 0
    else:
        dfdc = (APhTl - APhTh) / (APh - APl)
    if Th == Tl:
        dfdt = 0
    else:
        dfdt = (APlTh - APlTl) / (Th - Tl)
    trueproof = APlTl + (AP - APl)*dfdc + (Tf-Tl)*dfdt
    
    return trueproof
    

# Modified molefractions function
def molefractions(CP):
    CP_L = int(np.floor(CP))
    CP_H = int(np.ceil(CP))
    
    # Check indices to prevent out-of-bounds
    max_index = len(table6) - 1  # Maximum valid index
    CP_L = min(max(CP_L, 0), max_index)
    CP_H = min(max(CP_H, 0), max_index)
    
    # REMOVED +1 OFFSET
    abvL = table6[CP_L][1]
    abvH = table6[CP_H][1]
    wbvL = table6[CP_L][2]
    wbvH = table6[CP_H][2]
   
    Wvv = wbvL + ((CP - CP_L)/(CP_H - CP_L))*(wbvH - wbvL)
    Avv = CP/2
    molefret = molefret = (Avv * 0.79313/46.06844) / ((Avv * 0.79313/46.06844) + (Wvv * 0.99904/18.01528))
    return molefret


#EXPERIMENTAL DATA
#etden = np.mean([0.8164,0.8168,0.8151]) #g/cm @ 20c
#ettemp = np.mean([26.2,26.1,26.2])
#waden = np.mean([0.9879,0.9892,0.9886])
#watemp = np.mean([16.2,14.5,18.7])
#feden = np.mean([0.977,0.9786,0.9786])
#fetemp = np.mean([26.1,24.3,24.4])
 ## @ 60F
densities = [np.mean([0.9741,0.9747,0.9748]),np.mean([0.9892,0.9841,0.9844]),np.mean([0.8292,0.8275,0.8259])]
 ## 25C
temps = [np.mean([26.3,27.2,26.0]),np.mean([18.8,18.7,19.1]),np.mean([20.3,20.6,22.5])]

traydens = [np.mean([0.8375,0.8372,0.8326]),np.mean([0.8751,0.8717,0.8715]),np.mean([0.9565,0.9564,0.9561]),np.mean([0.9623,0.9639,0.9628])]
traytemps = [np.mean([20.1,20.7,22.8]),np.mean([21.5,21.8,21.9]),np.mean([21.4,20.1,20.5]),np.mean([19.8,19.6,20.2])]
#densities = np.array([0.921866667,
#0.8534,
#0.969866667,
#0.958766667,
#0.9743,
#0.827233333,
#0.874366667,
#0.904866667,
#0.978333333,
#])

#densuncertanties = np.array([0.003607468,
#0.005066667,
#0.002031668,
#0.00433757,
#0.000405658,
#0.000309826,
#0.000309826,
#0.00223419,
#0.000309826,
#0.006791996,
#])
#temps = np.array([23.66666667,
#22.73333333,
#26.2,
#25.9,
#25.23333333,
#23.66666667,
#24.4,
#25.93333333,
#25.46666667,
#])

trayproofs = []
proofs = []
for i in range(len(densities)):
    
    proofs.append(interpolate6(densities[i],temps[i]))
##remove sorted for experimental data
proofs = sorted(np.array(proofs,dtype=float), reverse=True)
for i in range(len(traydens)):
    
    trayproofs.append(interpolate6(traydens[i],traytemps[i]))
##remove sorted for experimental data
trayproofs = sorted(np.array(trayproofs,dtype=float), reverse=True)
#tempuncertainties = np.array([0.117103375,
#0.117103375,
#0.202828995,
#0.202828995,
#0.117103375,
#0.117103375,
#8.82543E-15,
#0.309826407,
#0.117103375,
#0.202828995,
#])

traytrueproofs = []
trueproofs = []
for i in range(len(trayproofs)):
    
    traytrueproofs.append(trueproof(trayproofs[i], traytemps[i]))
for i in range(len(proofs)):
    
    trueproofs.append(trueproof(proofs[i], temps[i]))

molefrs = []
traymolefrs = []

for i in range(len(traytrueproofs)):
  
    traymolefrs.append(molefractions(traytrueproofs[i]))
for i in range(len(trueproofs)):
  
    molefrs.append(molefractions(trueproofs[i]))
plt.figure('1')


plt.legend()
plt.grid(True)
plt.ylabel('Mole Fraction of Ethanol')
plt.xlabel('Density of Mixture (g/mL)')
plt.figure('2')
plt.ylabel('True Proof')
plt.xlabel('Density of Mixture (g/mL)')


plt.legend()
plt.grid(True)
