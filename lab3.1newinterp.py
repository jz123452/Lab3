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
   
    Wvv = wbvL + ((CP - CP_L)/(CP_H - CP_L+1e-12))*(wbvH - wbvL)
    Avv = CP/2
    molefret = molefret = (Avv * 0.79313/46.06844) / ((Avv * 0.79313/46.06844) + (Wvv * 0.99904/18.01528))
    return molefret



#EXPERIMENTAL DATA (3.2)
etden = np.mean([0.8164,0.8168,0.8151]) #g/cm @ 20c
ettemp = np.mean([26.2,26.1,26.2])
waden = np.mean([0.9879,0.9892,0.9886])
watemp = np.mean([16.2,14.5,18.7])
feden = np.mean([0.977,0.9786,0.9786])
fetemp = np.mean([26.1,24.3,24.4])
densities = [etden,waden,feden]
temps = [ettemp,watemp,fetemp]
traydens = [np.mean([0.8288,0.8312,0.8301]),np.mean([0.8716,0.8700,0.8688]),np.mean([0.9562,0.9560,0.9557]),np.mean([0.9683,0.9667,0.9681])]
traytemps =[np.mean([25.1,23.7,23.9]),np.mean([21.9,22.4,22.0]),np.mean([24.0,23.9,24.3]),np.mean([21.2,21.1,21.3]),np.mean([26.1,24.3,24.4])]
###THEORETICAL DATA

densities20 = np.linspace(0.7894,0.99819,100) ## 20C
densities25 = np.linspace(0.7850,0.99705,100) ## 25C
densities60 = np.linspace(0.7932,0.99898,100) ## @ 60F

temps20  = 20*np.ones(len(densities20))
temps60 = 15.55556*np.ones(len(densities25))
temps25 = 25**np.ones(len(densities60))


proofs = []
proofs60 = []
proofs25 = []
proofs20 = []
trayproofs = []
for i in range(len(densities)):
    proofs.append(interpolate6(densities[i],temps[i]))
for i in range(len(densities60)):
    proofs60.append(interpolate6(densities60[i],temps60[i]))
    proofs25.append(interpolate6(densities25[i],temps25[i]))
    proofs20.append(interpolate6(densities20[i],temps20[i]))
for i in range(len(traydens)):
    trayproofs.append(interpolate6(traydens[i],traytemps[i]))
proofs = sorted(np.array(proofs,dtype=float), reverse=True)##remove sorted for experimental data
proofs60 = sorted(np.array(proofs60,dtype=float), reverse=True)
proofs25 = sorted(np.array(proofs25,dtype=float), reverse=True)
proofs20 = sorted(np.array(proofs20,dtype=float), reverse=True)
trayproofs = sorted(np.array(trayproofs,dtype=float), reverse=True)

trueproofs = []
trueproofs60 = []
trueproofs25 = []
trueproofs20 = []
traytrueproofs = []
for i in range(len(proofs)):
    trueproofs.append(trueproof(proofs[i], temps[i]))
    
for i in range(len(proofs25)):
    trueproofs60.append(trueproof(proofs60[i], temps60[i]))
    trueproofs25.append(trueproof(proofs25[i], temps25[i]))
    trueproofs20.append(trueproof(proofs20[i], temps20[i]))
for i in range(len(trayproofs)):
    traytrueproofs.append(trueproof(trayproofs[i], traytemps[i]))
molefrs = []
molefrs60 = []
molefrs25 = []
molefrs20 = []
traymolefrs = []
for i in range(len(trueproofs)):
    molefrs.append(molefractions(trueproofs[i]))
for i in range(len(trueproofs25)):
    molefrs60.append(molefractions(trueproofs60[i]))
    molefrs25.append(molefractions(trueproofs25[i]))
    molefrs20.append(molefractions(trueproofs20[i]))
for i in range(len(traytrueproofs)):
    traymolefrs.append(molefractions(traytrueproofs[i]))
    
######### EQUATION 2
### 3.1 Data


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
water_density_table = [
    {"temp_C":4, "density":0.99997},
    {"temp_C":15.56, "density":0.99904},  # 60°F base
    {"temp_C":20, "density":0.99823},
    {"temp_C":25, "density":0.99705},
    # ... other temp/density pairs
]
def get_water_density(temp_C):
    return np.interp(temp_C, 
                    [x["temp_C"] for x in water_density_table],
                    [x["density"] for x in water_density_table])
pycnotemps = [24.4,
24.46666667,
25.33333333,
25.66666667,
25.26666667,
22.8,
23.73333333,
22.93333333,
24.93333333,
22.66666667,
]
m_e = 35.38 ##average empty pycno mass(from excel)
m_w = 60.80
m_f = [59.00,
57.11,
59.11,
59.89,
60.12,
56.44,
57.57,
58.31,
55.35,
59.00
]
pycnodens = []
def pycnopf(temp,m_e,m_w,m_f):
    d_w = get_water_density(temp)
    ### INPUTS, 
    ## d_water = f(t), plug in temp
    ## m_e = mass of empty pycnometer, grams
    ## m_w mass of pycnometer filled with water
    ## m_f final pycnometer mass
    d_f = ((m_f - m_e)/(m_w - m_e))*d_w
    return d_f
for i in range(len(m_f)):
    pycnodens.append(pycnopf(pycnotemps[i],m_e,m_w,m_f[i]))
### error propogation
import numpy as np

def propagate_error(m_f, m_e, m_w, d_m_f, d_m_e, d_m_w, T, delta_T):
    """Propagate uncertainties from mass and temperature (fixed)."""
    # Convert all inputs to numpy arrays
    m_f = np.asarray(m_f)
    m_e = np.asarray(m_e)
    m_w = np.asarray(m_w)
    T = np.asarray(T)
    d_m_f = np.asarray(d_m_f)
    d_m_e = np.asarray(d_m_e)
    d_m_w = np.asarray(d_m_w)
    
    # Calculate water density and its temperature derivative
    d_w = get_water_density(T)  # Assume get_water_density uses np.interp
    
    # Numerical derivative of d_w with respect to T
    epsilon = 0.01
    d_w_plus = get_water_density(T + epsilon)
    d_w_minus = get_water_density(T - epsilon)
    ddw_dT = (d_w_plus - d_w_minus) / (2 * epsilon)
    
    # Partial derivatives (ensure all terms are arrays)
    denominator = m_w - m_e
    df_dm_f = 1 / denominator
    df_dm_e = -(m_f - m_e) / denominator**2
    df_dm_w = -(m_f - m_e) / denominator**2
    df_dd_w = (m_f - m_e) / denominator
    
    # Temperature uncertainty contribution (scalar * array → array)
    delta_d_w = ddw_dT * np.asarray(delta_T)
    
    # Propagate uncertainties (all terms now arrays)
    uncertainties = np.sqrt(
        (df_dm_f * d_m_f)**2 +
        (df_dm_e * d_m_e)**2 +
        (df_dm_w * d_m_w)**2 +
        (df_dd_w * delta_d_w)**2  # Array multiplication
    )
    
    # Uncertainty of the mean density
    return np.sqrt(np.sum(uncertainties**2)) / len(uncertainties)

densuncs = []
dm_f = [0.005366351,
0.001171034,
0.001171034,
0.003098264,
0,
0.003513101,
0.00405658,
0.002342067,
0.021465403,
0.005366351
]
dm_e = 0.003098264
dm_w = 1.77 * 10**(-14)

for i in range(len(m_f)):
    densuncs.append(propagate_error(m_f[i],m_e,m_w,dm_f,dm_e,dm_w,pycnotemps[i],tempuncertainties[i]))
    
plt.figure('1')
plt.plot(densities60,molefrs60,'o', markersize=2,label='60°F')
plt.plot(densities25,molefrs25,'o', markersize=2,label='25°C')
plt.plot(densities20,molefrs20,'o', markersize=2,label='20°C')
plt.legend()
plt.grid(True)
plt.ylabel('Mole Fraction of Ethanol')
plt.xlabel('Density of Mixture (g/mL)')
plt.figure('2')
plt.ylabel('True Proof')
plt.xlabel('Density of Mixture (g/mL)')
plt.plot(densities20,trueproofs20,'o', markersize=2,label='20°C')
plt.plot(densities60,trueproofs60,'o', markersize=2,label='60°F')
plt.plot(densities25,trueproofs25,'o', markersize=2,label='25°C')
plt.legend()
plt.grid(True)
