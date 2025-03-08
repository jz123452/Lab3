# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 20:32:50 2025

@author: Lenovo
"""

import numpy as np 
import matplotlib.pyplot as plt
import scipy.optimize as so

antoine_et = [8.04494, 1554.30, 222.65] #https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Mask=4&Type=ANTOINE&Plot=on
antoine_h2o = [7.96681, 1668.21, 228.00]


def antoine(param, T):
    """
    Parameters
    ----------
    param : ARRAY
        Antoine parameters for the species.
    T : Float
        Temperature.

    Returns
    -------
    psat : Saturation Pressure.

    """
    A, B, C = param
    log_P_sat = A - B/(T + C)
    P_sat = 10**(log_P_sat)
    return P_sat

A_12 = 1.6798
A_21 = 0.9227
def van_laar(x_1):
    """
    Parameters
    ----------
    x_1 : Float
        Mole Fraction of species.

    Returns
    -------
    gamma1 : float
        Activity Coefficient (1)
    gamma2 : float
        Activity Coefficient (2)
    """
    x_2 = 1 - x_1
    log_gamma1 = A_12 * ((A_21 * x_2) / (A_12 * x_1 + A_21 * x_2))**2
    log_gamma2 = A_21 * ((A_12 * x_1) / (A_12 * x_1 + A_21 * x_2))**2
    
    gamma1 = np.exp(log_gamma1)
    gamma2 = np.exp(log_gamma2)
    return gamma1, gamma2

def vle(T, x_1, P):
    P1_sat = antoine(antoine_et,T)
    P2_sat = antoine(antoine_h2o,T)
    
    gamma_et, gamma_h2o = van_laar(x_1)
    
    y1 = (gamma_et * x_1 * P1_sat) / P
    y2 = (gamma_h2o * (1 - x_1) * P2_sat) / P
    
    return y1 + y2 -1
def calcy(T, x_1, P):
    P1_sat = antoine(antoine_et,T)
    gamma_et, gamma_h2o = van_laar(x_1)
    return (gamma_et * x_1 * P1_sat) / P

x_vals = np.linspace(0,1,50)
sol_T = []
sol_y = []
X1 = 0.5
P = 764.54 ##mmHg
T_guess = 50 #C
for x1 in x_vals:
    sol = so.root(vle, 70, args=(x1, P))  # Constrain T to 65-100°C)
    solved_T = sol.x[0]
    sol_T.append(solved_T)
    y1_value = calcy(solved_T,x1,P)
    sol_y.append(y1_value)
    
traytemps = [ 79.0,  80.1,  85.9,  87.1]
P = 764.54  # mmHg, assuming constant pressure for calculations

# Initialize lists to store results

tray_theoretical_x = []
# Calculate the theoretical mole fractions for each tray temperature
for temp in traytemps:
    # Use lambda to reverse argument order (x is now first)
    sol = so.root(
        lambda x, T_arg, P_arg: vle(T_arg, x, P_arg),  # Reverse args: x → T_arg → P_arg
        x0=[0.5],  # Initial guess for x1
        args=(temp, P),  # (T_arg, P_arg)
          # Restrict x1 to [0,1]
    )
    mol_frac = sol.x[0]
    tray_theoretical_x.append(mol_frac)

# ===== CORRECT SOLVER USE =====


        
reflux = 4
feed_frac_x = 0.065695 ##mol fr ethanol of feed stream
product_frac_x = 0.6733851
waste_frac_y = 0.026959
feed_dens = 0.974533
prod_dens = 0.8275333
waste_dens = 0.9859
##Tray Fracs
traydens = [np.mean([0.8375,0.8372,0.8326]),np.mean([0.8751,0.8717,0.8715]),np.mean([0.9565,0.9564,0.9561]),np.mean([0.9623,0.9639,0.9628])]
traytemps =[np.mean([25.1,23.7,23.9]),np.mean([21.9,22.4,22.0]),np.mean([24.0,23.9,24.3]),np.mean([21.2,21.1,21.3]),np.mean([26.1,24.3,24.4])]


##Mole wts
water_molwt = 18.015
ethanol_molwt = 46.068
feed_molwt = feed_frac_x * ethanol_molwt + (1-feed_frac_x)*water_molwt
waste_molwt = waste_frac_y * ethanol_molwt + (1-waste_frac_y)*water_molwt
prod_molwt = product_frac_x * ethanol_molwt + (1-product_frac_x)*water_molwt

molflux = 58.1/(feed_molwt/0.974533) ##ml/s -> mol/s
vol_prod = np.sum([7.5, 9.5, 6])
mol_prod = vol_prod/(prod_molwt/prod_dens)
F = molflux/mol_prod
#Molar Flow Rate = Volumetric Flow Rate / (Molar Mass / Density). 
oly =[]
olys = []


for i in range(len(x_vals)-1):
    oly.append(((reflux)/(reflux + 1))*x_vals[i] + product_frac_x/(reflux + 1))
    olys.append(((reflux + F)/(reflux+1))*x_vals[i] + waste_frac_y*((1-F/(reflux+1))))
olys.append(olys[-1])
oly.append(oly[-1])

trayfrac_exp = [0.649,
0.4403,
0.13105,
0.10917
]

ys = np.linspace(0,1,10)

plt.figure(0, figsize=(6, 4))
ax = plt.gca()
ax.set_ylim(0,1)
plt.plot(x_vals, sol_y, label='Vapor Mole Fraction (y₁)')
plt.plot(x_vals, x_vals, 'k-', label='y₁ = x₁ line')
#plt.plot(trayfracs,trayfracs,'bo',label='Tray 2, 4, 6 and 8 ethanol mole fraction')
plt.plot(product_frac_x*np.ones(len(ys)),ys,'k-.' )
plt.plot(waste_frac_y, waste_frac_y , 'b.')
plt.plot(feed_frac_x, feed_frac_x , 'b.')
plt.plot(product_frac_x, product_frac_x , 'b.')
plt.plot(waste_frac_y*np.ones(len(ys)),ys,'k--' )
plt.plot(feed_frac_x*np.ones(len(ys)),ys,'k--' )
plt.plot(x_vals[3:33], oly[3:33] , 'r-', label='operating line')
plt.plot(0, -10 , 'r-', label='operating line (stripping)')
plt.plot(0, -10 , 'g-', label='q - line')
plt.tight_layout()

#plt.title('VLE x-y Diagram for EtOH-Water at P = 764.54 mmHg')
plt.xlabel('Liquid Mole Fraction (x₁, Ethanol)')
plt.ylabel('Vapor Mole Fraction (y₁, Ethanol)')
plt.legend()
plt.grid(True)
# Plot T-x-y diagram
plt.figure(1, figsize=(6, 4))
plt.plot(x_vals, sol_T, label='Temperature vs. x₁')
plt.plot(sol_y, sol_T, label='Temperature vs. y₁')
plt.title('Temperature-Composition Diagram for EtOH-Water at P = 764.54 mmHg')
plt.xlabel('Mole Fraction')
plt.ylabel('Temperature (°C)')
plt.legend()
plt.grid(True)
plt.show()
plt.figure(0, figsize=(6, 4))
plt.plot(x_vals, sol_y, label='Vapor Mole Fraction (y₁)')
plt.plot(x_vals, x_vals, 'k--', label='y₁ = x₁ line')
plt.title('VLE x-y Diagram for EtOH-Water at P = 764.54 mmHg')
plt.xlabel('Liquid Mole Fraction (x₁, Ethanol)')
plt.ylabel('Vapor Mole Fraction (y₁, Ethanol)')
plt.legend()
plt.grid(True)
# Plot T-x-y diagram
plt.figure(1, figsize=(6, 4))
plt.plot(x_vals, sol_T, label='Temperature vs. x₁')
plt.plot(sol_y, sol_T, label='Temperature vs. y₁')
plt.title('Temperature-Composition Diagram for EtOH-Water at P = 764.54 mmHg')
plt.xlabel('Mole Fraction')
plt.ylabel('Temperature (°C)')
plt.legend()
plt.grid(True)
plt.show()