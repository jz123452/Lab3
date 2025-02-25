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