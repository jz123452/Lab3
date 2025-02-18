# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 00:05:59 2025
Lab 3,1
@author: Lenovo
"""
import pint 
import numpy as np
import matplotlib.pyplot as plt

u = pint.UnitRegistry()

def correctdensity(d,T):##calculate what density a glass hydrometer would have indicated were the glass calibrated at 60◦F but used at T
    alpha = 25 * 10**-6   ##degC^-1
    return d*(1+alpha*(T-15.5555))
def proof(d,T):
    SG = correctdensity(d,T)/0.99904##SG of water @ 60F
    proof = -991.99*SG + 1014 ## From Linear Interpolation of Table 6: Slope = [Proof/SG]
    return proof

def trueproof(P,T):
    proofdown = np.floor(P)
    tdown = np.floor(T)
    
densities = np.array([0.921866667,
0.8534,
0.969866667,
0.958766667,
0.9743,
0.827233333,
0.874366667,
0.904866667,
0.978333333,
0.786066667,
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
23.5,
])
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
proofs = proof(densities,temps)