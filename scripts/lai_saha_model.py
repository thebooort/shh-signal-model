#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  1 21:01:45 2018

@author: booort
"""


from math import *
import scipy as sp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# Definition of constants
# from Lai-Schaffer classic model

Shh = 15  # Shh quantity [0,30]
k_shh = 0.58  # dissociation constant shh-ptc bindings [0.58,2.0]
k_ptc = 8.3*10**-2  # half maximal concentration of ptc which inhibits smo signlaing
k_deg = 0.009  # degradation constant for all Gli related proteins
k_g3rc = 0.012  # rate constant for the conversion to signal strenGh
r_g3b = 1.6*10**-1  # basal rate of Gli3 synthesis
K_g3rc = 0.1  # sensitivity constant of the conversion to signal strenGh
k_deg_p = 0.09  # degradation rate constant for Ptc [0.045,0.071]
# --------------
K1 = 8.3*10**-1
K2 = 8.3*10**-1
c = 1
e = 0.5
r = 0.2
v_max = 2.4*10**-1 
r_bas = v_max/100
v_maxp = 7.5*10**-1
r_basp = v_maxp/100

def Signal(Ptc):

    return (1+Shh/k_shh)/(1+Shh/k_shh+Ptc/k_ptc)

def Promoter(Gli,Gli3,Gli3R):
    ((Gli3*K1+Gli*K2)*(3*e**2*K1**2*K2**2+3*c*e*K1*K2*(Gli3*K1+Gli*K2+2*e*Gli3R*K1*r)+c**2*(Gli3**2*K1**2+Gli**2*K2**2+3*e*Gli*Gli3R*K1*K2*r + 3*e**2*Gli3R**2*K1**2*r**2 + Gli3*K1*(2*Gli*K2 + 3*e*Gli3R*K1*r))))/(3*c*K1*K2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**2 + c**2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**3 +K1**2*K2**2*(3*Gli3*K1 + 3*Gli3R*K1 + (3*Gli + K1)*K2))
    
    
def Basal(Gli,Gli3,Gli3R):
    (3*c*K1*K2*(Gli3*K1+ Gli*K2 + Gli3R*K1*r)**2 + c**2*(Gli3*K1 + Gli*K2 + Gli3R*K1*r)**3 + K1**2*K2**2*(3*Gli3*K1 + 3*Gli*K2 + K1*(K2+ 3*Gli3R*r)))/ (3*c*K1*K2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**2 + c**2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**3 + K1**2*K2**2*(3*Gli3*K1 + 3*Gli3R*K1 + (3*Gli +K1)*K2))








def shh_evolution_system(X, t):
    Gli, Gli3, Gli3R, Ptc = X

    dGli_dt = v_max*Promoter(Gli, Gli3, Gli3R)+r_bas*Basal(Gli, Gli3, Gli3R)-k_deg*Gli
    dGli3_dt = r_g3b/Ptc-Gli3*(k_deg+k_g3rc/(K_g3rc+Signal(Ptc)))
    dGli3R_dt = Gli3*(k_g3rc/(K_g3rc+Signal(Ptc)))-k_deg*Gli3R
    dPtc_dt = v_maxp*Promoter(Gli, Gli3, Gli3R)+r_basp*Basal-k_deg_p*Ptc

    return dGli_dt, dGli3_dt, dGli3R_dt, dPtc_dt
