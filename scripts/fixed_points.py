#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 20:36:54 2018

@author: booort
"""


import scipy as sp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# Definition of constants
# from Lai-Schaffer classic model

Shh = 15*2  # Shh quantity [0,30]
k_shh = 1  # dissociation constant shh-ptc bindings [0.58,2.0]
k_ptc = 8.3*10**-2  # 1/2maximal concentration of ptc which inhibits smo signlaing
k_deg = 0.009*60  # degradation constant for all Gli related proteins
k_g3rc = 0.012*60  # rate constant for the conversion to signal strenGh
r_g3b = 1.6*10**-1*60  # basal rate of Gli3 synthesis
K_g3rc = 0.1  # sensitivity constant of the conversion to signal strenGh
k_deg_p = 0.09*60  # degradation rate constant for Ptc [0.045,0.071]
# --------------
K1 = 8.3*10**-1
K2 = 8.3*10**-1
c = 1
e = 0.5
r = 0.2
v_max = 2.4*10**-1*60
r_bas = v_max/100
v_maxp = 7.5*10**-1*60
r_basp = v_maxp/100
k_cc = v_maxp*k_deg/v_max*k_deg_p

def gli_curve(Gli):
    Ptc = k_cc*Gli
    Signal = (1+(Shh/k_shh))/(1+(Shh/k_shh)+(Ptc/k_ptc))
    Gli3 = (r_g3b*(K_g3rc+Signal))/(k_deg*(K_g3rc+Signal)+k_g3rc)
    Gli3R = (r_g3b)/(k_deg*Gli)-Gli3
    Basal = (1-1/c+1/c*(1+c*Gli*K1**-1+c*Gli3*K2**-1+r*c*Gli3R*K2**-1)**3)/(1-1/c+1/c*(1+c*Gli*K1**-1+c*Gli3*K2**-1+c*Gli3R*K2**-1)**3)
    Promoter = (1/c*(e+c*Gli*K1**-1+c*Gli3*K2**-1+e*r*c*Gli3R*K2**-1)**3-1/c*(e+e*r*c*Gli3R*K2**-1)**3)/(1-1/c+1/c*(1+c*Gli*K1**-1+c*Gli3*K2**-1+c*Gli3R*K2**-1)**3)
    return (Promoter+0.01*Basal)
    
    
mesh_size=0.001
Gli = sp.arange(0.01, 40.0, mesh_size)
Gli_1 = sp.arange(0.01, 40.0, mesh_size)
fig, ax = plt.subplots()
ax.plot(Gli_1, k_deg*Gli_1/v_max,label='Gli=Gli')
ax.plot(Gli, gli_curve(Gli), label='Gli_curve')
plt.show()



