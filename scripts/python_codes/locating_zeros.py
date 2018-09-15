#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 18:37:30 2018

This script is a visual representation of our search of zeros. It also contains a few tries to see times and basic results

@author: booort
"""

import scipy as sp
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import time
start_time = time.time()
# Definition of constants
c = 1  # positive constant, Greater than 1 implies cooperativity, less than 1 anti-cooperativity

a_Gli = 4.35  # transcriptional activation intensity for Gli
a_Gli3 = 4.35  # transcriptional activation intensity for Gli3
r_Gli3R = 5*10**-5  # transcriptional repression intensity for Gli

k_Gli = 9*10**1  # dissociation constant of activators for Gene enhancers
k_Gli3 = 9*10**1  # dissociation constant of activators for Gene enhancers
k_Gli3R = 9*10**1   # dissociation constant of repressors for Gene enhancers

k_RNAP = 1  # RNA polymerase binding affinity
RNAP = 1  # RNA polymerase concentration

c_b = 0.26*60  # BEWARE constant
c_b1 = 3.15*60
# from Lai-Schaffer classic model shared with BEWARE

Shh = 1.1*1  # Shh quantity [0,30]
k_shh = 1  # dissociation constant shh-ptc bindings [0.58,2.0]

k_ptc = 8.3*10**-2  # 1/2maximal concentration of ptc which inhibits smo signlaing
k_deg = 0.009*60  # degradation constant for all Gli related proteins

k_g3rc = 0.012*60  # rate constant for the conversion to signal strenGh
r_g3b = 60*0.16  # basal rate of Gli3 synthesis

K_g3rc = 0.1*10**0  # sensitivity constant of the conversion to signal strenGh
k_deg_p = 0.09*60  # degradation rate constant for Ptc [0.045,0.071]
# Exclusively from lai-saha
K1 = 8.3*10**-1
K2 = 8.3*10**-1
c = 2
e = 0.5
r = 0.2
v_max = 60*2.4*10**-1
r_bas = v_max/100
v_maxp = 60*7.5*10**-1
r_basp = v_maxp/100
k_cc = (v_maxp*k_deg)/(v_max*k_deg_p)


def gli_curve(Gli):
    Ptc = k_cc*Gli
    Signal = (1+(Shh/k_shh))/(1+(Shh/k_shh)+(Ptc/k_ptc))
    Gli3 = (r_g3b*(K_g3rc+Signal))/((k_deg*(K_g3rc+Signal)+k_g3rc)*Gli)
    Gli3R = (r_g3b)/(k_deg*Gli)-Gli3
    Basal = (3*c*K1*K2*(Gli3*K1 + Gli*K2 + Gli3R*K1*r)**2 + c**2*(Gli3*K1 + Gli*K2 + Gli3R*K1*r)**3 + K1**2*K2**2*(3*Gli3*K1 + 3*Gli*K2 + K1*(K2 + 3*Gli3R*r))
             ) / (3*c*K1*K2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**2 + c**2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**3 + K1**2*K2**2*(3*Gli3*K1 + 3*Gli3R*K1 + (3*Gli + K1)*K2))
    Promoter = ((Gli3*K1+Gli*K2)*(3*e**2*K1**2*K2**2+3*c*e*K1*K2*(Gli3*K1+Gli*K2+2*e*Gli3R*K1*r)+c**2*(Gli3**2*K1**2+Gli**2*K2**2+3*e*Gli*Gli3R*K1*K2*r + 3*e**2*Gli3R**2*K1**2*r**2 +
                                                                                                       Gli3*K1*(2*Gli*K2 + 3*e*Gli3R*K1*r))))/(3*c*K1*K2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**2 + c**2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**3 + K1**2*K2**2*(3*Gli3*K1 + 3*Gli3R*K1 + (3*Gli + K1)*K2))
    return v_max*(Promoter+0.01*Basal)/k_deg


def gli_curve_1(Gli):
    Ptc = (0.89*60*k_deg)/(k_deg_p*c_b)*Gli
    Signal = (1+(Shh/k_shh))/(1+(Shh/k_shh)+(Ptc/k_ptc))
    Gli3 = (r_g3b*(K_g3rc+Signal))/((k_deg*(K_g3rc+Signal)+k_g3rc)*Gli)
    Gli3R = (r_g3b)/(k_deg*Gli)-Gli3
    F_reg_nt_coop = (1-1/c+1/c*(1+a_Gli*c*Gli/k_Gli+a_Gli3*c*Gli3/k_Gli3+r_Gli3R*c *
                                Gli3R/k_Gli3R)**3)/(1-1/c+1/c*(1+c*Gli/k_Gli+c*Gli3/k_Gli3+c*Gli3R/k_Gli3R)**3)
    beware = c_b/(1+k_RNAP/(RNAP*F_reg_nt_coop))
    return (beware/k_deg)


mesh_size = 0.1
Gli = sp.arange(0.01, 29, mesh_size)
Gli_1 = sp.arange(0.01, 29, mesh_size)
fig = plt.figure()
ax = fig.add_subplot(1, 2, 1)
plt.title('Lai-Saha(Gli)-Gli')
ax.scatter(Gli_1, gli_curve(Gli)-Gli_1, s=0.5)
ax.grid(True, which='both', ls=':')
ax.axhline(y=0, color='k')

mesh_size = 0.01
Gli = sp.arange(0.01, 29.0, mesh_size)
Gli_1 = sp.arange(0.01, 29.0, mesh_size)

ax = fig.add_subplot(1, 2, 2)
plt.title('NewBeware(Gli)-Gli')
ax.scatter(Gli_1, gli_curve_1(Gli)-Gli_1, s=0.5)
ax.grid(True, which='both', ls=':')
ax.axhline(y=0, color='k')
plt.show()


def count_zeros(vector):
    count = 0
    for i in range(len(vector)-2):
        if vector[i] > 0 and vector[i+1] < 0:
            count += 1
        elif vector[i] < 0 and vector[i+1] > 0:
            count += 1
    return count


print(count_zeros(gli_curve(Gli)-Gli_1), count_zeros(gli_curve_1(Gli)-Gli_1))
print("--- %s seconds ---" % (time.time() - start_time))
