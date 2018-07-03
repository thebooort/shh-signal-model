#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 09:52:05 2018

@author: booort
"""


from math import *
import scipy as sp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# Definition of constants

# from beware new model:

c = 1  # positive constant, Greater than 1 implies cooperativity, less than 1 anti-cooperativity
a_Gli = 4.35  # transcriptional activation intensity for Gli
a_Gli3 = 4.35  # transcriptional activation intensity for Gli
r_Gli3R = 5*10**-5  # transcriptional repression intensity for Gli
k_Gli = 9*10**1  # dissociation constant of activators for Gene enhancers
k_Gli3 = 9*10**1   # dissociation constant of activators for Gene enhancers
k_Gli3R = 9*10**1   # dissociation constant of repressors for Gene enhancers
k_RNAP = 1  # RNA polymerase binding affinity
RNAP = 1  # RNA polymerase concentration
c_b = 1  # BEWARE constant

# from Lai-Schaffer classic model

Shh = 0  # Shh quantity [0,30]
k_shh = 0.58  # dissociation constant shh-Ptc bindings [0.58,2.0]
k_Ptc = 8.3*10**-2  # half maximal concentration of Ptc which inhibits smo signlaing
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
# regulation function with non/total cooperativity

def F_reg_nt_coop(Gli, Gli3, Gli3R):

    return (1-1/c+1/c*(1+a_Gli*c*Gli/k_Gli+a_Gli3*c*Gli3/k_Gli3+r_Gli3R*c*Gli3R/k_Gli3R)**3)/(1-1/c+1/c*(1+c*Gli/k_Gli+c*Gli3/k_Gli3+c*Gli3R/k_Gli3R)**3)

# Signal function


def Signal(Ptc):

    return (1+Shh/k_shh)/(1+Shh/k_shh+Ptc/k_Ptc)

# Beware operator


def BEWARE(Gli, Gli3, Gli3R):

    return c_b/(1+k_RNAP/(RNAP*F_reg_nt_coop(Gli, Gli3, Gli3R)))

def Promoter(Gli, Gli3, Gli3R):
    return ((Gli3*K1+Gli*K2)*(3*e**2*K1**2*K2**2+3*c*e*K1*K2*(Gli3*K1+Gli*K2+2*e*Gli3R*K1*r)+c**2*(Gli3**2*K1**2+Gli**2*K2**2+3*e*Gli*Gli3R*K1*K2*r + 3*e**2*Gli3R**2*K1**2*r**2 + Gli3*K1*(2*Gli*K2 + 3*e*Gli3R*K1*r))))/(3*c*K1*K2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**2 + c**2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**3 +K1**2*K2**2*(3*Gli3*K1 + 3*Gli3R*K1 + (3*Gli + K1)*K2))


def Basal(Gli, Gli3, Gli3R):
    return (3*c*K1*K2*(Gli3*K1+ Gli*K2 + Gli3R*K1*r)**2 + c**2*(Gli3*K1 + Gli*K2 + Gli3R*K1*r)**3 + K1**2*K2**2*(3*Gli3*K1 + 3*Gli*K2 + K1*(K2+ 3*Gli3R*r)))/ (3*c*K1*K2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**2 + c**2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**3 + K1**2*K2**2*(3*Gli3*K1 + 3*Gli3R*K1 + (3*Gli +K1)*K2))


# Frist we define our range and values
Gli = sp.arange(0.0, 30.0, 0.1)
Gli3 = 0
Gli3R_values = [0, 1, 2, 5, 10]

# Plotting configuration
fig, ax = plt.subplots()
ax.set_ylabel(r"$BEWARE_{value}$")
ax.set_xlabel(r'$Gli[nM]$')
plt.title(r'Variacion de BEWARE con [Gli3]=0')
for Gli3R in Gli3R_values:
    ax.plot(Gli, BEWARE(Gli, Gli3, Gli3R), label='G3R= '+str(Gli3R)+' nM')
ax.hlines(y=BEWARE(len(Gli)/10-1, 0, 0), xmin=0, xmax=len(Gli)/10, linewidth=1.5, color='grey', linestyles='dotted', label=str(BEWARE(len(Gli)-1, 0, 0)))
ax.legend(loc='lower right', fancybox=True, framealpha=0.5)
plt.show()

# Plotting configuration
fig, ax = plt.subplots()
ax.set_ylabel(r"$Promoter_{value}$")
ax.set_xlabel(r'$Gli[nM]$')
plt.title(r'Variacion de Promoter con [Gli3]=0')
for Gli3R in Gli3R_values:
    ax.plot(Gli, Promoter(Gli, Gli3, Gli3R), label='G3R= '+str(Gli3R)+' nM')
ax.hlines(y=BEWARE(len(Gli)/10-1, 0, 0), xmin=0, xmax=len(Gli)/10, linewidth=1.5, color='grey', linestyles='dotted', label=str(BEWARE(len(Gli)-1, 0, 0)))
ax.legend(loc='lower right', fancybox=True, framealpha=0.5)
plt.show()

fig, ax = plt.subplots()
ax.set_ylabel(r"$Basal_{value}$")
ax.set_xlabel(r'$Gli[nM]$')
plt.title(r'Variacion de Basal con [Gli3]=0')
for Gli3R in Gli3R_values:
    ax.plot(Gli, Basal(Gli, Gli3, Gli3R), label='G3R= '+str(Gli3R)+' nM')
ax.hlines(y=BEWARE(len(Gli)/10-1, 0, 0), xmin=0, xmax=len(Gli)/10, linewidth=1.5, color='grey', linestyles='dotted', label=str(BEWARE(len(Gli)-1, 0, 0)))
ax.legend(loc='lower right', fancybox=True, framealpha=0.5)
plt.show()

