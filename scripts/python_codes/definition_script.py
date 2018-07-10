#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 17:48:12 2018

This script is created to define all the functions needed in my study

@author: booort
"""

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
c_b = 0.26  # BEWARE constant
c_b1 = 3.15
# from Lai-Schaffer classic model

Shh = 0.1  # Shh quantity [0,30]
k_shh = 0.58  # dissociation constant shh-ptc bindings [0.58,2.0]
k_ptc = 8.3*10**-2  # half maximal concentration of ptc which inhibits smo signlaing
k_deg = 0.009  # degradation constant for all Gli related proteins
k_g3rc = 0.012*1000  # rate constant for the conversion to signal strenGh
r_g3b = 1.6*10**-1  # basal rate of Gli3 synthesis
K_g3rc = 0.1  # sensitivity constant of the conversion to signal strenGh
k_deg_p = 0.09  # degradation rate constant for Ptc [0.045,0.071]


# regulation function with non/total cooperativity

def F_reg_nt_coop(Gli, Gli3, Gli3R):

    return (1-1/c+1/c*(1+a_Gli*c*Gli/k_Gli+a_Gli3*c*Gli3/k_Gli3+r_Gli3R*c*Gli3R/k_Gli3R)**3)/(1-1/c+1/c*(1+c*Gli/k_Gli+c*Gli3/k_Gli3+c*Gli3R/k_Gli3R)**3)


def Signal(Ptc):

    return (1+Shh/k_shh)/(1+Shh/k_shh+Ptc/k_ptc)


def BEWARE(Gli, Gli3, Gli3R):

    return c_b/(1+k_RNAP/(RNAP*F_reg_nt_coop(Gli, Gli3, Gli3R)))


def shh_evolution_system(X, t):
    Gli, Gli3, Gli3R, Ptc = X

    dGli_dt = BEWARE(Gli, Gli3, Gli3R)-k_deg*Gli
    dGli3_dt = r_g3b/Gli-Gli3*(k_deg+k_g3rc/(K_g3rc+Signal(Ptc)))
    dGli3R_dt = Gli3*(k_g3rc/(K_g3rc+Signal(Ptc)))-k_deg*Gli3R
    dPtc_dt = c_b1*BEWARE(Gli, Gli3, Gli3R)-k_deg_p*Ptc

    return dGli_dt, dGli3_dt, dGli3R_dt, dPtc_dt


# Frist we define our temporal range
t = sp.arange(0.0, 2400.0, 0.1)

# definition of odeint for solve the system numerically

vector_solution = odeint(shh_evolution_system, [0.01, 0, 0, 0], t)

# Extraction of Gli,gli3,gli3r,ptc numerical values of the solution

evol_gli_1 = vector_solution[:, 0]
evol_ptc_1 = vector_solution[:, 3]
evol_gli3_1 = vector_solution[:, 1]
evol_gli3r_1 = vector_solution[:, 2]

# Plotting the results (scaling them previously)
fig, ax = plt.subplots()
ax.plot(t, evol_gli_1, label=r'Gli' )
ax.plot(t, evol_ptc_1, label=r'Ptc' )
ax.plot(t, evol_gli3_1, label=r'Gli3' )
ax.plot(t, evol_gli3r_1, label=r'Gli3R' )
scale_x = 600
ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_x))
ax.xaxis.set_major_formatter(ticks_x)

ax.set_xlabel(r"$time(hr)$")
ax.set_ylabel(r'Concentration quantity [nM]')
ax.hlines(y=evol_gli_1[-1], xmin=0, xmax=len(evol_gli_1)/10, linewidth=1.5 ,color='blue', linestyles='dotted', label=str(evol_gli_1[-1]))
ax.hlines(y=evol_ptc_1[-1], xmin=0, xmax=len(evol_ptc_1)/10, linewidth=1.5 ,color='orange', linestyles='dotted', label=str(evol_ptc_1[-1]))
ax.hlines(y=evol_gli3_1[-1], xmin=0, xmax=len(evol_gli3_1)/10, linewidth=1.5 ,color='green', linestyles='dotted', label=str(evol_gli3_1[-1]))
ax.hlines(y=evol_gli3r_1[-1], xmin=0, xmax=len(evol_gli3r_1)/10, linewidth=1.5 ,color='red', linestyles='dotted', label=str(evol_gli3r_1[-1]))
ax.legend(loc='right', fancybox=True, framealpha=0.5)
plt.title('New Model')
ax.grid(True, which='both',ls=':')
plt.show()

print(evol_gli_1[2000])








