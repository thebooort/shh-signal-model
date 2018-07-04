#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  1 21:38:51 2018

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
k_Gli = 90 # dissociation constant of activators for Gene enhancers
k_Gli3 = 90  # dissociation constant of activators for Gene enhancers
k_Gli3R = 90  # dissociation constant of repressors for Gene enhancers
k_RNAP = 1  # RNA polymerase binding affinity
RNAP = 1  # RNA polymerase concentration
c_b = 0.26  # BEWARE constant
c_b1 = 3.15
# from Lai-Schaffer classic model

Shh = 15  # Shh quantity [0,30]
k_shh = 0.58  # dissociation constant shh-ptc bindings [0.58,2.0]
k_ptc = 8.3*10**-2  # 1/2maximal concentration of ptc which inhibits smo signlaing
k_deg = 0.009  # degradation constant for all Gli related proteins
k_g3rc = 0.012  # rate constant for the conversion to signal strenGh
r_g3b = 1.6*10**-1  # basal rate of Gli3 synthesis
K_g3rc = 0.1  # sensitivity constant of the conversion to signal strenGh
k_deg_p = 0.09  # degradation rate constant for Ptc [0.045,0.071]
# --------------
K1 = 8.3*10**-1
K2 = 8.3*10**-1

e = 0.5
r = 0.2
v_max = 2.4*10**-1
r_bas = v_max/100
v_maxp = 7.5*10**-1
r_basp = v_maxp/100


def Signal(Ptc):

    return (1+Shh/k_shh)/(1+Shh/k_shh+Ptc/k_ptc)


def Promoter(Gli, Gli3, Gli3R):
    return ((Gli3*K1+Gli*K2)*(3*e**2*K1**2*K2**2+3*c*e*K1*K2*(Gli3*K1+Gli*K2+2*e*Gli3R*K1*r)+c**2*(Gli3**2*K1**2+Gli**2*K2**2+3*e*Gli*Gli3R*K1*K2*r + 3*e**2*Gli3R**2*K1**2*r**2 + Gli3*K1*(2*Gli*K2 + 3*e*Gli3R*K1*r))))/(3*c*K1*K2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**2 + c**2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**3 +K1**2*K2**2*(3*Gli3*K1 + 3*Gli3R*K1 + (3*Gli + K1)*K2))


def Basal(Gli, Gli3, Gli3R):
    return (3*c*K1*K2*(Gli3*K1+ Gli*K2 + Gli3R*K1*r)**2 + c**2*(Gli3*K1 + Gli*K2 + Gli3R*K1*r)**3 + K1**2*K2**2*(3*Gli3*K1 + 3*Gli*K2 + K1*(K2+ 3*Gli3R*r)))/ (3*c*K1*K2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**2 + c**2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**3 + K1**2*K2**2*(3*Gli3*K1 + 3*Gli3R*K1 + (3*Gli +K1)*K2))


def lai_saha_model(X, t):
    Gli, Gli3, Gli3R, Ptc = X

    dGli_dt = v_max*Promoter(Gli, Gli3, Gli3R)+r_bas*Basal(Gli, Gli3, Gli3R)-k_deg*Gli
    dGli3_dt = r_g3b/Gli-Gli3*(k_deg+k_g3rc/(K_g3rc+Signal(Ptc)))
    dGli3R_dt = Gli3*(k_g3rc/(K_g3rc+Signal(Ptc)))-k_deg*Gli3R
    dPtc_dt = v_maxp*Promoter(Gli, Gli3, Gli3R)+r_basp*Basal(Gli, Gli3, Gli3R)-k_deg_p*Ptc

    return dGli_dt, dGli3_dt, dGli3R_dt, dPtc_dt


# Frist we define our temporal range
t = sp.arange(0.0, 1200.0, 0.1)

# definition of odeint for solve the system numerically

vector_solution = odeint(lai_saha_model, [0.01, 0, 0, 0], t)

# Extraction of Gli,gli3,gli3r,ptc numerical values of the solution

evol_gli = vector_solution[:, 0]
evol_ptc = vector_solution[:, 3]
evol_gli3 = vector_solution[:, 1]
evol_gli3r = vector_solution[:, 2]
# regulation function with non/total cooperativity


def F_reg_nt_coop(Gli, Gli3, Gli3R):

    return (1-1/c+1/c*(1+a_Gli*c*Gli/k_Gli+a_Gli3*c*Gli3/k_Gli3+r_Gli3R*c*Gli3R/k_Gli3R)**3)/(1-1/c+1/c*(1+c*Gli/k_Gli+c*Gli3/k_Gli3+c*Gli3R/k_Gli3R)**3)


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
t = sp.arange(0.0, 1200.0, 0.1)

# definition of odeint for solve the system numerically

vector_solution = odeint(shh_evolution_system, [0.01, 0, 0, 0], t)

# Extraction of Gli,gli3,gli3r,ptc numerical values of the solution

evol_gli_1 = vector_solution[:, 0]
evol_ptc_1 = vector_solution[:, 3]
evol_gli3_1 = vector_solution[:, 1]
evol_gli3r_1 = vector_solution[:, 2]

# Plotting the results (scaling them previously)

fig, ax = plt.subplots()
ax.plot(t, evol_gli, label=r'Gli Lai-Saha',color='purple'  )
ax.plot(t, evol_gli_1, label=r'Gli New Model',color='orange')
plt.title('Gli Lai-Saha vs. Gli New Model')
scale_x = 60
ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_x))
ax.xaxis.set_major_formatter(ticks_x)
ax.set_xlabel(r"$time(hr)$")
ax.set_ylabel(r'Concentration quantity [nM]')
ax.hlines(y=evol_gli_1[-1], xmin=0, xmax=len(evol_gli_1)/10, linewidth=1.5 ,color='orange', linestyles='dotted', label=str(evol_gli_1[-1]))
ax.hlines(y=evol_gli[-1], xmin=0, xmax=len(evol_gli)/10, linewidth=1.5 ,color='purple', linestyles='dotted', label=str(evol_gli[-1]))
ax.legend(loc='best', fancybox=True, framealpha=0.5)
plt.show()

fig, ax = plt.subplots()
ax.plot(t, evol_ptc, label=r'Ptc Lai-Saha',color='purple'  )
ax.plot(t, evol_ptc_1, label=r'Ptc New Model',color='orange' )
plt.title('Ptc Lai-Saha vs. Ptc New Model')
scale_x = 60
ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_x))
ax.xaxis.set_major_formatter(ticks_x)
ax.set_xlabel(r"$time(hr)$")
ax.set_ylabel(r'Concentration quantity [nM]')
ax.hlines(y=evol_ptc_1[-1], xmin=0, xmax=len(evol_ptc_1)/10, linewidth=1.5 ,color='orange', linestyles='dotted', label=str(evol_ptc_1[-1]))
ax.hlines(y=evol_ptc[-1], xmin=0, xmax=len(evol_ptc)/10, linewidth=1.5 ,color='purple', linestyles='dotted', label=str(evol_ptc[-1]))
ax.legend(loc='best', fancybox=True, framealpha=0.5)
plt.show()

fig, ax = plt.subplots()
ax.plot(t, evol_gli3, label=r'Gli3 Lai-Saha',color='purple'  )
ax.plot(t, evol_gli3_1, label=r'Gli3 New Model',color='orange' )
plt.title('Gli3 Lai-Saha vs. Gli3 New Model')
scale_x = 60
ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_x))
ax.xaxis.set_major_formatter(ticks_x)
ax.set_xlabel(r"$time(hr)$")
ax.set_ylabel(r'Concentration quantity [nM]')
ax.hlines(y=evol_gli3_1[-1], xmin=0, xmax=len(evol_gli3_1)/10, linewidth=1.5 ,color='orange', linestyles='dotted', label=str(evol_gli3_1[-1]))
ax.hlines(y=evol_gli3[-1], xmin=0, xmax=len(evol_gli3)/10, linewidth=1.5 ,color='purple', linestyles='dotted', label=str(evol_gli3[-1]))
ax.legend(loc='best', fancybox=True, framealpha=0.5)
plt.show()

fig, ax = plt.subplots()
ax.plot(t, evol_gli3r, label=r'Gli3R Lai-Saha',color='purple'  )
ax.plot(t, evol_gli3r_1, label=r'Gli3R New Model',color='orange' )
plt.title('Gli3R Lai-Saha vs. Gli3R New Model')
scale_x = 60
ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/scale_x))
ax.xaxis.set_major_formatter(ticks_x)
ax.set_xlabel(r"$time(hr)$")
ax.set_ylabel(r'Concentration quantity [nM]')
ax.hlines(y=evol_gli3r_1[-1], xmin=0, xmax=len(evol_gli3r_1)/10, linewidth=1.5 ,color='orange', linestyles='dotted', label=str(evol_gli3r_1[-1]))
ax.hlines(y=evol_gli3r[-1], xmin=0, xmax=len(evol_gli3r)/10, linewidth=1.5 ,color='purple', linestyles='dotted', label=str(evol_gli3r[-1]))
ax.legend(loc='best', fancybox=True, framealpha=0.5)
plt.show()

print(evol_gli_1[-1],evol_gli[-1],evol_gli_1[-1]-evol_gli[-1])
print(evol_ptc_1[-1],evol_ptc[-1],evol_ptc_1[-1]-evol_ptc[-1])
print(evol_gli3_1[-1],evol_gli3[-1],evol_gli3_1[-1]-evol_gli3[-1])
print(evol_gli3r_1[-1],evol_gli3r[-1],evol_gli3r_1[-1]-evol_gli3r[-1])
