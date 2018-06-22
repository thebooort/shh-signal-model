#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 17:48:12 2018

This script is created to define all the functions needed in my study

@author: booort
"""

from math import *
import scipy as sp
from scipy.integrate import odeint

# Definition of constants

# from beware new model:

c = 1  # positive constant, greater than 1 implies cooperativity, less than 1 anti-cooperativity
a_gli = 1  # transcriptional activation intensity for gli
a_gli3 = 1  # transcriptional activation intensity for gli
r_gli3R = 1  # transcriptional repression intensity for gli
k_gli = 1  # dissociation constant of activators for gene enhancers
k_gli3 = 1  # dissociation constant of activators for gene enhancers
k_gli3R = 1  # dissociation constant of repressors for gene enhancers
k_RNAP = 1  # RNA polymerase binding affinity
RNAP = 1  # RNA polymerase concentration
c_b = 1  # BEWARE constant

# from Lai-Schaffer classic model

Shh = 1  # Shh quantity [0,30]
k_shh = 0.58  # dissociation constant shh-ptc bindings [0.58,2.0]
k_ptc = 8.3*10**-11  # half maximal concentration of ptc which inhibits smo signlaing
k_deg = 0.009  # degradation constant for all gli related proteins
k_g3rc = 0.012  # rate constant for the conversion to signal strengh
r_g3b = 1.6*10**-19  # basal rate of Gli3 synthesis
K_g3rc = 0.1  # sensitivity constant of the conversion to signal strengh
k_deg_p = 0.09 # Degradation rate constant for Ptc [0.045,0.071]
 



#Regulation function with non/total cooperativity

def F_reg_nt_coop(gli, gli3, gli3R):
    
    return (1-1/c+1/c*(1+a_gli*c*gli/k_gli+a_gli3*c*gli3/k_gli3R+r_gli3R*c*gli3R/k_gli3R)**3)/(1-1/c+1/c*(1+c*gli/k_gli+c*gli3/k_gli3R+c*gli3R/k_gli3R)**3)


def Signal(Ptc):
    
    return (1+Shh/k_shh)/(1+Shh/k_shh+Ptc/k_ptc)


def BEWARE(gli, gli3, gli3R):

    return c_b/(1+k_RNAP/(RNAP*F_reg_nt_coop(gli,gli3,gli3R)))


def shh_evolution_system(X,t):
    Gli, Gli3, Gli3R, Ptc = X
    
    dGli_dt = BEWARE(Gli, Gli3, Gli3R)-k_deg*Gli
    dGli3_dt = r_g3b/Ptc-Gli3*(k_deg+k_g3rc/(K_g3rc+Signal(Ptc)))
    dGli3R_dt = Gli3*(k_g3rc/(K_g3rc+Signal(Ptc)))-k_deg*Gli3R
    dPtc_dt = BEWARE(Gli, Gli3, Gli3R)-k_deg_p*Ptc
    
    return dGli_dt, dGli3_dt, dGli3R_dt, dPtc_dt

#Frist we define our temporal range
t = sp.arange(0.0, 200.0, 0.1)

# definition of odeint for solve the system numerically 

vector_solution = odeint(shh_evolution_system,[0,0,0,0], t)












