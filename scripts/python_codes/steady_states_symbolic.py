#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 10:41:46 2018

this script is created to solve the steady states

@author: booort
"""

from sympy import *
init_printing()
c = Symbol('c', real=True)  # positive constant, greater than 1 implies cooperativity, less than 1 anti-cooperativity
a_gli = Symbol('a_gli', real=True)   # transcriptional activation intensity for gli
a_gli3 = Symbol('a_gli3', real=True)   # transcriptional activation intensity for gli
r_gli3R = Symbol('r_gli3R', real=True)   # transcriptional repression intensity for gli
k_gli = Symbol('k_gli', real=True)   # dissociation constant of activators for gene enhancers
k_gli3 = Symbol('k_gli3', real=True)   # dissociation constant of activators for gene enhancers
k_gli3R = Symbol('k_gli3R', real=True)   # dissociation constant of repressors for gene enhancers
k_RNAP = Symbol('k_RNAP', real=True)   # RNA polymerase binding affinity
RNAP = Symbol('RNAP', real=True)   # RNA polymerase concentration
c_b = Symbol('c_b', real=True)   # BEWARE constant

# from Lai-Schaffer classic model

Shh = Symbol('Shh', real=True)   # Shh quantity [0,30]
k_shh = Symbol('k_shh', real=True)   # dissociation constant shh-ptc bindings [0.58,2.0]
k_ptc = Symbol('k_ptc', real=True)   # half maximal concentration of ptc which inhibits smo signlaing
k_deg = Symbol('k_deg', real=True)   # degradation constant for all gli related proteins
k_g3rc = Symbol('k_g3rc', real=True)   # rate constant for the conversion to signal strengh
r_g3b = Symbol('r_g3b', real=True)   # basal rate of Gli3 synthesis
K_g3rc = Symbol('K_g3rc', real=True)   # sensitivity constant of the conversion to signal strengh
k_deg_p = Symbol('k_deg_p', real=True)  # Degradation rate constant for Ptc [0.045,0.071]
 
gli = Symbol('gli', real=True)
gli3 = Symbol('gli3', real=True)
gli3R = Symbol('gli3R', real=True)
# Equation system

freg = (1-1/c+1/c*(1+a_gli*c*gli/k_gli+a_gli3*c*gli3/k_gli3R+r_gli3R*c*gli3R/k_gli3R)**3)/(1-1/c+1/c*(1+c*gli/k_gli+c*gli3/k_gli3R+c*gli3R/k_gli3R)**3)
beware = c_b/(1+k_RNAP/(RNAP*factor(expand(freg))))


print(latex(factor(expand(freg)),mode='plain'))
print(factor(expand(beware)))
