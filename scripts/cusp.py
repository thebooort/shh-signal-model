#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  8 21:15:01 2018

@author: booort
"""

import scipy as sp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# Definition of constants

a=1
b=1
def f(Ptc,t):
    return a+b-Ptc**3
t = sp.arange(0.0, 20.0, 0.1)

# definition of odeint for solve the system numerically

vector_solution = odeint(f,0, t)
evol_gli = vector_solution[:, 0]
fig, ax = plt.subplots()
ax.plot(t, evol_gli, label=r'Gli' )