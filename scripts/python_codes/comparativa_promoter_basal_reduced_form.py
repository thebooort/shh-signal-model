#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 11:03:14 2018

@author: booort

This scripts runs a comparative between ol define beware operators and the new forms that are created in Cambon Sanchez.
Our aim is to prove that their behaviour is the same, plotting old and new at the same time

"""
from math import *
import scipy as sp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# Definition of constants

# from Lai-Schaffer classic model (same used in new model)

Shh = 0  # Shh quantity [0,30]
k_shh = 0.58  # dissociation constant shh-Ptc bindings [0.58,2.0]
k_Ptc = 8.3*10**-2  # half maximal concentration of Ptc which inhibits smo signlaing
k_deg = 0.009  # degradation constant for all Gli related proteins
k_g3rc = 0.012  # rate constant for the conversion to signal strenGh
r_g3b = 1.6*10**-1  # basal rate of Gli3 synthesis
K_g3rc = 0.1  # sensitivity constant of the conversion to signal strenGh
k_deg_p = 0.09  # degradation rate constant for Ptc [0.045,0.071]
# -------------- (Specifically used lai-saha model)
K1 = 8.3*10**-1
K2 = 8.3*10**-1
c = 1
e = 0.5
r = 0.2
v_max = 2.4*10**-1
r_bas = v_max/100
v_maxp = 7.5*10**-1
r_basp = v_maxp/100

#Promoter definition

def Promoter(Gli, Gli3, Gli3R):
    return ((Gli3*K1+Gli*K2)*(3*e**2*K1**2*K2**2+3*c*e*K1*K2*(Gli3*K1+Gli*K2+2*e*Gli3R*K1*r)+c**2*(Gli3**2*K1**2+Gli**2*K2**2+3*e*Gli*Gli3R*K1*K2*r + 3*e**2*Gli3R**2*K1**2*r**2 + Gli3*K1*(2*Gli*K2 + 3*e*Gli3R*K1*r))))/(3*c*K1*K2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**2 + c**2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**3 +K1**2*K2**2*(3*Gli3*K1 + 3*Gli3R*K1 + (3*Gli + K1)*K2))

#Basal operator definition

def Basal(Gli, Gli3, Gli3R):
    return (3*c*K1*K2*(Gli3*K1+ Gli*K2 + Gli3R*K1*r)**2 + c**2*(Gli3*K1 + Gli*K2 + Gli3R*K1*r)**3 + K1**2*K2**2*(3*Gli3*K1 + 3*Gli*K2 + K1*(K2+ 3*Gli3R*r)))/ (3*c*K1*K2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**2 + c**2*(Gli3*K1 + Gli3R*K1 + Gli*K2)**3 + K1**2*K2**2*(3*Gli3*K1 + 3*Gli3R*K1 + (3*Gli +K1)*K2))


#Reduced form of promoter definition
def Reduced_Promoter(Gli, Gli3, Gli3R):
    return (1/c*(e+c*Gli*K1**-1+c*Gli3*K2**-1+e*r*c*Gli3R*K2**-1)**3-1/c*(e+e*r*c*Gli3R*K2**-1)**3)/(1-1/c+1/c*(1+c*Gli*K1**-1+c*Gli3*K2**-1+c*Gli3R*K2**-1)**3)

#Reduced form of basal definition
def Reduced_Basal(Gli, Gli3, Gli3R):
    return (1-1/c+1/c*(1+c*Gli*K1**-1+c*Gli3*K2**-1+r*c*Gli3R*K2**-1)**3)/(1-1/c+1/c*(1+c*Gli*K1**-1+c*Gli3*K2**-1+c*Gli3R*K2**-1)**3)




# Frist we define our range and values
Gli = sp.arange(0.0, 30.0, 0.1)
Gli3 = 0
Gli3R_values = [0, 1, 2, 5, 10]
fig = plt.figure()


# Plotting configuration
ax=fig.add_subplot(3,2,1)
ax.set_ylabel(r"$Promoter_{value}$")
ax.set_xlabel(r'$Gli[nM]$')
plt.title(r'Variacion de Promoter con [Gli3]=0')
for Gli3R in Gli3R_values:
    ax.plot(Gli, Promoter(Gli, Gli3, Gli3R), label='G3R= '+str(Gli3R)+' nM')
    ax.plot(Gli, Reduced_Promoter(Gli, Gli3, Gli3R),'--', label='G3R= '+str(Gli3R)+' nM', color='black')
ax.legend(loc='lower right', fancybox=True, framealpha=0.5)
#plt.show()

ax=fig.add_subplot(3,2,2)
ax.set_ylabel(r"$Basal_{value}$")
ax.set_xlabel(r'$Gli[nM]$')
plt.title(r'Variacion de Basal con [Gli3]=0')
for Gli3R in Gli3R_values:
    ax.plot(Gli, Basal(Gli, Gli3, Gli3R), label='G3R= '+str(Gli3R)+' nM')
    ax.plot(Gli, Reduced_Basal(Gli, Gli3, Gli3R),'--', label='G3R= '+str(Gli3R)+' nM',color='black')
ax.legend(loc='lower right', fancybox=True, framealpha=0.5)
#plt.show()

# Frist we define our range and values
Gli3 = sp.arange(0.0, 30.0, 0.1)
Gli = 0
Gli3R_values = [0, 1, 2, 5, 10]

# Plotting configuration
ax=fig.add_subplot(3,2,3)
ax.set_ylabel(r"$Promoter_{value}$")
ax.set_xlabel(r'$Gli3[nM]$')
plt.title(r'Variacion de Promoter con [Gli]=0')
for Gli3R in Gli3R_values:
    ax.plot(Gli3, Promoter(Gli, Gli3, Gli3R), label='G3R= '+str(Gli3R)+' nM')
    ax.plot(Gli3, Reduced_Promoter(Gli, Gli3, Gli3R),'--', label='G3R= '+str(Gli3R)+' nM', color='black')
ax.legend(loc='lower right', fancybox=True, framealpha=0.5)
#plt.show()

ax=fig.add_subplot(3,2,4)
ax.set_ylabel(r"$Basal_{value}$")
ax.set_xlabel(r'$Gli3[nM]$')
plt.title(r'Variacion de Basal con [Gli]=0')
for Gli3R in Gli3R_values:
    ax.plot(Gli3, Basal(Gli, Gli3, Gli3R), label='G3R= '+str(Gli3R)+' nM')
    ax.plot(Gli3, Reduced_Basal(Gli, Gli3, Gli3R),'--', label='G3R= '+str(Gli3R)+' nM',color='black')
ax.legend(loc='lower right', fancybox=True, framealpha=0.5)
#plt.show()

# Frist we define our range and values
Gli3R = sp.arange(0.0, 30.0, 0.1)
Gli3 = 2
Gli_values = [0, 1, 2, 5, 10]
print(Promoter(10, 2, 30))
# Plotting configuration
ax=fig.add_subplot(3,2,5)
ax.set_ylabel(r"$Promoter_{value}$")
ax.set_xlabel(r'$Gli3R[nM]$')
plt.title(r'Variacion de Promoter con [Gli3]=0')
for Gli_1 in Gli_values:
    ax.plot(Gli3R, Promoter(Gli_1, Gli3, Gli3R), label='Gli= '+str(Gli_1)+' nM')
    ax.plot(Gli3R, Reduced_Promoter(Gli_1, Gli3, Gli3R),'--', label='Gli= '+str(Gli_1)+' nM', color='black')
ax.legend(loc='lower right', fancybox=True, framealpha=0.5)
#plt.show()

ax=fig.add_subplot(3,2,6)
ax.set_ylabel(r"$Basal_{value}$")
ax.set_xlabel(r'$Gli3R[nM]$')
plt.title(r'Variacion de Basal con [Gli3]=0')
for Gli_1 in Gli_values:
    ax.plot(Gli3R, Basal(Gli_1, Gli3, Gli3R), label='Gli= '+str(Gli_1)+' nM')
    ax.plot(Gli3R, Reduced_Basal(Gli_1, Gli3, Gli3R),'--', label='Gli= '+str(Gli_1)+' nM',color='black')
ax.legend(loc='lower right', fancybox=True, framealpha=0.5)
plt.show()
