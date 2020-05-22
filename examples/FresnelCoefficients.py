#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 13:16:08 2020

@author: mathieu
"""

import numpy as np
import S4 as S4
import matplotlib.pyplot as plt

c_const = 3e8

#%% 

a = 1 ## period, normalized units. This defines the length scale
lbda = 1e-6 ## say we compute at 1µm wavelength
f = c_const/lbda ## ferquency in SI units
f0 = f/c_const*1e-6 ## so the reduced frequency is f/c_const*a[SI units]

Dimension = 3 # 2 or 3

if Dimension == 2:
    ### 2 dimensional simulation, 1D lattice
    S = S4.New(Lattice = a, 
                 NumBasis = 1)
elif Dimension == 3:
    ### 3 dimensional simulation, 2D lattice
    S = S4.New(Lattice = ((a,0), (0,a)), 
                 NumBasis = 1)

theta = np.arange(0,89,2) ## angle of incidence in degree


n_superstrate = 1 ## incident medium index
eps_superstrate = n_superstrate**2 ## incicent medium permittivity

n_substrate = 2 ## substrate index
eps_substrate = n_substrate**2 ## substrate permittivity 

S.SetMaterial(Name='Air', Epsilon = eps_superstrate)
S.SetMaterial(Name='Substrate', Epsilon=eps_substrate)

AirThick = 1 ## in reduced units, same scale as a
SubThick = 1 ## in reduced units, same scale as a

S.AddLayer(Name='Air', Thickness=AirThick, Material='Air')
S.AddLayer(Name='Substrate', Thickness=SubThick, Material='Substrate')


######################
## Too technical for tutorial. See later
# S.SetOptions(
#     Verbosity=0, ## How much feedback you get from C++ code
#     DiscretizedEpsilon=False, ## True as soon as you do something complicated
#     DiscretizationResolution=8, # at least 8 if near field calculation involved
#     LanczosSmoothing = True, ## Mabe ?? especially if metals
#     SubpixelSmoothing=True, 
#     PolarizationDecomposition=True, # mandatory if metals. Improves convergence
#     PolarizationBasis='Jones') # Normal should help

#%%
Rte = np.zeros(len(theta))
Rtm = np.zeros(len(theta))
Tte = np.zeros(len(theta))
Ttm = np.zeros(len(theta))

S.SetFrequency(f0)

for ii, thi in enumerate(theta):
    ########
    ## TM excitation
    S.SetExcitationPlanewave(
            IncidenceAngles=(thi,0), #S4 names are reversed (phi,theta)
            sAmplitude=0.,
            pAmplitude=1.,
            Order=0)
    inc, r = S.GetPowerFlux(Layer='Air',zOffset=0)
    fw, _ = S.GetPowerFlux(Layer='Substrate', zOffset=SubThick)
    Rtm[ii] = np.abs(-r/inc)
    Ttm[ii] = np.abs(fw/inc)
    
    #########
    ## TE excitation
    S.SetExcitationPlanewave(
            IncidenceAngles=(thi,0), #S4 names are reversed (phi,theta)
            sAmplitude=1.,
            pAmplitude=0.,
            Order=0)
    inc, r = S.GetPowerFlux(Layer='Air',zOffset=0)
    fw, _ = S.GetPowerFlux(Layer='Substrate', zOffset=SubThick)
    Rte[ii] = np.abs(-r/inc)
    Tte[ii] = np.abs(fw/inc)
    
#%%
#####################
## Fresnel coefficients 
th1 = np.deg2rad(theta)
th2 = np.arcsin(n_superstrate/n_substrate*np.sin(th1))
Fresnel_Rte = ((n_superstrate*np.cos(th1)-n_substrate*np.cos(th2))/(n_superstrate*np.cos(th1)+n_substrate*np.cos(th2)))**2
Fresnel_Rtm = ((n_superstrate*np.cos(th2)-n_substrate*np.cos(th1))/(n_superstrate*np.cos(th2)+n_substrate*np.cos(th1)))**2
Fresnel_Tte = (n_substrate*np.cos(th2)/(n_superstrate*np.cos(th1)))*((2*n_superstrate*np.cos(th1))/(n_superstrate*np.cos(th1)+n_substrate*np.cos(th2)))**2
Fresnel_Ttm = (n_substrate*np.cos(th2)/(n_superstrate*np.cos(th1)))*((2*n_superstrate*np.cos(th1))/(n_superstrate*np.cos(th2)+n_substrate*np.cos(th1)))**2

th_brewster = np.rad2deg(np.arctan2(n_substrate, n_superstrate))
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(theta, Rte, 'ob', fillstyle='none', label='Rte')
ax.plot(theta, Tte, 'sb', fillstyle='none', label='Tte')
ax.plot(theta, Rtm, 'or', fillstyle='none', label='Rtm')
ax.plot(theta, Ttm, 'sr', fillstyle='none', label='Ttm')
ax.plot(theta, Fresnel_Rte, 'b')
ax.plot(theta, Fresnel_Tte, 'b')
ax.plot(theta, Fresnel_Rtm, 'r')
ax.plot(theta, Fresnel_Ttm, 'r')
ax.plot([th_brewster, th_brewster],[0,1], ':k', label='Brewster')
ax.set_ylim([0,1])
ax.set_xlabel('Angle (degree)')
ax.set_ylabel('Reflection/Transmission')
ax.legend(bbox_to_anchor=(1,1))
fig.tight_layout()
fig.show()

