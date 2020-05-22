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

Dimension = 3 # or 3
if Dimension == 2:
    ### 2 dimensional simulation, 1D lattice
    S = S4.New(Lattice = a, 
                 NumBasis = 1)
elif Dimension == 3:
    ### 3 dimensional simulation, 2D lattice
    S = S4.New(Lattice = ((a,0), (0,a)), 
                 NumBasis = 1)


n_superstrate = 1 ## incident medium index
eps_superstrate = n_superstrate**2 ## incicent medium permittivity

n_slab = 2 ## substrate index
eps_slab = n_slab**2 ## substrate permittivity 

S.SetMaterial(Name='Air', Epsilon = 1)
S.SetMaterial(Name='Slab', Epsilon=eps_slab)

AirThick = 1
SlabThick = 1

S.AddLayer(Name='AirTop', Thickness=AirThick, Material='Air')
S.AddLayer(Name='Slab', Thickness=SlabThick, Material='Slab')
S.AddLayer(Name='AirBottom', Thickness=AirThick, Material='Air')


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

#%%  Sweep over angle of incidence, fixed wavelength
lbda = 1e-6 ## say we compute at 1µm wavelength
f = c_const/lbda ## ferquency in SI units
f0 = f/c_const*1e-6 ## so the reduced frequency is f/c_const*a[SI units]
theta = np.arange(0,89,1)

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
    inc, r = S.GetPowerFlux(Layer='AirTop',zOffset=0)
    fw, _ = S.GetPowerFlux(Layer='AirBottom', zOffset=AirThick)
    Rtm[ii] = np.abs(-r/inc)
    Ttm[ii] = np.abs(fw/inc)
    
    #########
    ## TE excitation
    S.SetExcitationPlanewave(
            IncidenceAngles=(thi,0), #S4 names are reversed (phi,theta)
            sAmplitude=1.,
            pAmplitude=0.,
            Order=0)
    inc, r = S.GetPowerFlux(Layer='AirTop',zOffset=0)
    fw, _ = S.GetPowerFlux(Layer='AirBottom', zOffset=AirThick)
    Rte[ii] = np.abs(-r/inc)
    Tte[ii] = np.abs(fw/inc)
    
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(theta, Rte, 'o-b', fillstyle='none', label='Rte')
ax.plot(theta, Tte, 's-b', fillstyle='none', label='Tte')
ax.plot(theta, Rtm, 'o-r', fillstyle='none', label='Rtm')
ax.plot(theta, Ttm, 's-r', fillstyle='none', label='Ttm')
ax.set_ylim([0,1])
ax.set_xlabel('Angle (degree)')
ax.set_ylabel('Reflection/Transmission')
ax.legend(bbox_to_anchor=(1,1))
fig.tight_layout()
fig.show()

#%%  Sweep over frequency, fixed angle of incidence
lbda = np.linspace(400e-9,1e-6,200) ## 400nm to 1µm
f = c_const/lbda ## ferquency in SI units
f0 = f/c_const*1e-6 ## so the reduced frequency is f/c_const*a[SI units]
theta = 15 ## fixed angle of incidence

Rte = np.zeros(len(f0))
Rtm = np.zeros(len(f0))
Tte = np.zeros(len(f0))
Ttm = np.zeros(len(f0))


for ii, fi in enumerate(f0):
    S.SetFrequency(fi)
    ########
    ## TM excitation
    S.SetExcitationPlanewave(
            IncidenceAngles=(theta,0), #S4 names are reversed (phi,theta)
            sAmplitude=0.,
            pAmplitude=1.,
            Order=0)
    inc, r = S.GetPowerFlux(Layer='AirTop',zOffset=0)
    fw, _ = S.GetPowerFlux(Layer='AirBottom', zOffset=AirThick)
    Rtm[ii] = np.abs(-r/inc)
    Ttm[ii] = np.abs(fw/inc)
    
    #########
    ## TE excitation
    S.SetExcitationPlanewave(
            IncidenceAngles=(theta,0), #S4 names are reversed (phi,theta)
            sAmplitude=1.,
            pAmplitude=0.,
            Order=0)
    inc, r = S.GetPowerFlux(Layer='AirTop',zOffset=0)
    fw, _ = S.GetPowerFlux(Layer='AirBottom', zOffset=AirThick)
    Rte[ii] = np.abs(-r/inc)
    Tte[ii] = np.abs(fw/inc)
    
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(lbda, Rte, 'o-b', fillstyle='none', label='Rte')
ax.plot(lbda, Tte, 's-b', fillstyle='none', label='Tte')
ax.plot(lbda, Rtm, 'o-r', fillstyle='none', label='Rtm')
ax.plot(lbda, Ttm, 's-r', fillstyle='none', label='Ttm')
ax.set_ylim([0,1])
ax.set_xlabel('Wavelength (m)')
ax.set_ylabel('Reflection/Transmission')
ax.legend(bbox_to_anchor=(1,1))
fig.tight_layout()
fig.show()

