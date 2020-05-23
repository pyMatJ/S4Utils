#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 11:56:06 2020

@author: mathieu
"""

import numpy as np
import matplotlib.pyplot as plt
import S4 as S4

f = 1 ## single frequency in reduced units
a = 1 ## period
S = S4.New(Lattice=a, NumBasis=1) ## Simple lattice

SlabThick = 2 ## thickness of absorbing slab
AirThick = 1 ## 0 should work also...

n = 3 ## slab index, real part
k = 0.01 ## slab index, imaginary part
epsSlab = (n+1.0j*k)**2 ## slab permittivity
epsMirr = -1e10 ## mirror, to mimic a metal

theta = np.arange(0,89,2)

S.SetMaterial(Name='Air', Epsilon=(1.0+1.0j*0))
S.SetMaterial(Name='Slab', Epsilon=epsSlab)
# S.SetMaterial(Name='Mirror', Epsilon=epsMirr)

S.AddLayer(Name='AirTop', Thickness=AirThick, Material='Air')
S.AddLayer(Name='Slab', Thickness=SlabThick, Material='Slab')
# either mirror or air for substrate
# S.AddLayer(Name='AirBottom', Thickness=AirThick, Material='Mirror')
S.AddLayer(Name='AirBottom', Thickness=AirThick, Material='Air')


#%%
Rtm = np.empty(len(theta)) ## reflectivity
Ttm = np.zeros_like(Rtm) ## transmission
Atm = np.zeros_like(Rtm) ## absorption


for ii, thi in enumerate(theta):
    S.SetFrequency(f)
    S.SetExcitationPlanewave(
            IncidenceAngles=(thi,0), #S4 names are reversed (phi,theta)
            sAmplitude=0.,
            pAmplitude=1.,
            Order=0)
    inc, r = S.GetPowerFlux('AirTop', 0.)
    fw, _ = S.GetPowerFlux('AirBottom', 0.)
    Rtm[ii] = np.abs(-r/inc)
    Ttm[ii] = np.abs(fw/inc)
    fw1, bw1 = S.GetPowerFlux('Slab', 0)
    fw2, bw2 = S.GetPowerFlux('Slab', SlabThick)
    Atm[ii] = np.abs((fw2-fw1-(bw1-bw2))/inc)
        
#%%
chsum = Rtm+Ttm+Atm
figsp = plt.figure()
ax = figsp.add_subplot(111)
ax.plot(theta, Rtm, 'b', label='R')
ax.plot(theta, Ttm, 'r', label='T')
ax.plot(theta, Atm, 'o', c='grey', fillstyle='none', label='A')
ax.plot(theta, 1-(Rtm+Ttm), '--', label='1-(R+T)')
ax.plot(theta, chsum, ':k', label='sum')
ax.legend(bbox_to_anchor=(1,1))
ax.set_xlabel('Angle (degree)')
ax.set_ylabel('Reflection / Transmission / Absorption')
figsp.tight_layout()
plt.show()