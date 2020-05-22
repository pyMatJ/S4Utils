#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 13:29:48 2020

@author: mathieu
"""
import sys, time
import S4 as S4
import numpy as np
import matplotlib.pyplot as plt
dirS4Utils = '/home/mathieu/Softwares/S4Utils/'
sys.path.append(dirS4Utils)
import MaterialFunctions as mat
import S4Utils as S4Utils

tstart = time.time()

c_const = 3e8
#%%
# General parameters
fmin = 12.0*1e12 ## 400cm-1 en Hz
fmax = 50.0*1e12
f = np.linspace(fmin, fmax, 200)
f0 = f/c_const*1e-6

px = 3.6
ff = 0.8
s = ff*px
S = S4.New(Lattice = px,
           NumBasis = 41) ### NumBasis <=> halfnpw in RCWA

ISBOn = False ## whether or not using a doped active region

theta = np.arange(0,90,5) ### for a dispersion plot
# theta = [0] ## for a single spectrum

epsAu = mat.epsAu(f)
epsGaAs = mat.epsGaAs(f)

#### doped region parameters
fisb = 31*1e12 ## isb frequency before plasma shift, Hz
omega_isb = 2*np.pi*fisb ## pulsation
gamma_isb = 0.1*fisb ## broadening
N2D = 7e11*1e4 ## doping
w_well=7.5e-9 ## well thickness
w_barr = 23.75e-9 ## barrier thickness 
nQW = 32 ## number of wells
fw = nQW*w_well/(nQW*(w_well+w_barr)) ## filling factor
eps_w = mat.epsGaAs(f) ## well background material
eps_b = mat.epsAlGaAs(f, xAl=0.25) ## barrier background material
omega_p = mat.omegaP_2D(N2D, 0.063, 10.89, w_well) ## plasma frequency 
epsARxx, epsARzz = np.conj(mat.epsZal(f, eps_w, eps_b, omega_isb, gamma_isb, omega_p, fw))
epsAR = np.array([[epsARxx, np.zeros(len(f)), np.zeros(len(f))],
                     [np.zeros(len(f)), epsARxx, np.zeros(len(f))],
                     [np.zeros(len(f)), np.zeros(len(f)), epsARzz]])


DisplayThick = 0.5 ## Incident medium thickness
ARThick = 1.0 ## Active region thickness
AuThick = 0.1 ## Gold thickness


### Materials
S.SetMaterial(Name='Air', Epsilon=(1.0 + 0.0*1.0j))
S.SetMaterial(Name='Au', Epsilon=(epsAu[0]))
# material list for the update function
Mat_list = ['Au', 'AR']
Eps_list = [epsAu]
if ISBOn:
    print('Doped Active Region in')
    S.SetMaterial(Name='AR', Epsilon=S4Utils.totuple(epsAR[:,:,0]))
    Eps_list.append(epsAR)
else:
    S.SetMaterial(Name='AR', Epsilon=(epsGaAs[0]))
    Eps_list.append(epsGaAs)
    
### Layers
S.AddLayer(Name='top', Thickness = DisplayThick, Material = 'Air')
S.AddLayer(Name='TopGrating', Thickness = AuThick, Material = 'Air')
S.AddLayer(Name='AR', Thickness = ARThick, Material = 'AR')
S.AddLayer(Name='Bulk', Thickness = 2*AuThick, Material = 'Au')

### Geometry
S.SetRegionRectangle(
    Layer='TopGrating',
    Material='Au',
    Center=(0.0,0.0),
    Angle=0.0,
    Halfwidths=(s/2.,0)) # 1D

### Simulation options
S.SetOptions(
    Verbosity=0, ## verbosity of the C++ code
    DiscretizedEpsilon=True, ## Necessary for high contrast 
    DiscretizationResolution=8,  ## at least 8 if near field calculations
    LanczosSmoothing = True, ## Mabe ?? especially for metals, for near fields
    SubpixelSmoothing=True, ## definitely
    PolarizationDecomposition=True, ## Along with 'normal', should help convergence
    PolarizationBasis='Normal')
#%%
R = np.empty((len(theta),len(f)))

ProgOrder = 0
for ii, thi in enumerate(theta):
    currOrder = int((ii*10)/len(theta))
    if currOrder>ProgOrder:
        print(currOrder)
        ProgOrder = currOrder
    ProgF = 0
    for jj, fj in enumerate(f0):
        currF = int((jj*10)/len(f0))
        if currF>ProgF:
            print('\t %d'%currF)
            ProgF = currF
        
        S.SetFrequency(fj)
        S4Utils.UpdateMaterials(S, Mat_list, Eps_list, fj, f0)
        S.SetExcitationPlanewave(
                IncidenceAngles=(thi,0.),
                sAmplitude=0.,
                pAmplitude=1.,
                Order=0)
        inc, r = S.GetPowerFlux('top', 0.)
        R[ii,jj] = np.abs(-r/inc)
    
#%%
if len(theta)==1:
    figsp = plt.figure()    
    ax = figsp.add_subplot(111)
    ax.plot(f, R.T)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Reflectivity')
else:
    figdisp = plt.figure()
    axm = figdisp.add_subplot(111)
    thm, fm = np.meshgrid(theta,f)
    km = (px*1e-6)*(2*np.pi*fm)/c_const*np.sin(np.deg2rad(thm))/np.pi
    cax = axm.pcolormesh(km, fm, R.T,
                         vmin=0, vmax=1,
                         shading='None')
    axm.set_xlabel('pk$_{\parallel}$/$\pi$')
    axm.set_ylabel('Frequency (cm$^{-1}$)')
    figdisp.tight_layout()
plt.show()

#%%
tstop = time.time()
tscript = tstop-tstart
print('Execution time %.2fs'%(tscript))