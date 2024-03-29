#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 13:29:48 2020

@author: mathieu
"""
import time ## to time script execution
import S4 as S4
import numpy as np
import matplotlib.pyplot as plt
import S4Utils.S4Utils as S4Utils
import S4Utils.MaterialFunctions as mat

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
NBasis = 41
S = S4.New(Lattice = px,
           NumBasis = NBasis) ### NumBasis <=> halfnpw in RCWA

theta = np.arange(0,90,5) ### for a dispersion plot
# theta = [0] ## for a single spectrum

epsAu = mat.epsAu(f)
epsGaAs = mat.epsGaAs(f)


DisplayThick = 0.5 ## Incident medium thickness
ARThick = 1.0 ## Active region thickness
AuThick = 0.1 ## Gold thickness


### Materials
S.SetMaterial(Name='Air', Epsilon=(1.0 + 0.0*1.0j))
S.SetMaterial(Name='Au', Epsilon=(epsAu[0]))
S.SetMaterial(Name='AR', Epsilon=(epsGaAs[0]))
    # material list for the update function
Mat_list = ['Au', 'AR']
Eps_list = [epsAu, epsGaAs]
    
### Layers
S.AddLayer(Name='top', Thickness = DisplayThick, Material = 'Air') ## incident medium, air
S.AddLayer(Name='TopGrating', Thickness = AuThick, Material = 'Air') ## grating layer, will be patterned
S.AddLayer(Name='AR', Thickness = ARThick, Material = 'AR') ## active region layer 
S.AddLayer(Name='Bulk', Thickness = 2*AuThick, Material = 'Au') ## bottom mirror
TotalThick = DisplayThick+ARThick+3*AuThick ## total thickness of the simulation

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

ProgAngle = 0 ## progress in angle sweep
for ii, thi in enumerate(theta): ## angle sweep
    currAngle = int((ii*10)/len(theta)) ## current angle by 10% steps
    if currAngle>ProgAngle:  
        print(currAngle) # print progress every 10%
        ProgAngle = currAngle
    ProgF = 0  ## progress in frequency sweep
    for jj, fj in enumerate(f0):  ## frequency sweep 
        currF = int((jj*10)/len(f0)) ## current frequency by 10% steps
        if currF>ProgF:
            print('\t %d'%currF) # print progress every 10%
            ProgF = currF
        
        S.SetFrequency(fj)  # set the current frequency 
        S4Utils.UpdateMaterials(S, Mat_list, Eps_list, fj, f0) # set epsilons
        S.SetExcitationPlanewave(
                IncidenceAngles=(thi,0.),
                sAmplitude=0.,
                pAmplitude=1.,  ## p-pol plane wave
                Order=0)
        inc, r = S.GetPowerFlux('top', 0.)
        R[ii,jj] = np.abs(-r/inc) ## reflectivity
    
#%% Plotting the results
if len(theta)==1:
    figsp = plt.figure()    
    ax = figsp.add_subplot(111)
    ax.plot(f, R.T, label='R')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Reflectivity')
    ax.set_title('Theta=%d, NumBasis=%d'%(theta[0], NBasis))
    ax.legend(bbox_to_anchor=(1,1))
    figsp.tight_layout()
else:
    figdisp = plt.figure()
    axm = figdisp.add_subplot(111)
    thm, fm = np.meshgrid(theta,f)
    km = (px*1e-6)*(2*np.pi*fm)/c_const*np.sin(np.deg2rad(thm))/np.pi
    cax = axm.pcolormesh(km, fm, R.T,
                         vmin=0, vmax=1,
                         shading='None')
    axm.set_xlabel('pk$_{\parallel}$/$\pi$')
    axm.set_ylabel('Frequency (Hz)')
    axm.set_title('NumBasis=%d'%NBasis)
    figdisp.colorbar(cax,label='Reflectivity')
    figdisp.tight_layout()
plt.show()

#%% timing
tstop = time.time()
tscript = tstop-tstart
print('Execution time %.2fs'%(tscript))

#%% Plotting electric field
########
## Setting up the simulation
fplot = 28.8e12 ## frequency at which we plot (SI)
f0plot = fplot/c_const*1e-6 # reduced units
S.SetFrequency(f0plot)
S4Utils.UpdateMaterials(S, Mat_list, Eps_list, f0plot, f0) # update materials
S.SetExcitationPlanewave(
    IncidenceAngles=(0.,0.), ## normal incidence
    sAmplitude=0.,
    pAmplitude=1.,  ## p-pol plane wave
    Order=0)
########
## Spatial coordinates for the plot
resx = 150 ## x-resolution
resz = 150 ## z-resolution
x = np.linspace(-px/2, px/2, resx) ## x coordinates
z = np.linspace(0, TotalThick, resz) ## z coordinates

### compute electric and magnetic field
E, H = S4Utils.GetSlice(S, ax1=x, ax2=z, plane='xz', mode='Field')
FigEz = S4Utils.SlicePlot(x, z ,np.real(E[:,:,2]), hcoord=s/4)
FigEz.ax_2D.set_ylabel('z')
FigEz.ax_2D.set_xlabel('x')
FigEz.axcbar.set_ylabel('Ez')

## retrieve the recomposed permittivity 
Eps_plot = S4Utils.GetSlice(S, ax1=x, ax2=z, plane='xz', mode='Epsilon')
FigEps = S4Utils.SlicePlot(x, z, np.real(Eps_plot), cmap=plt.cm.bone,
                           sym=False)
FigEps.ax_2D.set_ylabel('z')
FigEps.ax_2D.set_xlabel('x')
FigEps.ax_a1slice.set_xlabel(r'$\varepsilon$')
FigEps.ax_a2slice.set_ylabel(r'$\varepsilon$')
FigEps.axcbar.set_ylabel(r'$\varepsilon$')
