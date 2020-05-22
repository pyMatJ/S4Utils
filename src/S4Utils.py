#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 12:05:13 2020

@author: mathieu
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pyvista as pv

def totuple(a):
    """
    Converts lists and array to a tuple form. Usefull since S4 mostly accepts 
    tuples instead of numpy arrays.
    """
    try:
        return tuple(totuple(i) for i in a)
    except TypeError:
        return a
    
def UpdateMaterials(S, Mat_list, Eps_list, fi, f0):
    """
    Updates a Simulation to set all refractive indices at given frequency
    
    Parameters
    ----------
    S: S4.Simulation
    Mat_list: list
        list of string containing the name of the materials to update. These 
        are typically the ones passes in the 'Material' argument for the 
        S4.Simulation.SetMaterial() function
    Eps_list: list
        list of arrays containing the permittivities.
        each element is a 1D array or a 3x3xlen(f) array
    fi: float
        frequency at which the simulation is set (reduced units)
    f0: 1Darray
        array of reduced frequencies over which the simulation runs
    """
    fidx = np.abs(f0-fi).argmin()
    for ii, mat_i in enumerate(Mat_list):
        if Eps_list[ii].shape==f0.shape:
            epsi = Eps_list[ii][fidx]
        elif Eps_list[ii].shape==(3,3,len(f0)):
            epsi = totuple(Eps_list[ii][:,:,fidx])
        elif Eps_list[ii].shape==(len(f0),3,3):
            epsi = totuple(Eps_list[ii][fidx])
        else:
            print('Wrong shape for the permittivities of '+mat_i)
        S.SetMaterial(Name=mat_i, Epsilon=epsi)
        
        

def GetField_xy(S, res, zs):
    """
    Get E and H fields in an xy slice. Built-in S4 function that 
    requires some care and a tweak, but extra efficient !!
    
    Parameters
    -----------
    S: S4.Simulation
    res: tuple
        (resx, resy) tuple of number of points along x and y
        Note that it corresponds to *resy lines* and *resx columns*
    zs: float
        z-coordinate at which to slice the field
    """
    resx, resy = res[0], res[1] 
    
    Eplane, Hplane = S.GetFieldsOnGrid(z=zs, NumSamples=(resx,resy), # resy lines x resx columns !
                                       Format='Array')
    Eplane = np.array(Eplane)
    Hplane = np.array(Hplane)
    ## These are actually off by 1/2 grid on each direction 
    ## this step rolls them back to center
    # **not corrected**  still off-centered by a small amount (probably mistake in C++) 
    Eplane = np.roll(np.roll(Eplane, Eplane.shape[0]//2,axis=0),Eplane.shape[1]//2, axis=1)
    Hplane = np.roll(np.roll(Hplane, Hplane.shape[0]//2,axis=0),Hplane.shape[1]//2, axis=1)
    return Eplane, Hplane

def GetField_xz(S, x, y, z, ys=0.0):
    """
    Get E and H fields in an xz slice.
    Wraps around GetField_xy which is hard-coded in the C code and thus much faster.
    *BUT* you should check that it matches with GetSlice at least once because 
    the C function GetFieldsOnGrid tends to be touchy.
    """
    E_slice = np.empty((len(z), len(x), 3), dtype=np.complex128)
    H_slice = np.empty((len(z), len(x), 3), dtype=np.complex128)
    xm_xz, zm_xz = np.meshgrid(x,z)
    res = (len(x), len(y))
    if len(y)==1:
        print('''Warning. Wrapping around GetFieldsOnGrid in a 2D 
              simulation might result in false output enforcing symmetry''')
    ys_idx = np.abs(y-ys).argmin()
    for jj, zj in enumerate(z):
        e, h = GetField_xy(S, res, zj)
        E_slice[jj,:,:] = np.array(e)[ys_idx,:,:]
        H_slice[jj,:,:] = np.array(h)[ys_idx,:,:]
    return E_slice, H_slice

def GetField_yz(S, x, y, z, xs=0.0):
    """
    Get E and H fields in an yz slice.
    Wraps around GetField_xy which is hard-coded in the C code and thus much faster.
    *BUT* you should check that it matches with GetSlice at least once because 
    the C function GetFieldsOnGrid tends to be touchy.
    """
    E_slice = np.empty((len(z), len(y), 3), dtype=np.complex128)
    H_slice = np.empty((len(z), len(y), 3), dtype=np.complex128)
    ym_yz, zm_yz = np.meshgrid(y,z)
    res = (len(x), len(y))
    xs_idx = np.abs(x-xs).argmin()
    for jj, zj in enumerate(z):
        e, h = GetField_xy(S, res, zj)
        E_slice[jj,:,:] = np.array(e)[:,xs_idx,:]
        H_slice[jj,:,:] = np.array(h)[:,xs_idx,:]
    return E_slice, H_slice


def GetSlice(S, ax1, ax2, axs=0, plane='xz', mode='Both'):
    """
    Slice field and/or epsilon in an arbitrary slice along 2 axis at a given 
    position on the third.
    Point-by-point computation version.
    There is no built-in function for this, so it is much slower. However it is 
    less ambiguous than the GetField_ methods and can be used to check/debug.
    
    Parameters
    ----------
    S: S4.Simulation
    ax1: 1Darray
        first axis to define the slice plane
    ax2: 1Darray
        second axis to define the slice plane
    axs: float
        coordinate along thirs axis on which to slice
    plane: str
        'xz', 'yz' or 'xz' to choose which slice it is
    mode: str
        either 'Field', 'Epsilon', 'All' to choose which quantity to slice
    """
    if mode=='Field' or mode=='All':
        E_slice = np.empty((len(ax2), len(ax1), 3), dtype=np.complex128)
        H_slice = np.empty((len(ax2), len(ax1), 3), dtype=np.complex128)
    elif mode=='Epsilon' or mode=='All':
        Eps_slice = np.empty((len(ax2), len(ax1)), dtype=np.complex128)
    a1m, a2m = np.meshgrid(ax1,ax2)
    for ii, a1i in enumerate(ax1):
        for jj, a2j in enumerate(ax2):
            if plane=='xz':
                setcoord = (a1i, axs, a2j)
            elif plane=='yz':
                setcoord = (axs, a1i, a2j)
            elif plane=='xy':
                setcoord = (a1i, a2j, axs)

            if mode=='Field' or mode=='All':
                E_slice[jj,ii,:],H_slice[jj,ii,:] = S.GetFields(*setcoord)
            elif mode=='Epsilon' or mode=='All':
                Eps_slice[jj,ii] = S.GetEpsilon(*setcoord)
    if mode=='Field':
        return E_slice, H_slice
    elif mode=='Epsilon':
        return Eps_slice
    elif mode=='All':
        return E_slice, H_slice, Eps_slice


class SlicePlot():
    """
    Helper to plot a slice with linecuts at arbitrary positions
    
    Parameters
    ----------
    a1: 1Darray
        1Darray for the *horizontal* axis of the 2D plot (typically x)
    a2: 1Darray
        1Darray for the *vertical* axis of the 2D plot (typically y or z)
    Field: 2Darray
        2Darray to be plotted (len(a2) x len(a1))
    hcoord: float
        horizontal coordinate for the vertical slice
    vcoord: float
        vertical coordinate for the horizontal slice
    sym: bool
        whether to force symmetry in the colormap or leav it automated
    Notes
    -----
    For a 3D stack of slices one should move to Plotly instead of 
    matplotlib
    https://nbviewer.jupyter.org/github/empet/Plotly-plots/blob/master/Plotly-Slice-in-volumetric-data.ipynb
    """
    def __init__(self,a1, a2, Field, cmap = plt.cm.bwr, 
                 hcoord=None, vcoord=None, sym=True):
        a1m, a2m = np.meshgrid(a1,a2) ## len(a2) lines and len(a1) columns !
        
        figratio = 5 ## general aspect ratio of the graph
        
        self.figslice = plt.figure() ## figure instance
        gs = gridspec.GridSpec(2,2, width_ratios=[figratio,1], 
                               height_ratios=[1,figratio])
        gs2D = gridspec.GridSpecFromSubplotSpec(1,2, subplot_spec = gs[1,0],
                                                width_ratios=[15,1])
        
        ## the 2D plot with the plane slice
        self.ax_2D = self.figslice.add_subplot(gs2D[0,0]) 
        if sym:
            absmax = np.abs(Field).max() # for normalization purpose
            vmin = -absmax
            vmax = absmax
        else:
            vmin = Field.min()
            vmax = Field.max()
        self.cax = self.ax_2D.pcolormesh(a1m, a2m, Field,
                               vmin=vmin, vmax=vmax, # symmetrize
                               cmap = cmap)
        axcbar = self.figslice.add_subplot(gs2D[0,1])
        
        ## horizontal slice at a given position along a2
        gs_a2slice = gridspec.GridSpecFromSubplotSpec(1,2, subplot_spec = gs[0,0],
                                                      width_ratios=[15,1])
        self.ax_a2slice = self.figslice.add_subplot(gs_a2slice[0,0], 
                                                    sharex=self.ax_2D)
        if vcoord is None: # if None, target center
            a2cen = (a2.min()+a2.max())/2.
            a2slice = np.abs(a2-a2cen).argmin()
        else:
            a2slice = np.abs(a2-vcoord).argmin()
        Fieldslice = Field[a2slice,:]
        self.ax_a2slice.plot(a1, Fieldslice)
        self.ax_a2slice.set_ylabel('Field')
        if sym:
            absmax = np.abs(Fieldslice).max()
            self.ax_a2slice.set_ylim([-absmax,absmax])
        
        ## vertical slice at a given position along a1
        self.ax_a1slice = self.figslice.add_subplot(gs[1,1], 
                                                    sharey=self.ax_2D)
        if hcoord is None: # if None target center
            a1cen = (a1.min()+a1.max())/2.
            a1slice = np.abs(a1-a1cen).argmin()
        else:
            a1slice = np.abs(a1-hcoord).argmin()
        Fieldslice = Field[:,a1slice]
        self.ax_a1slice.plot(Fieldslice, a2)
        self.ax_a1slice.set_xlabel('Field')
        if sym:
            absmax = np.abs(Field[:,a1slice]).max()
            self.ax_a1slice.set_xlim([-absmax,absmax])
            
        ## show the slices on 2D graph
        self.ax_2D.plot([a1[a1slice], a1[a1slice]], [a2.min(), a2.max()], 'k--')
        self.ax_2D.plot([a1.min(), a1.max()], [a2[a2slice], a2[a2slice]], 'k--')
        self.ax_2D.invert_yaxis()
        self.ax_2D.axis('tight')
        
        self.figslice.colorbar(self.cax, axcbar)
        self.figslice.tight_layout()


def Plot3DSlices(axes, points, Slices):
    """
    3D rendering of 3 slices in orthogonal planes.
    
    Parameters
    ----------
    axes: 3-tuple
        3-tuple of 1Darrays x,y,z for coordinates
    points: 3-tuple
        3-tuple of floats for the xs,ys,zs coordinates of the slices
    Slices: 3-tuple
        3-tuple of 2Darrays of the field sliced in xy, yz and xz planes
    """
    ### unpack
    x, y, z = axes
    xs, ys, zs = points
    Fxy, Fyz, Fxz = Slices
    
    xx, yy = np.meshgrid(x,y) ## in-plane coordinates (not reused)
    zzs = np.ones(xx.shape)*zs # altitude of the slice
    grid1 = pv.StructuredGrid(xx,yy,zzs) # make a mesh
    grid1['values'] = (Fxy.T).ravel() # set the color. why the ravel() ?
    xx, zz = np.meshgrid(x,z) ## vertical slice 1
    yys = np.ones(xx.shape)*ys ## y-coord of the slice
    grid2 = pv.StructuredGrid(xx,yys,zz)
    grid2['values'] = (Fxz).T.ravel()
    yy, zz = np.meshgrid(y,z)
    xxs = np.ones(yy.shape)*xs
    grid3 = pv.StructuredGrid(xxs,yy,zz)
    grid3['values'] = (Fyz).T.ravel()
    
    absmax = max(np.abs(Fxy).max(), 
                 np.abs(Fyz).max(),
                 np.abs(Fxz).max())

    plotter = pv.BackgroundPlotter()
    plotter.add_mesh(grid1, cmap=plt.cm.bwr,
                     clim=(-absmax,absmax))
    plotter.add_mesh(grid2, cmap=plt.cm.bwr,
                     clim=(-absmax,absmax))
    plotter.add_mesh(grid3, cmap=plt.cm.bwr,
                     clim=(-absmax,absmax))
    plotter.set_viewup((1,1,-1))
    isocam = plotter.camera.GetPosition()
    plotter.camera.SetPosition((isocam[0], isocam[1], 0))
    plotter.show_axes()
    plotter.show_grid()
