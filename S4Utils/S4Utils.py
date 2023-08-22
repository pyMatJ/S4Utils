#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 12:05:13 2020

@author: mathieu
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pyvista as pv ## might not be installed
from matplotlib.widgets import Slider, CheckButtons

#%%  General utilities
def totuple(a):
    """
    Converts lists and array to a tuple form. Usefull since S4 mostly accepts 
    tuples instead of numpy arrays.
    
    Parameters
    ----------
    a: array-like, list or float
        Array-like parameter to cast to a tuple
    
    Returns
    -------
    ta: tuple
        elements of a cast into a tuple of same 'shape'.
        
    Raises
    ------
    TypeError:
        If the type of a cannot be cast into a tuple form, returns a
    """
    try:
        ta = tuple(totuple(i) for i in a)
        return ta
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
    
    Returns
    -------
        None.
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
        
#%% Extracting values from simulations        

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
    
    Returns
    -------
        None.
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
    
    Parameters
    ----------
    S: S4.Simulation Object
        Simulation to run
    x: 1Darray
        1Darray of x-coordinates
    y: 1Darray
        1Darray of y-coordinates
    z: 1Darray
        1Darray of z-coordinates
    ys: float, optional
        y-coordinate of the slice. Defaults to 0 (center of the unit cell).
    
    Returns
    -------
    E_slice: len(z) by len(x) by 3 array
        Array of the electric field vector on the slice
    H_slice: len(z) by len(x) by 3 array
        Array of the magnetic field vector on the slice

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

    Parameters
    ----------
    S: S4.Simulation Object
        Simulation to run
    x: 1Darray
        1Darray of x-coordinates
    y: 1Darray
        1Darray of y-coordinates
    z: 1Darray
        1Darray of z-coordinates
    xs: float, optional
        x-coordinate of the slice. Defaults to 0 (center of the unit cell).
    
    Returns
    -------
    E_slice: len(z) by len(y) by 3 array
        Array of the electric field vector on the slice
    H_slice: len(z) by len(y) by 3 array
        Array of the magnetic field vector on the slice
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



def GetSlice(S, ax1, ax2, axs=0, plane='xz', mode='All'):
    """
    Slice field and/or epsilon in an arbitrary slice along 2 axis at a given 
    position on the third.
    Point-by-point computation version.
    There is no built-in function for this, so it is much slower. However it is 
    less ambiguous than the GetField methods and can be used to check/debug.
    
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


def GetFieldVolume(S, x, y, z):
    """
    Get the electric field vector in the whole simulation volume

    Parameters
    ----------
    S : S4.Simulation obect
        Simulation object to run
    x : 1Darray
        1Darray of x-coordinates
    y : 1Darray
        1Darray of y-coordinates
    z : 1Darray
        1Darray of z-coordinates

    Returns
    -------
    E_vol : len(x) by len(y) by len(z) by 3 array
        Array of electric field in the simulation volume
    H_vol : len(x) by len(y) by len(z) by 3 array
        Array of magnetic field in the simulation volume

    """
    E_vol = np.empty((len(x), len(y), len(z), 3), dtype=np.complex128)
    H_vol = np.empty((len(x), len(y), len(z), 3), dtype=np.complex128)
    res = (len(x), len(y))
    for jj, zj in enumerate(z):
        progress = jj/len(z)*100
        print('\r z-altitude, Progress : {:.2f}'.format(progress), end='')
        e, h = GetField_xy(S, res, zj)
        E_vol[:,:,jj,:] = np.array(e)
        H_vol[:,:,jj,:] = np.array(h)
    return E_vol, H_vol


def GetEpsilonVolume(S, x, y, z):
    """
    Get the reconstructed permittivity in the whole simulation volume

    Parameters
    ----------
    S : S4.Simulation obect
        Simulation object to run
    x : 1Darray
        1Darray of x-coordinates
    y : 1Darray
        1Darray of y-coordinates
    z : 1Darray
        1Darray of z-coordinates

    Returns
    -------
    Eps_vol : len(x) by len(y) by len(z) array
        Array of the permittivity in the simulation volume
    
    """
    Eps_vol = np.empty((len(x), len(y), len(z)), dtype=np.complex128)
    for jj, zj in enumerate(z):
        progress = jj/len(z)*100
        print('\r z-altitude, Progress : {:.2f}'.format(progress), end='')
        Eps_slice = GetSlice(S, ax1=x, ax2=y, axs=zj, plane='xy', mode='Epsilon')
        Eps_vol[:,:,jj] = np.array(Eps_slice)
    return Eps_vol
    

#%% 2D Slices + GUI

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
    cmap: `matplotib.colors.Colormap` object
        The colormap to pass to `matplotlib.pcolormesh` for plotting. Defaults 
        to matplotlib.cm.bwr.
    hcoord: float, optional
        Horizontal coordinate for the vertical slice. Defaults to the center 
        of the axe coordinate array given (a1).
    vcoord: float, optional
        Vertical coordinate for the horizontal slice. Defaults to the center 
        of the axe coordinate array given (a2).
    sym: bool, optional
        Flag to force symmetry in the colormap or leave it automated
    
    Attributes
    ----------
    a1: 1Darray
        1Darray for the *horizontal* axis of the 2D plot (typically x)
    a2: 1Darray
        1Darray for the *vertical* axis of the 2D plot (typically y or z)
    Field: 2Darray
        2Darray to be plotted (len(a2) x len(a1))
    hcoord: float, optional
        Horizontal coordinate for the vertical slice. Defaults to the center 
        of the axe coordinate array given (a1).
    vcoord: float, optional
        Vertical coordinate for the horizontal slice. Defaults to the center 
        of the axe coordinate array given (a2).
    sym: bool, optional
        Flag to force symmetry in the colormap or leave it automated
    figslice: matplotlib.figure 
        The matplotlib figure in which the plot is drawn
    fignum: int
        The number of the figure as retrieved by matplotlib.figure.number 
    ax_2D: matplotlib.axes
        The axes in which the 2D slice is drawn
    cax: matplotlib.axes
        The axes in which the colorbar for ax_2D drawn
    ax_cbar: matplotlib.colorbar
        The matplotlib.colorbar instance
    ax_a1slice: matplotlib.axes
        The axes in which the slice along a2 at a fixed a1 position is drawn
    lv: matplotlib.line
        The line of the *vertical* cut, e.g. along a2
    ax_a2slice: matplotlib.axes
        The axes in which the slice along a1 at a fixed a2 position is drawn
    lh: matplotlib.line
        The line of the *horizontal* cut, e.g. along a1
    lslice_h: matplotlib.line
        The line showing the position of the horizontal slice in ax_2D
    lslice_v: matplotlib.line
        The line showing the position of the vertical slice in ax_2D
    GUIButton: matplotlib.widgets.CheckButtons
        The button to hide/show the GUI
        
    Notes
    -----
    For a 3D stack of slices one should move to Plotly instead of 
    matplotlib
    https://nbviewer.jupyter.org/github/empet/Plotly-plots/blob/master/Plotly-Slice-in-volumetric-data.ipynb
    """
    def __init__(self,a1, a2, Field, cmap = plt.cm.bwr, 
                 hcoord=None, vcoord=None, sym=True):
        
        ### initial values for the 1D slices
        self.hcoord = hcoord
        self.vcoord = vcoord
        ## axes
        self.a1 = a1 ## abscissa
        self.a2 = a2 ## ordinate
        self.Field = Field ## 2D "field" to plot
        self.sym = sym ## whether to force symmetrical colormap or not
        
        ## len(a2) lines and len(a1) columns !
        a1m, a2m = np.meshgrid(self.a1,self.a2) 
        figratio = 5 ## general aspect ratio of the graph
        
        ###############
        ## Main figure instance containing all plots
        self.figslice = plt.figure(constrained_layout=True) 
        ## keep the number to associate to the GUI figure
        self.fignum = self.figslice.number
        gs = gridspec.GridSpec(2,3, width_ratios=[figratio/20, figratio, 1], 
                               height_ratios=[1,figratio],
                               figure=self.figslice)
        
        ##############
        ## the 2D plot with the plane slice
        self.ax_2D = self.figslice.add_subplot(gs[1,1]) 
        ## set min and max with or without symmetry
        if self.sym:
            absmax = np.abs(self.Field).max() # for normalization purpose
            vmin = -absmax
            vmax = absmax
        else:
            vmin = self.Field.min()
            vmax = self.Field.max()
        self.cax = self.ax_2D.pcolormesh(a1m, a2m, self.Field,
                               vmin=vmin, vmax=vmax, # symmetric or not
                               cmap = cmap,
                               shading = 'auto')
        self.axcbar = self.figslice.add_subplot(gs[1,0]) ## colorbar
        
        ###############
        ## *horizontal slice* at a given position along a2
        self.ax_a2slice = self.figslice.add_subplot(gs[0,1], 
                                                    sharex=self.ax_2D)
        if self.vcoord is None: # if None, target center
            a2cen = (self.a2.min()+self.a2.max())/2.
            self.vcoord = a2cen
        a2slice = np.abs(self.a2-self.vcoord).argmin() ## index for the slice
        Fieldslice = self.Field[a2slice,:] ## slice data
        ## plot the data, keep reference to the line for GUI
        self.lh, = self.ax_a2slice.plot(self.a1, Fieldslice)  
        self.ax_a2slice.set_ylabel('Field')
        self._UpdateHSliceLims() ## properly sets the axes limits
        
        ##############
        ## vertical slice at a given position along a1
        self.ax_a1slice = self.figslice.add_subplot(gs[1,2], 
                                                    sharey=self.ax_2D)
        if self.hcoord is None: # if None target center
            a1cen = (self.a1.min()+self.a1.max())/2.
            self.hcoord = a1cen
        a1slice = np.abs(self.a1-self.hcoord).argmin() ## index for the slice
        Fieldslice = self.Field[:,a1slice] ## slice data
        ## plot the data, keep reference to the line for GUI
        self.lv, = self.ax_a1slice.plot(Fieldslice, self.a2)
        self.ax_a1slice.set_xlabel('Field')
        self._UpdateVSliceLims() ## properly sets the axes limits
        
        ###############
        ## show the slices on 2D graph
        self.lslice_h, = self.ax_2D.plot([self.a1[a1slice], self.a1[a1slice]], 
                                         [self.a2.min(), self.a2.max()], 'k--')
        self.lslice_v, = self.ax_2D.plot([self.a1.min(), self.a1.max()], 
                                         [self.a2[a2slice], self.a2[a2slice]], 'k--')
        ## flipping the vertical axis, mainly usefull for vertical slices
        self.ax_2D.invert_yaxis() 
        self.ax_2D.axis('tight')
        self.ax_2D.set_xlabel('ax1') ## can be changed externally
        self.ax_2D.set_ylabel('ax2') ## can be changed externally

        ## put the colorbar on the left of the figure
        self.figslice.colorbar(self.cax, cax=self.axcbar)
        self.axcbar.tick_params(axis='y', labelleft=True,
                                labelright=False,
                                right=False,
                                left=True)
        self.axcbar.set_ylabel('Field')
        self.axcbar.yaxis.set_label_position('left')
        
        #######
        ### Check button to activate/deactivate the GUI
        axGuiBut = self.figslice.add_subplot(gs[0,2])
        self.GUIButton = CheckButtons(axGuiBut, labels=['GUI'])
        self.GUIButton.on_clicked(self._ShowGUI)
        
        plt.show()
    
    def _UpdateHSliceLims(self):
        """
        Properly sets the axes limits of the horizontal slice in ax_a2slice 
        axes to the data plotted inside the axes

        Returns
        -------
        None.

        """
        slicedata = self.lh.get_data()[1]
        datalims = np.array([slicedata.min(), slicedata.max()])
        if self.sym:
            absmax = np.abs(datalims).max()
            datamin = -absmax
            datamax = absmax
        else:
            datamin = datalims.min()
            datamax = datalims.max()
        self.ax_a2slice.set_ylim([datamin,datamax])
        
    def _UpdateVSliceLims(self):
        """
        Properly sets the axes limits of the vertical slice in ax_a1slice 
        axes to the data plotted inside the axes

        Returns
        -------
        None.

        """
        slicedata = self.lv.get_data()[0]
        datalims = np.array([slicedata.min(), slicedata.max()])
        if self.sym:
            absmax = np.abs(datalims).max()
            datamin = -absmax
            datamax = absmax
        else:
            datamin = datalims.min()
            datamax = datalims.max()
        self.ax_a1slice.set_xlim([datamin,datamax])
        
    def _ShowGUI(self, label):
        """
        Show/Hide the GUI side figure which enables changing the slice locations.
        The figure number is used to differenciate between GUI of different 
        figures.
        
        Parameters
        ----------
        label : List
            List of the labels of the CheckButtons.

        Returns
        -------
        None.

        """
        ### check if the figure exists. If yes, hide it
        if plt.fignum_exists('GUIFig_%i'%self.fignum): 
            plt.close('GUIFig_%i'%self.fignum)
        ## if not, instanciate the GUISlicePlot
        else:
            self.GUI = GUISlicePlot(self)
    def savefig(self, fname, **kwargs):
        self.figslice.savefig(fname, **kwargs)
        
class GUISlicePlot():
    """
    Generates an auxiliary figure to interact with a SlicePlot.
    The GUI side plot allows to move the vertical and horizontal sliders 
    to show slices at different locations.

    Parameters
    ----------
    parentfig : :py:class:`S4Utils.S4Utils.SlicePlot` instance
        The SlicePlot figure to which the GUI applies.
    
    Attributes
    ----------
    parentfig : :py:class:`S4Utils.S4Utils.SlicePlot` instance
        The SlicePlot figure to which the GUI applies.
    figGUI: matplotlib.figure instance
        The matplotlib.figure in which the GUI is drawn
    ax1: matplotlib.axes
        The axes for the first slider linked to parentfig.a1
    ax2: matplotlib.axes
        The axes for the first slider linked to parentfig.a2
    Sl_ax1: matplotlib.widgets.Slider
        The slider linked to the parentfig.a1 axis
    Sl_ax2: matplotlib.widgets.Slider
        The slider linked to the parentfig.a2 axis
    
    
    Returns
    -------
        None. 
        
    .. todo:: 
        Add a frequency input textbox to change the frequency and update plot
    """
    def __init__(self, parentfig):
        ## parent figure to interact with (SlicePlot instance)
        self.parentfig = parentfig
        ## the parent figure number to associate a GUI to each figure
        parentnum = self.parentfig.fignum
        ## The actual GUI figure
        self.figGUI = plt.figure(figsize=(4,1), num='GUIFig_%i'%parentnum)
        gs = gridspec.GridSpec(2,1)
        self.ax1 = self.figGUI.add_subplot(gs[0]) ## axes for the slider1
        self.ax2 = self.figGUI.add_subplot(gs[1]) ## axes for the slider2
        
        ### slider 1 along abscissa
        ax1lims = self.parentfig.ax_2D.get_xlim() ## get axis limits
        ax1init = self.parentfig.hcoord ## get current slice coordinate
        ax1label = self.parentfig.ax_2D.get_xlabel() ## get current axe label
        ## create the slider
        self.Sl_ax1 = Slider(self.ax1, ax1label, valmin=min(ax1lims), valmax=max(ax1lims),
                             valinit=ax1init)
        ## callback function
        self.Sl_ax1.on_changed(self._update_vslice)
        ### slider 2 along ordinate
        ax2lims = self.parentfig.ax_2D.get_ylim() ## get axis limits
        ax2init = self.parentfig.vcoord ## get current slice coordinate
        ax2label = self.parentfig.ax_2D.get_ylabel() ## get current axe label
        ## create the slider
        self.Sl_ax2 = Slider(self.ax2, ax2label, valmin=min(ax2lims), valmax=max(ax2lims),
                       valinit=ax2init)
        ## callback function
        self.Sl_ax2.on_changed(self._update_hslice)
        
    def _update_hslice(self,val):
        """
        Callback function for the slider along the abscissa coordinate.

        Parameters
        ----------
        val : Float
            Value returned by the slider object

        Returns
        -------
        None.

        """
        coord = self.Sl_ax2.val ## cordinate from the slider value
        self.parentfig.vcoord=coord ## set this coordinate to the parent figure
        a2slice = np.abs(self.parentfig.a2-coord).argmin() ## get the index along the axe array
        Fieldslice = self.parentfig.Field[a2slice,:] ## slice data
        self.parentfig.lh.set_ydata(Fieldslice) ## set the new data to the referenced line
        self.parentfig.figslice.canvas.draw() ## update plot
        ### check if axes limits should be updated and update if necessary
        axlims = self.parentfig.ax_a2slice.get_ylim() 
        if (Fieldslice.min()<min(axlims)) or (Fieldslice.max()>max(axlims)):
            self.parentfig._UpdateHSliceLims() 
        ### plot the slice position on 2D axis        
        self.parentfig.lslice_v.set_ydata([self.parentfig.a2[a2slice], self.parentfig.a2[a2slice]])
        
    def _update_vslice(self,val):
        """
        Callback function for the slider along the abscissa coordinate.

        Parameters
        ----------
        val : Float
            Value returned by the slider object

        Returns
        -------
        None.

        """
        coord = self.Sl_ax1.val ## cordinate from the slider value
        self.parentfig.hcoord = coord ## set this coordinate to the parent figure
        a1slice = np.abs(self.parentfig.a1-coord).argmin() ## get the index along the axe array
        Fieldslice = self.parentfig.Field[:,a1slice] ## slice data
        self.parentfig.lv.set_xdata(Fieldslice) ## set the new data to the referenced line
        self.parentfig.figslice.canvas.draw() ## update plot
        ### check if axes limits should be updated and update if necessary
        axlims = self.parentfig.ax_a1slice.get_xlim()
        if (Fieldslice.min()<min(axlims)) or (Fieldslice.max()>max(axlims)):
            self.parentfig._UpdateVSliceLims()
        ## plot the slice position on 2D axis
        self.parentfig.lslice_h.set_xdata([self.parentfig.a1[a1slice], self.parentfig.a1[a1slice]])
        
#%% 3D Plotting
def Plot3DSlices(axes, DataVol, points, InteractiveSlices=False):
    """
    3D rendering of 3 slices in orthogonal planes.
    
    Parameters
    ----------
    axes:
        x, y and z 1Darray vectors of coordinates on which the data is 
        calculated. Usually passed as (x,y,z) or np.array((x,y,z))
        
    DataVol: len(x) by len(y) by len(z) numpy array
        Volumetric data in numpy array to be plotted
        Usually obtained using :py:func:`S4Utils.S4Utils.GetFieldVolume` or 
        :py:func:`S4Utils.S4Utils.GetEpsilonVolume` functions.
    
    points: tuple, list or array
        point at which to slice the volume data. Usually (xs, ys, zs), in 
        real space coordinates
    
    InteractiveSlices: bool (opt)
        Boolean flag to activate PyVista's defaults 3D slicing. Defaults to 
        False. 
        
    Returns
    -------
        None.
        
    Notes
    -----
        Follows https://docs.pyvista.org/examples/00-load/create-uniform-grid.html
        and https://docs.pyvista.org/plotting/widgets.html#pyvista.WidgetHelper.add_mesh_slice_orthogonal
    """
    xs, ys, zs = points ## coordinates for each slice
    x, y, z = axes ## axes
    xi = np.abs(x-xs).argmin() # index of x-coordinate of the yz slice
    yi = np.abs(y-ys).argmin() # index of y-coordinate of the xz slice
    zi = np.abs(z-zs).argmin() # index of z-coordinate of the xy slice

    grid = pv.UniformGrid() # initialize the 3D grid
    grid.dimensions = DataVol.shape 
    grid.point_arrays['values']=DataVol.flatten(order='F') 
    plotter = pv.BackgroundPlotter()
    if not InteractiveSlices:
        slices =  grid.slice_orthogonal(x=xi, y=yi, z=zi)
        plotter.add_mesh(slices, cmap='bwr')
    else:
        ### this alone does the trick but I find it hard to manipulate
        ## maybe leave as an option
        plotter.add_mesh_slice_orthogonal(grid, cmap='bwr') 
    ##########
    ## plotter options
    plotter.set_viewup((1,1,-1)) ## flip z-axis view to match S4 coordinates
    isocam = plotter.camera.GetPosition() ## get initial camera location
    plotter.camera.SetPosition((isocam[0], isocam[1], 0)) ## set new camera location
    plotter.show_axes() ## show axes orientation widget
    plotter.show_grid() ## show the grid and axes
    
