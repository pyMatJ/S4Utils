#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 17:50:31 2020

@author: Mathieu, Paul

Various material definitiosn are available

**Warning** 08/05/2020
all definitions are not consistent in terms of imaginary part
You have to ensure that you have a *positive* imaginary part for a *lossy* 
material. np.conj solves your problems.
"""

import numpy as np

##################
### physical constants
eV = 1.6e-19
m0 = 9.1e-31
h_eV = 4.135e-15
hbar_eV = h_eV/(2*np.pi)
eps0 = 8.854e-12
c_const = 3e8

#%% Material-specific permittivity functions

def epsAu(f):
    """
    Au permittivity assuming a Drude metal
    Parameters
    ----------
    f : 1D array
        frequency (Hz)
    Returns
    -------
    eps : 1D array
        complex permittivity (len(f))
    """
    wp = 2*np.pi*2.17e15
    gamma = 2*np.pi*6.48e12
    eps = 1.0-(wp**2)/((2*np.pi*f)**2+1.0j*gamma*(2*np.pi*f))
    return eps

def epsGaAs(f):
    """
    GaAs permittivity assuming a Lorentz model with 2 phonons
    Parameters
    ----------
    f : 1D array
        frequency (Hz)
    Returns
    -------
    eps : 1D array
        complex permittivity (len(f))
    """    
    wL =  292.1*2.99792458e10*2*np.pi
    wT = 267.8*2.99792458e10*2*np.pi
    gamma = 1.86*2.99792458e10*2*np.pi
    eps_inf = 10.89
    eps = epsLorentzPhonon(f, wT, wL, gamma, eps_inf)
    return eps

def epsAlGaAs(f, xAl):
    """
    AlGaAs permittivity assuming a double 2 phonons model
    following Adachi, J. Appl. Phys. 58 R1-R29 (1985), table I
    Parameters
    ----------
    f : 1D array
        frequency (Hz)
    Returns
    -------
    eps : 1D array
        complex permittivity (len(f))
    """
    eps_i = 10.89-2.73*xAl # ioffe
    wL_GaAs = (36.25-6.55*xAl+1.79*xAl**2)*1e-3/hbar_eV
    wT_GaAs = (33.29-0.64*xAl-1.16*xAl**2)*1e-3/hbar_eV
    wL_AlAs = (44.63+8.78*xAl-3.32*xAl**2)*1e-3/hbar_eV
    wT_AlAs = (44.63+0.55*xAl-0.3*xAl**2)*1e-3/hbar_eV
    gamma = 1.86*2.99792458e10*2*np.pi*3 ## PRL 56 1842 (1986)
    eps = xAl*epsLorentzPhonon(f, wT_AlAs, wL_AlAs, gamma, eps_i)
    eps = eps+(1-xAl)*epsLorentzPhonon(f, wT_GaAs, wL_GaAs, gamma, eps_i)
    return eps
    
def epsGaNx(f): 
    """
    Gallium nitride (GaN) ordinary axis permittivity (without excitons)
    Parameters
    ----------
    f : 1D array
        frequency (Hz)
    Returns
    -------
    eps : 1D array
        complex permittivity (len(f))
    """
    wlGaN      = 742.1*c_const*1e2*2*np.pi
    wtGaN     = 560.1*c_const*1e2*2*np.pi
    GtGaN     = 4.0*c_const*1e2*2*np.pi
    GlGaN     = 4.0*c_const*1e2*2*np.pi
    epsinfGaN = 5.04
    eps =  eps2Phonons(f,wtGaN,wlGaN,GtGaN,GlGaN,epsinfGaN)
    return eps
def epsGaNz(f): # GaN extraordinary epsilon
    """
    Gallium nitride (GaN) extraordinary axis permittivity (without excitons)

    :param array f: frequency (array or float)
    :return: permittivity (float or len(f)-array)
    """
    wlGaNz     = 732.5*c_const*1e2*2*np.pi
    wtGaNz     = 537*c_const*1e2*2*np.pi
    GtGaNz     = 4.0*c_const*1e2*2*np.pi
    GlGaNz     = 4.0*c_const*1e2*2*np.pi
    epsinfGaNz = 5.01
    eps =  eps2Phonons(f,wtGaNz,wlGaNz,GtGaNz,GlGaNz,epsinfGaNz)
    return eps

def epsAlNx(f):
    """
    Aluminium nitride (AlN) ordinary axis permittivity (without excitons)

    :param array f: frequency (array or float)
    :return: permittivity (float or len(f)-array)
    """
    wlAlN     = 909.6*c_const*1e2*2*np.pi
    wtAlN     = 667.2*c_const*1e2*2*np.pi
    GtAlN     = 2.2*c_const*1e2*2*np.pi
    epsinfAlN = 4.160
    eps = epsLorentzPhonon(f, wtAlN, wlAlN, GtAlN, epsinfAlN)
    return eps

def epsAlNz(f):
    """
    Aliuminium nitride (AlN) extraordinary axis permittivity (without excitons)

    :param array f: frequency (array or float)
    :return: permittivity (float or len(f)-array)
    """
    wlAlNz     = 888.9*c_const*1e2*2*np.pi
    wtAlNz     = 608.5*c_const*1e2*2*np.pi
    GtAlNz     = 2.2*c_const*1e2*2*np.pi
    epsinfAlNz = 4.350
    eps = epsLorentzPhonon(f, wtAlNz, wlAlNz, GtAlNz, epsinfAlNz)
    return eps

#%% Model permittivity functions

def epsZal(f, eps_w, eps_b, omega_isb, gamma_isb, omega_p, fw, f12=0.96):
    """
    Zaluzny model for an active region composed of square QWs (eps_w) and 
    barriers (eps_b)
    Parameters
    ----------
    f : 1D array or float
        frequency (Hz)
    eps_w : len(f) array
        complex permittivity array of the well permittivity
    eps_barr : len(f) array
        complex permittivity array of the well permittivity
    omega_isb : float
        ISB transition pulsation (rad/s) *without plasma shift*
    gamma_isb : float
        ISB broadening
    omega_p : float
        plasma pulsation (rad/s)
    fw : float
        filling fraction (QW thick/period thick)
    f12 : float
        Oscillator strength. The default is 0.96 (infinite square well)

    Returns
    -------
    eps_xx : len(f) array
        complex in-plane permittivity
    eps_zz : len(f) array
        complex out-of-plane permittivity

    """
    omega = 2*np.pi*f
    ezz_inv = fw/eps_w+(1-fw)/eps_b # fractional inverse permittivity, term 1 in zaluzny
    exx = fw*eps_w+(1-fw)*eps_b
    omega_tilde = np.sqrt(omega_isb**2+omega_p**2)
    
    epszz_inv = ezz_inv-1./eps_w*omega_p**2*fw*f12/(omega_tilde**2-omega**2+2*1.0j*omega*gamma_isb)
    eps_zz = 1./epszz_inv
    
    eps_xx = exx - omega_p**2*fw*eps_w/(omega**2+1.0j*omega*gamma_isb)
    return eps_xx, eps_zz
    
def epsPQW(f):
    """
    Doped Parabolic Quantum Well permittivity
    Paul?
    """
    
def epsDrude(f, fp, gammap, epsinf=1.0):
    """
    Drude model for a material with free electrons
    Parameters
    ----------
    f : 1D array
        frequency (Hz)
    fp : float
        plasma frequency (Hz)
    gammap : float
        damping (collision rate) 
    epsinf : float
        high-frenquency permittivity. Defaults to 1.0 (metal).
    Returns
    -------
    eps : len(f) array
        complex permittivity 
    """    
    eps = epsinf-fp**2/(f**2-1.0j*gammap*f)
    return eps


def eps2Phonons(f, wT, wL, gammaT, gammaL, eps_inf):
    """
    2 phonons model for a semiconductor material 
    Parameters
    ----------
    f : 1D array
        frequency (Hz)
    wT : float
        Transverse phonon pulsation (rad/s)
    wL : float
        Longitudinal phonon pulsation (rad/s)
    gammaT : float
        Transverse phonon damping rate (rad/s)
    gammaL : float
        Longitudinal phonon damping rate (rad/s)
    eps_inf : float
        high-frenquency permittivity. Defaults to 1.0 (metal).
    Returns
    -------
    eps : len(f) array
        complex permittivity 
    """    
    eps = eps_inf*((2*np.pi*f)**2-wL**2+1.0j*gammaL*(2*np.pi*f))/((2*np.pi*f)**2-wT**2+1.0j*gammaT*(2*np.pi*f))
    return eps

def epsLorentzPhonon(f, wT, wL, gamma, eps_inf):
    """
    Lorentz phonon model for a semiconductor material 
    Parameters
    ----------
    f : 1D array
        frequency (Hz)
    wT : float
        Transverse phonon pulsation (rad/s)
    wL : float
        Longitudinal phonon pulsation (rad/s)
    gamma : float
        Phonon damping rate (rad/s)
    eps_inf : float
        high-frenquency permittivity. Defaults to 1.0 (metal).
    Returns
    -------
    eps : len(f) array
        complex permittivity 
    """    
    eps = (eps_inf*(1+(wL**2-wT**2)/(wT**2-(2*np.pi*f)**2-1.0j*gamma*(2*np.pi*f))))
    return eps
    

#%% Utility functions


def omegaP_2D(N2D, meff, eps_i, w_QW):
    """
    Plasma pulsation (in units of rad/s) for an equivalent 2D doping N2D
    Parameters
    ----------
    N2D : float
        equivalent sheet carrier density (m-2)
    meff : float
        effective mass coefficient (no units)
    eps_i : float
        high-frequency permittivity
    w_QW : float
        quantum well thickness (m)
    Returns
    -------
    wp : float
        plasma pulsation (rad/s)
    """
    wp = np.sqrt(N2D*eV**2/(meff*m0*eps0*eps_i*w_QW))
    return wp

def omegaP_3D(N3D, meff, eps_i):
    """
    Plasma pulsation (in units of rad/s) for an equivalent 3D doping N3D
    Parameters
    ----------
    N3D : float
        carrier concentration (m-3)
    meff : float
        effective mass coefficient (no units)
    eps_i : float
        high-frequency permittivity
    Returns
    -------
    wp : float
        plasma pulsation (rad/s)
    """
    wp = np.sqrt(N3D*eV**2/(meff*m0*eps0*eps_i))
    return wp

