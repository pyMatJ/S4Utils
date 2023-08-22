#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 17:50:31 2020

@author: Mathieu, Paul

Various material definitions are available
"""

import numpy as np
from math import pi
from pathlib import Path
##################
### physical constants

eV = 1.6e-19
m0 = 9.1e-31
h_eV = 4.135e-15
hbar_eV = h_eV/(2*pi)
eps0 = 8.854e-12
c_const = 3e8

#%% Material-specific permittivity functions

def epsAu(f):
    """
    Au permittivity assuming a Drude metal
    
    Parameters
    ----------
    f: 1D array
        frequency (Hz)
    
    Returns
    -------
    eps: 1D array
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
    f: 1D array
        frequency (Hz)
    Returns
    -------
    eps: 1D array
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
    https://doi.org/10.1063/1.336070
    
    Parameters
    ----------
    f: 1D array
        frequency (Hz)
    
    Returns
    -------
    eps: 1D array
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

def epsInP(f):
    """
    Parameters from Lockwood - SSCommun 136 404-409 (2005)
    
    
    Parameters
    ----------
    f : TYPE
        DESCRIPTION.

    Returns
    -------
    eps: 1D array
        complex permittivity (len(f))

    """
    wT_InP = 303.62 * c_const * 1e2 * 2*np.pi
    wL_InP = 345.32 * c_const * 1e2 * 2*np.pi
    gammaL_InP = 0.95 * c_const * 1e2 * 2*np.pi
    gammaT_InP = 2.8 * c_const * 1e2 * 2*np.pi
    eps_inf_InP = 9.61
    
    eps = eps2Phonons(f, wT_InP, wL_InP, gammaT_InP, gammaT_InP, eps_inf_InP)
    return eps
    
def epsInGaAs(f):
    """
    InGaAs permittivity 
    NOTE: check parameter values. So far it matches my COMSOL parameters
    
    Parameters
    ----------
    f: 1D array
        frequency (Hz)
    
    Returns
    -------
    eps: 1D array
        complex permittivity (len(f))
    """
    wL1_InGaAs = 5.1549e+13
    wT1_InGaAs = 4.7690e+13
    g1_InGaAs = 4.1020e+11
    wL2_InGaAs = 4.4271e+13
    wT2_InGaAs = 4.2813e+13
    g2_InGaAs = 3.4943e+11
    epsInf_InGaAs = 11.6108

    eps1 = eps2Phonons(f, wT1_InGaAs, wL1_InGaAs, g1_InGaAs, g1_InGaAs, epsInf_InGaAs)
    eps2 = eps2Phonons(f, wT2_InGaAs, wL2_InGaAs, g2_InGaAs, g2_InGaAs, epsInf_InGaAs)
    eps = eps1*eps2/epsInf_InGaAs
    return eps

def epsAlInAs(f):
    """
    InGaAs permittivity 
    NOTE: check parameter values. So far it matches my COMSOL parameters
    
    Parameters
    ----------
    f: 1D array
        frequency (Hz)
    
    Returns
    -------
    eps: 1D array
        complex permittivity (len(f))
    """

    wL1_AlInAs = 6.7972e+13
    wT1_AlInAs = 6.5693e+13
    g1_AlInAs = 5.4694e+11
    wL2_AlInAs = 4.5061e+13
    wT2_AlInAs = 4.2539e+13
    g2_AlInAs = 3.6462e+11
    epsInf_AlInAs = 10.2868
    
    eps1 = eps2Phonons(f, wT1_AlInAs, wL1_AlInAs, g1_AlInAs, g1_AlInAs, epsInf_AlInAs)
    eps2 = eps2Phonons(f, wT2_AlInAs, wL2_AlInAs, g2_AlInAs, g2_AlInAs, epsInf_AlInAs)
    eps = eps1*eps2/epsInf_AlInAs
    return eps    

def epsGaNx(f): 
    """
    Gallium nitride (GaN) ordinary axis permittivity (without excitons)
    
    Parameters
    ----------
    f: 1D array
        frequency (Hz)
    Returns
    -------
    eps: 1D array
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

    Parameters
    ------------
    f: 1D array
        frequency (Hz)
    Returns
    -------
    eps: 1D array
        complex permittivity (len(f))
    
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

    Parameters
    ------------
    f: 1D array
        frequency (Hz)
    Returns
    -------
    eps: 1D array
        complex permittivity (len(f))
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

    Parameters
    ------------
    f: 1D array
        frequency (Hz)
    Returns
    -------
    eps: 1D array
        complex permittivity (len(f))
    """
    wlAlNz     = 888.9*c_const*1e2*2*np.pi
    wtAlNz     = 608.5*c_const*1e2*2*np.pi
    GtAlNz     = 2.2*c_const*1e2*2*np.pi
    epsinfAlNz = 4.350
    eps = epsLorentzPhonon(f, wtAlNz, wlAlNz, GtAlNz, epsinfAlNz)
    return eps

def epsSiN(f):
    """
    Dielectric properties of low stress SiN (PECVD)
    from Cataldo et al., Optics Letters, Vol. 37 No. 20, pp. 4200-4202 (2012)

    Parameters
    ----------
    f : 1D array
        Frequency (Hz).

    Returns
    -------
    eps: 1D array
        complex permittivity (len(f))
    """
    eps_prime_j = np.array([7.582, 6.754, 6.601, 5.430, 4.601, 4.562])
    eps_sec_j = np.array([0, 0.3759, 0.0041, 0.1179, 0.2073, 0.0124])
    omega_T_j = 2*np.pi*np.array([13.913, 15.053, 24.521, 26.440, 31.724])
    Gamma_T_j = 2*np.pi*np.array([5.810, 6.436, 2.751, 3.482, 5.948])
    alpha_j = np.array([0.0001, 0.3427, 0.0006, 0.0002, 0.008])
    Deps = -np.diff(eps_prime_j) - 1.0j*np.diff(eps_sec_j)
    
    omega = 2*np.pi*f*1e-12 ### THz is the base unit in the paper
    
    Gamma_j = np.zeros((len(Gamma_T_j), len(omega)), dtype=np.complex128)
    eps = (eps_prime_j[-1] + 1.0j*eps_sec_j[-1]) * np.ones(len(omega), dtype=np.complex128)
    for jj, gamma_j in enumerate(Gamma_T_j):
        Gamma_j[jj] = gamma_j*np.exp(-alpha_j[jj]*((omega_T_j[jj]**2-omega**2)/(omega*gamma_j))**2)
        eps += Deps[jj]*omega_T_j[jj]**2/(omega_T_j[jj]**2 - omega**2 - 1.0j*omega*Gamma_j[jj])
    
    return eps
        
def epsCaF2(f):
    """
    CaF2 permittivity in the 150nm - 12 µm region
    From https://refractiveindex.info/?shelf=main&book=CaF2&page=Li

    Parameters
    ----------
    f : 1D array
        Frequency (in Hz).

    Returns
    -------
    eps : 1D array
        Complex permittivity.

    """
    lbda_um = c_const/f*1e6 ## Sellmeier equation below defined for lbda in µm
    n=(1+0.33973+0.69913/(1-(0.09374/lbda_um)**2)+0.11994/(1-(21.18/lbda_um)**2)+4.35181/(1-(38.46/lbda_um)**2))**.5
    eps = n**2
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
    
    epszz_inv = ezz_inv-1./eps_w*omega_p**2*fw*f12/(omega_tilde**2-omega**2-2*1.0j*omega*gamma_isb)
    eps_zz = 1./epszz_inv
    
    eps_xx = exx - omega_p**2*fw*eps_w/(omega**2-1.0j*omega*gamma_isb)
    return eps_xx, eps_zz
    
def epsPQW(f):
    """
    Doped Parabolic Quantum Well permittivity
    Paul?
    
    Parameters
    ----------
    f: 1D array
    
    Returns
    -------
    None
    
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
    eps = eps_inf*(wL**2-(2*np.pi*f)**2-1.0j*gammaL*(2*np.pi*f))/(wT**2-(2*np.pi*f)**2-1.0j*gammaT*(2*np.pi*f))
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

