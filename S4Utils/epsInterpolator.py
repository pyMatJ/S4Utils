#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 09:12:22 2020

@author: mathieu
"""
import numpy as np

matdata = '/mnt/Data/ODIN/Calculs/MaterialData/'


def eps_interpolator(f, fp, x1, x2, nk=False):
    """
    Returns the complex permittivity interpolated from tabulated data, either 
    in (n,k) files or in (eps1,eps2) files.

    Parameters
    ----------
    f : 1D array
        1D array of frequency at which to interpolate the data
    fp : 1D array
        Array of frequency for the tabulated data
    x1 : 1D array
        Array of tabulated data (n or eps1)
    x2 : 1D array
        Array of tabulated data (k or eps2)
    nk : Bool, optional
        Flag. True if the data is in (n,k) form, False if in (eps1, eps2) form. 
        Defaults to False

    Returns
    -------
    eps : 1D array
        1D array of (complex) permittivity (len(f))

    """
    if nk==False:
        return np.interp(f, fp, x1)+1.0j*np.interp(f, fp, x2)
    else:
        e1 = x1**2-x2**2
        e2 = 2*x1*x2
        return np.interp(f, fp, e1)+1.0j*np.interp(f, fp, e2)

def eps_Al2O3_o(f):
    fAl2O3_o, nAl2O3_o = np.genfromtxt(matdata+'/Bulk/Al2O3_no_Querry.txt', 
                                       skip_header=1, unpack=True)
    fAl2O3_o, kAl2O3_o = np.genfromtxt(matdata+'/Bulk/Al2O3_ko_Querry.txt', 
                                       skip_header=1, unpack=True)
    return eps_interpolator(f, fAl2O3_o, nAl2O3_o, kAl2O3_o, nk=True)

def eps_Al2O3_e(f):
    fAl2O3_e, nAl2O3_e = np.genfromtxt(matdata+'/Bulk/Al2O3_ne_Querry.txt', 
                                       skip_header=1, unpack=True)
    fAl2O3_e, kAl2O3_e = np.genfromtxt(matdata+'/Bulk/Al2O3_ke_Querry.txt', 
                                       skip_header=1, unpack=True)
    return eps_interpolator(f, fAl2O3_e, nAl2O3_e, kAl2O3_e, nk=True)
    
def eps_AlN(f):
    fAlN, nAlN = np.genfromtxt(matdata+'/Bulk/AlN_n_Kischkat.txt', 
                               skip_header=1, unpack=True)
    fAlN, kAlN = np.genfromtxt(matdata+'/Bulk/AlN_k_Kischkat.txt', 
                               skip_header=1, unpack=True)
    return eps_interpolator(f, fAlN, nAlN, kAlN, nk=True)

def eps_GaN(f):
    fGaN, nGaN = np.genfromtxt(matdata+'/Bulk/GaN_n_interp.txt', 
                               skip_header=1, unpack=True)
    fGaN, kGaN = np.genfromtxt(matdata+'/Bulk/GaN_k_interp.txt', 
                               skip_header=1, unpack=True)
    return eps_interpolator(f, fGaN, nGaN, kGaN, nk=True)

def eps_GaAs(f):
    fGaAs, nGaAs = np.genfromtxt(matdata+'/Bulk/GaAs_n.txt', 
                                 skip_header=1, unpack=True)
    fGaAs, kGaAs = np.genfromtxt(matdata+'/Bulk/GaAs_k.txt', 
                                 skip_header=1, unpack=True)
    return eps_interpolator(f, fGaAs, nGaAs, kGaAs, nk=True)

def eps_Au(f):
    """
    Ordal tabulated values for gold

    Parameters
    ----------
    f : 1D array
        Frequency

    Returns
    -------
    eps : 1D array
        Complex permittivity of Au (Ordal tabulated data)

    """
    fAu, nAu = np.genfromtxt(matdata+'/Bulk/Au_n_Ordal.txt', 
                             skip_header=1, unpack=True)
    fAu, kAu = np.genfromtxt(matdata+'/Bulk/Au_k_Ordal.txt', 
                             skip_header=1, unpack=True)
    return eps_interpolator(f, fAu, nAu, kAu, nk=True)

