# -*- coding: utf-8 -*-
"""
A few functions for calculating erosion using the Universal Soil Loss Equation
(USLE)
"""

from numpy import product, exp

def usle(R, K, C, P, LS):
    '''Simplest form of the Universal Soil Loss Equation
    

    Parameters
    ----------
    R : Rainfall Erosivity Parameter
        
    K : Soil erodibility factor: between 0 and 1
        
    C : Cover management factor: between 0 and 1
        
    P : Practice factor: between 0 and 1
        DESCRIPTION.
    LS : Length-slope factor.
        Can be calculated with function calc_ls

    Returns
    -------
    Soil Erosion (metric tons/acre/yr)

    '''
    
    return product([R, K, C, P, LS])





def calc_CFRG(perc_coarse_frag):
    '''
    Calculate the CFRG factor for the MUSLE, based on the % coarse rock fragments.
    
    From the SWAT Technical docs, page 259

    Parameters
    ----------
    (perc_coarse_frag : numeric
        percent rock fragments in the soil.

    Returns
    -------
    float
        DESCRIPTION.

    '''
    cfrg = exp(-0.053 * perc_coarse_frag)
    return cfrg


def musle(R_factor, K, C, P , LS, CFRG):
   '''Modified Universal soil loss equation. 
    

    Parameters
    ----------
    R_factor : float
        rainfall factor.
    
    K : float
       USLE K factor for soil erodibility. (between 0 and 1)
    C : float
        USLE C factor for crop-type. (between 0 and 1)
    P : float
        USLE P factor for support practice (between 0 and 1)
    LS : float
        USLE Length-slope factor (>0)
    CFRG : float
        Course fragment factor 

    Returns
    -------
    sediment yield. 
        
        
 Page 252 SWAT technical docs
   '''
   sed_yield = product([11.8, R_factor, K, C, P, LS, CFRG ])
   return sed_yield
   






def calc_LS(slope, slope_length, constant = 22.1):
    '''
    Calculate the LS factor for a field based on average slope 
    and the length of slope.

    Parameters
    ----------
    slope : numeric
        average slope in % grade.
    slope_length : slope length in meters
        DESCRIPTION.
    constant : TYPE, optional
        DESCRIPTION. The default is 22.1.

    Returns
    -------
    LS : USLE LS factor (unitless)

    '''
    
    if slope<1:
        exponent = .2
    elif slope<3:
        exponent = .3
    elif slope<5:
        exponent = .4
    else:
        exponent = .5
    LS = (0.065 + 0.0456*slope + 0.006541*slope**2) * (slope_length / constant)**exponent
    return LS