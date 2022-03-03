#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for converting results from Curve-Number simulations to Green-Ampt parameters
and vice-versa.



Taken From:
Baiamonte, G. (2019). SCS curve number and green-Ampt infiltration models. 
Journal of Hydrologic Engineering, 24(10), 04019034.
https://ascelibrary.org/doi/pdf/10.1061/%28ASCE%29HE.1943-5584.0001838?casa_token=W7b6GRxE6j0AAAAA:t5kZ935hvintrt5mUAcvaFVIwcAuTG085rj6fFi2peXxQJUpWJjPIwS1BAkfN1uMoB5sxjmfRw

"""

import numpy as np
import math
from math import sqrt

def GA_t_scale(delt_vmc, matric_pot, Ksat ):
    '''
    Calculate the sorbitivity timescale for the soil. 

    Parameters
    ----------
  delt_vmc :
      Change in volumetric moisture content across the wetting front (mm/mm)

      
  matric_pot :  matric potential of a soil in pascals
  
  Ksat : saturated hydraulic conductivity (mm/hr)
 

    Returns
    -------
    tc_ga : Characteristic timescale of the infiltration process.
   
    
    Equation 12B

    '''
    
    return delt_vmc * matric_pot /Ksat

#def CN_est1(total_rfall, ro_coef_GA, lamb = .2):
#    num = 2 * total_rfall * (1-ro_coef_GA)
#    denom = ro_coef_GA *(1-lamb) + 2 * lamb + sqrt(ro_coef_GA) + sqrt(ro_coef_GA * (1-lamb)**2 + 4 * lamb)
#    return 25400 / (254 + (num/denom))

def CN_est(t, Ksat, rho, norm_duration, ro_coef_GA, lamb = .2):
    '''
    

    Parameters
    ----------
    t : time (hrs)
    Ksat : saturated hydraulic conductivity (mm/hr)
        
    rho: rainfall intensity noramlized to the 
    Saturated Hydraulic Conductivity of the soil (unitless)
        
    norm_duration :  normalized duration of rainfall (in hours)
    calculated with function calc_norm_duration
        
    runoff_coef : fraction of rainfall lost as runoff. (between 0 and 1)
    
    lamb : ratio of initial abstractions to the S parameter of the SCS CN equation.
     The default is .2.

    Returns
    -------
    Estimated Curve Number
        
        Equation #22b

    '''
    
    num = t* Ksat * (1-ro_coef_GA) * rho * (rho * norm_duration +1 - lamb)
    denom = rho * norm_duration * (1+lamb) - lamb**2 
    
        
    
    return 25400 / (254 +(num/denom))


def calc_rho(Ksat, r_intens):
    '''
    Return the rainfall intensity noramlized to the 
    Saturated Hydraulic Conductivity of the soil

    Parameters
    ----------
    Ksat : saturated hydraulic conductivity (mm/hr)
        

        
    r_intens : rainfall intensity (mm/hr)
        

    Returns
    -------
    rho: rainfall intensity noramlized to the 
    Saturated Hydraulic Conductivity of the soil (unitless)
    
    Immediately after equation #13.
    '''
    
    return r_intens/Ksat

def calc_norm_duration(ro_coef_GA, rho, lamb = .2):
    '''
    Calculate the normalized duration of rainfall in hours.      
    
    
    Parameters
    ----------
    ro_coef_GA : runoff as a fraction of rainfall in the storm, 
                as simulated using the green-ampt-model
                
      rho: rainfall intensity noramlized to the 
      Saturated Hydraulic Conductivity of the soil (unitless). 
      Can be calculated with the function calc_rho.
      
    lamb : ratio of initial abstractions to the S parameter of the SCS CN equation.
     The default is .2.

    Returns
    -------
    normalized duration of rainfall (in hours)
    
    Equation #16

    '''
    
    num = ro_coef_GA * (1- lamb) + 2*lamb + sqrt(ro_coef_GA) * sqrt(ro_coef_GA * (1-lamb)+4*lamb)
    denom = 2*rho * (1-ro_coef_GA)
    return num/denom

def KsMax(GA_t_scale, r_intense, t):
    '''
    Calculate maximum Ksat value compatible with any runoff from the soil.
    

    Parameters
    ----------
    GA_t_scale :  Characteristic timescale of the infiltration process for the soil.
    
    
     r_intense : intensity of rainfall (mm/hr)
        DESCRIPTION.
        
    t : time of rainfall (hrs)

    Returns
    -------
    Theoretical maximum value of Ksat for the soil.
        
        Equation 23b

    '''
    
    return (r_intense * (sqrt(t)* sqrt(t+ 4* GA_t_scale)-t )) / 2*GA_t_scale
    
#%%
def CN_min(Ksat_max, hrs, rho, lamb =.2):
    '''
    Calculate the minimum estimated Curve Number for the soil based on the 
    rainfall simulation. 

    Parameters
    ----------
    Ksat_max : Max value for saturated Hydraulic Conductivity (mm/mm)
    calculated using Ks_max
    
   
    
   rho: rainfall intensity noramlized to the 
   Saturated Hydraulic Conductivity of the soil (unitless). 
   Can be calculated with the function calc_rho.
   
 lamb : ratio of initial abstractions to the S parameter of the SCS CN equation.
  The default is .2.

    Returns
    -------
    Estimate of the minimum value of the CN for the soil.

    Equation #24
    '''
    
    return 25400 / \
        (254 + (hrs *Ksat_max * rho)/ (lamb * (1+lamb)- lamb**2))
        
        
def estimate_CN(delt_vmc, Ksat, matric_pot, hrs, r_intens, ro_coef_GA, lamb = .2):
    '''
    Run all steps for estimating the effective curve-number of a soil 
    based on Green-Ampt Simlations and parameters.

    Parameters
    ----------
    delt_vmc :
        Change in volumetric moisture content across the wetting front (mm/mm)
    Ksat : saturated hydraulic conductivity (mm/hr)
        
    matric_pot :  matric potential of a soil in pascals
    hrs : time of storm event (in hrs)
        
    r_intens : rainfall intensity (mm/hr)
        
    ro_coef_GA : runoff as a fraction of rainfall in the storm, 
                as simulated using the green-ampt-model
    lamb : TYPE, optional
        DESCRIPTION. The default is .2.

    Returns
    -------
    An estimated equivalent Curve Number

    '''
    
    rho = calc_rho(Ksat, r_intens)
    norm_duration = calc_norm_duration(ro_coef_GA, rho, lamb)
    return CN_est(hrs, Ksat, rho, norm_duration, ro_coef_GA, lamb )


    
    