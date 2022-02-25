#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 08:02:33 2021

@author: bdube

Equations from:
https://ascelibrary-org.ezproxy.uvm.edu/doi/pdf/10.1061/%28ASCE%29HE.1943-5584.0001838
"""

import numpy as np
import math
from math import sqrt

def GA_t_scale(delt_vmc, matric_pot, Ksat, ):
    return delt_vmc * matric_pot /Ksat

def CN_est1(total_rfall, ro_coef_GA, lamb = .2):
    num = 2 * total_rfall * (1-ro_coef_GA)
    denom = ro_coef_GA *(1-lamb) + 2 * lamb + sqrt(ro_coef_GA) + sqrt(ro_coef_GA * (1-lamb)**2 + 4 * lamb)
    return 25400 / (254 + (num/denom))

def CN_est(t, Ksat, rho, norm_duration, ro_coef_GA, lamb = .2):
    '''
    

    Parameters
    ----------
    t : TYPE
        DESCRIPTION.
    Ksat : TYPE
        DESCRIPTION.
    rho : TYPE
        DESCRIPTION.
    norm_duration : TYPE
        DESCRIPTION.
    runoff_coef : TYPE
        DESCRIPTION.
    lamb : TYPE, optional
        DESCRIPTION. The default is .2.

    Returns
    -------
    TYPE
        DESCRIPTION.
        
        Equation #22

    '''
    
    num = t* Ksat * (1-ro_coef_GA) * rho * (rho * norm_duration +1 - lamb)
    denom = rho * norm_duration * (1+lamb) - lamb**2 
    
        
    
    return 25400 / (254 +(num/denom))


def calc_rho(Ksat, r_intens):
    return r_intens/Ksat

def calc_norm_duration(ro_coef_GA, rho, lamb = .2):
    num = ro_coef_GA * (1- lamb) + 2*lamb + sqrt(ro_coef_GA) * sqrt(ro_coef_GA * (1-lamb)+4*lamb)
    denom = 2*rho * (1-ro_coef_GA)
    return num/denom

def KsMax(GA_t_scale, i, t):
    
    return (i * (sqrt(t)* sqrt(t+ 4* GA_t_scale)-t )) / 2*GA_t_scale
    
#%%
def CN_min(Ksat_max, hrs, rho, lamb =.2):
    return 25400 / \
        (254 + (hrs *Ksat_max * rho)/ (lamb * (1+lamb)- lamb**2))
        
        
def estimate_CN(delt_vmc, Ksat, matric_pot, hrs, r_intens, ro_coef_GA, lamb = .2):
    rho = calc_rho(Ksat, r_intens)
    norm_duration = calc_norm_duration(ro_coef_GA, rho, lamb)
    return CN_est(hrs, Ksat, rho, norm_duration, ro_coef_GA, lamb )


    
    