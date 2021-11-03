#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 12:09:39 2021

@author: bdube
"""
import math
import numpy as np
from math import exp, log, nan
from data_utils import soil_tex_pyramid, soil_data
from rosetta import rosetta, SoilData
from numbers import Number
import pandas as pd
import paths
from numpy import vectorize
from sympy.solvers import solve
from sympy import Symbol, Eq
from functools import partial
from SCS_GA_links import KsMax, GA_t_scale, CN_min, calc_norm_duration, calc_rho, CN_est, estimate_CN
#%%

def calc_LS(slope, slope_length, constant = 22.1):
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

def soil_perm_class(Ksat):
    '''
    Soil Permeability class for wishmeir K
    page 255 swat technical docs.

    Parameters
    ----------
    Ksat : TYPE
        DESCRIPTION.

    Returns
    -------
    int
        DESCRIPTION.

    '''
    
    if Ksat> 150:
        return 1
    elif Ksat>50:
        return 2
    elif Ksat>15:
        return 3
    elif Ksat>5:
        return 4
    elif Ksat>1:
        return 5
    return 6


def run_on_df(func):
    def wrapper(df, columns):
        return vectorize(func)(*[df[col] for col in columns])
    return wrapper

def run_df(func):
    def wrapper(df, standard_args = {}):
        n_args = func.__code__.co_argcount
        columns = func.__code__.co_varnames[:n_args]
        
        if standard_args:
            columns = [c for c in columns if c not in standard_args]
            f = partial(func, **standard_args)
            try:
                return vectorize(f)(**{col: df[col] for col in columns})
            except:
                print(columns)
                #print(standard_args)
                print(f)
                raise
                
        return vectorize(func)(*[df[col] for col in columns])
    
    return wrapper

def rosetta3(perc_clay, perc_sand, perc_silt, bulk_dens):
    '''
    

    Parameters
    ----------
    perc_clay : TYPE
        DESCRIPTION.
    perc_sand : TYPE
        DESCRIPTION.
    perc_silt : TYPE
        DESCRIPTION.
    bulk_dens : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    if isinstance(perc_clay, Number):
        single = True
        perc_clay = [perc_clay]
        perc_sand = [perc_sand]
        perc_silt = [perc_silt]
        bulk_dens = [bulk_dens]
    else:
        single = False
    data = SoilData.from_array(np.array([perc_sand,perc_silt, perc_clay, bulk_dens]).T)
    means, std, codes = rosetta(3, data)
    
    residual_h2o_cont = means[:,0]
    sat_h2o_cont = means[:,1]
    ret_param1 = np.power(10, means[:,2])
    ret_param2 = np.power(10, means[:,3])
    Ksat = np.power(10, means[:,4]) * 10/24 #convert cm/d to mm/hr
    if single:
        out_array =  np.array([a[0] for a in [residual_h2o_cont, sat_h2o_cont, ret_param1, ret_param2, Ksat]])
    else:
        out_array = np.array([residual_h2o_cont, sat_h2o_cont, ret_param1, ret_param2, Ksat])
    return  out_array



#%%
def drainExcessSW(soil_water, field_capacity, TT_perc, timestep = 1):
    '''
    SWAT docs 151

    Parameters
    ----------
    soil_water : TYPE
        DESCRIPTION.
    field_capacity : TYPE
        DESCRIPTION.
    TT_perc : TYPE
        DESCRIPTION.
    timestep : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    None.
    
    

    '''
    if soil_water<field_capacity:
        return 0
    
    else:
        return soil_water - field_capacity * (1-exp(*-1 * timestep/TT_perc))



def cacl_TTperc(saturated_wc, field_capacity, Ksat):
    '''
    SWAT docs page 151

    Parameters
    ----------
    saturated_wc : TYPE
        DESCRIPTION.
    field_capacity_wc : TYPE
        DESCRIPTION.
    Ksat : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    pass

def relativeFieldCapacity(saturated_wc, field_capacity_wc):
    '''
    From Reynolds, et al 2009
    doi:10.1016/j.geoderma.2009.06.009
    

    Parameters
    ----------
    saturated_wc : TYPE
        DESCRIPTION.
    field_capacity_wc : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    
    rel_field_capac =  field_capacity_wc / saturated_wc
    return rel_field_capac


def airCapacity(saturated_wc, field_capacity_wc):
    '''
    From Reynolds, et al 2009
    doi:10.1016/j.geoderma.2009.06.009

    Parameters
    ----------
    saturated_wc : TYPE
        DESCRIPTION.
    field_capacity_wc : TYPE
        DESCRIPTION.

    Returns
    -------
    air_capacity : TYPE
        DESCRIPTION.

    '''
    
    air_capacity = saturated_wc - field_capacity_wc
    return air_capacity

def strucStabilityIndex(perc_silt, perc_clay, orgC):
    '''
    From Reynolds, et al 2009
    doi:10.1016/j.geoderma.2009.06.009

    Parameters
    ----------
    perc_silt : TYPE
        DESCRIPTION.
    perc_clay : TYPE
        DESCRIPTION.
    orgC : TYPE
        DESCRIPTION.

    Returns
    -------
    SI : TYPE
        DESCRIPTION.

    '''
    
    SI = 1.724 * orgC / (perc_silt + perc_clay)
    
    return SI

    
def fieldCapacity(perc_sand, OM, saturated_wc):
    '''Equation to calculate water content of soil at field capacity in %.
    From doi:10.1016/j.ecolmodel.2008.07.009, page  7

    Parameters
    ----------
    perc_sand : TYPE
        DESCRIPTION.
    OM : TYPE
        DESCRIPTION.
    saturated_wc : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    field_capacity_wc = -.33 + .32 * perc_sand - 1.53 * OM + saturated_wc
    return field_capacity_wc


def fieldCapacityRawls1(perc_sand,  OM, pwp_capacity):
    '''
    Calculate the water content % at field capacity using the regression equation reported by Rawls (1982)

    Parameters
    ----------
    perc_sand : TYPE
        DESCRIPTION.
    OM : TYPE
        DESCRIPTION.
    pwp_capacity : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    field_capacity_wc = .2391 -.0019 * perc_sand +.021 * OM + pwp_capacity * .72 
    return field_capacity_wc


def permwiltPointRawls(perc_clay, OM):
    '''
    Caculate the pwp water content using the regression equation reported by Rawls (1982)

    Parameters
    ----------
    perc_clay : int
        DESCRIPTION.
    OM : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    pwp_capacity = .0260  + perc_clay * .005 + .0158 * OM 
    return pwp_capacity



def fieldCapacityRawls2(perc_sand, perc_clay, OM):
    '''
    

    Parameters
    ----------
    perc_sand : TYPE
        DESCRIPTION.
    perc_clay : TYPE
        DESCRIPTION.
    OM : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    field_capacity_wc = .2576 - perc_sand * .002 + perc_clay * .0036 + OM * .0299
    return field_capacity_wc

def permWiltPoint(perc_clay, OM, field_capacity):
    '''Equation to calculate water content of soil at permanant wilting point in %.
    From doi:10.1016/j.ecolmodel.2008.07.009, page  8
    May want a better method

    Parameters
    ----------
    perc_clay : TYPE
        DESCRIPTION.
    OM : TYPE
        DESCRIPTION.
    field_capacity : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    
    pwp_capacity = .23 * perc_clay + .33 * OM + field_capacity * .18
    
    return pwp_capacity




def particle_density(perc_clay, OM):
    '''Predict soil Particle Density based on clay fraction and SOM.
    equation from https://doi.org/10.1016/j.geoderma.2016.10.020
    '''
    particle_dens = 2.652 + .216 * perc_clay/100 - 2.237 * OM/100
    return particle_dens


def calc_porosity(bulk_dens, particle_dens):
    '''Calculate porosity cm3/cm3 based on particle density and bulk density.'''
    porosity =  1- (bulk_dens/particle_dens)
    return porosity


def williamsK(perc_clay, perc_silt, perc_sand,  orgC):
    '''
    Williams equation for calculating the USLE K factor of soil

    Parameters
    ----------
    perc_clay : numeric
        percent clay particles in soil. Range: 0-100
    perc_silt : numeric
        percent silt particles in soil. Range: 0-100
    perc_sand : numeric
        percent sand particles in soil. Range: 0-100
    orgC : numeric
        % organic carbon of soil. Range: 0-100
    Returns
    -------
    float
    USLE K factor: soil erodibility. Range: 0-1

    '''
    
    
    csand = .2 + .3 * exp(-.0256 * perc_sand * (1- perc_silt/100))
    cl_si = (perc_silt/(perc_clay+perc_silt)) **.3
    org_c_f = 1 - (.25 * orgC)/(orgC+ exp(3.72- 2.95 * orgC))
    sandfac = (1-perc_sand/100)
    hisand = 1 - (.7 * sandfac) / (sandfac + exp(-5.51 + 22.9 * sandfac))
    k_factor_williams = np.product([csand, cl_si, org_c_f, hisand])
    return k_factor_williams

#williamsK(60, 20, 20, 3)#
#%%



    
    
#%%
def wishmeirK(perc_silt, perc_vfs, perc_clay, perm_code, soil_str, OM):
    '''
    Wishmeir equation for calculating USLE K factor.

    Parameters
    ----------
    perc_silt : TYPE
         percent silt particles in soil. Range: 0-100
    vfs : TYPE
        percent very fine sand particles in soil. Range: 0-100
    perc_clay : TYPE
        percent clay particles in soil. Range: 0-100
    perm_code : int
        soil permeability class, range:  1-6
    soil_str : int
        soil structure class: range 1-4
    OM : TYPE
        percent organic matter in soil. Range: 0-100
        

    Returns
    -------
    float
    USLE Soil erodibility factor K, between 0 and 1

    '''
    
    M = (perc_silt + perc_vfs) * (100-perc_clay)
    k_factor_wishmeir = ( .00021 * M**1.14 * (12-OM) + 3.25*soil_str + 2.5*(perm_code-3))/100
    return k_factor_wishmeir

wishmeirK(15, 2, 60, 3, 3, 2)


#%%


def rational_eq(Qsurf, Rday, i):
    '''page 110 of docs
    C: a coefficient
    i: rainfall intensity mm/hr
    
    note, original equation includes area, and is scaled to give 
    m**3 water based on km**2.
        '''
    
    Qpeak = (Qsurf/Rday)*i/3.6
    return Qpeak



def calc_tconc(lenslope, avg_slope, flow_rate, manning_n):
    '''Caculate time of concentration for a subbasin.
    Page 111
    lenslope: length of slope
    avg_slope: avg slope gradient
    flow_rate = rate of flow
    manning_n = mannings roughness coefficient
    '''
    
    
    tconc = lenslope **.6 * manning_n **.6 / (18 *avg_slope**.3)
    return tconc


def modified_rational(frac_rainfall, Qsurf, tconc):
    '''Peak runoff rate. Page 115
    Original equation yields m**3 per km**2, 
    we divide by 100 to yield M**3 per ha
    '''
    return frac_rainfall * Qsurf/ (3.6 * tconc)


def calc_frac_rainfall(tconc, max_hlf_hr):
    '''Calculate fraction_rainfall falling in 1/2 hr max'''
    frac_rainfall = 1 - exp(2 * tconc * math.log(1-max_hlf_hr))
    return frac_rainfall

def peak_runoff(Qsurf, tconc, max_hlf_hr):
    frac_rainfall = calc_frac_rainfall(tconc, max_hlf_hr)
    Qpeak = modified_rational(frac_rainfall, Qsurf,  tconc)
    return Qpeak




#%%


def SCN_alt(SW, w1, w2, Smax):
    '''
    Alternate method for curve number, based on 
    soil water content 

    Parameters
    ----------
    SW : numeric
        soil water content in mm
    w1 : TYPE
        shape parameter generated by calc_alt_params
    w2 : TYPE
        shape parameter generated by calc_alt_params
    Smax : TYPE
        Maximum value that S, the retention parameter, can take for this soil.

    Returns
    -------
    An SCS curve number between 0 and 100

    
    
    page 103'''
    S = Smax * (1- (SW/ (SW + exp(w1-w2 * SW))))
    SCNe =25400 /(S +254)
    return SCNe
    
def calc_alt_params(curve_num, field_cap_wat, sat_wat):
    '''Calculate the needed parameters for the alternate SCN equation.
    curve_number = the curve number for average soil conditions
    field_cap_wat = soil water at field capacity
    sat_wat = soil water at field capacity
    
    Page 104
    '''
    
    Smax = S_calc(CN1(curve_num))
    S3 = S_calc(CN3(curve_num))
    num = log(field_cap_wat/(1-S3*(1/Smax))-field_cap_wat) - log(sat_wat / (1-2.54 * (1/Smax)) - sat_wat)
    denom = sat_wat - field_cap_wat
    w2 = num/denom
    w1 = log(field_cap_wat/ (1-S3 * 1/Smax) - field_cap_wat ) + w2 * field_cap_wat
    
    return w1, w2, Smax

def S_calc(curve_num):
        
    S = 25.4 * (1000/curve_num -10)
    return S

def inv_Scalc(S):
    curve_num = 1000/(S/25.4 + 10)
    return curve_num


def effective_S_infl(infl_to_runoff):
    if infl_to_runoff ==0:
        return nan
    
    else:
        return 1000/(infl_to_runoff/10.8 + 10)


def effective_S(Qsurf, Rday):
    '''
    Estimate the effective S value from the SCN equation based on 
    the surface runoff and Rainfall
    

    Parameters
    ----------
    Qsurf : float
    surface runoff (mm)
    Rday : numeric
        rainfall for the day (mm).

    Returns
    -------
    TYPE
        DESCRIPTION.
    str
        DESCRIPTION.

    '''
    S = Symbol('S')
    if Qsurf == 0:
        Smin = Rday * 5
        return  Smin, 'min'
    else:
        S_vals = solve((Rday-(.2*S))**2/(Rday+(.8*S))-Qsurf, S)
        return min(S_vals), 'est'
    

def SCN_eq(Rday, curve_num, Ia = None):
    '''SCS curve number equation.
    Rday = rainfall for day (mm)
    Ia  = intial abstractions (mm H2O)
    CN = retention parameter
    Page 122
    '''

    S = S_calc(curve_num)
    if not Ia:
        Ia = .2*S
    if Rday > Ia:
        Qsurf = (Rday- Ia)**2/(Rday -Ia +S)
    else:
        Qsurf = 0
    
    return Qsurf
    
def CN1(CN2):
    cn_adj = CN2 - (20 * (100 - CN2)/(100-CN2 + exp(2.533 -.0636 * (100-CN2))))
    return cn_adj

def CN3(CN2):
    cn_adj = CN2 * exp(.00673*(100-CN2))
    return cn_adj


def infl_to_runoff(curve_num):
    S = S_calc(curve_num)
    return S *2
    
def CN(moisture_condition, CN2):
    '''Adjust Curve Number for moisture condition.
    page 103
    conditions: 
        1: dry
        2: average
        3: field capacity
    '''
    if moisture_condition ==2:
        SCNe = CN2
    elif moisture_condition == 1:
        SCNe = CN1(CN2)
    elif moisture_condition == 3:
        SCNe = CN3(CN2)
    return SCNe
    
def slope_adj_CN(CN2, slope):
    '''Adjust the standard CN2 for slope.
    Page 106
    '''
    
    C3 = CN3(CN2)
    C1 = CN1(CN2)
    SCNe = (C3 - C1)/3 * (1- 2* exp(-13.86 *slope)) + CN2 
    return SCNe

def musle(R_factor,  K, C, P , LS, CFRG):
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
   sed_yield = np.product([11.8, R_factor, K, C, P, LS, CFRG ])
   return sed_yield
   
def calc_CFRG(perc_coarse_frag):
    '''
    

    Parameters
    ----------
    perc_rock : float
        DESCRIPTION.

    Returns
    -------
    float
        DESCRIPTION.

    '''
    cfrg = exp(-0.053 * perc_coarse_frag)
    return cfrg




def greenAmpt_timestep(Ke, matricPot, delt_vmc, cum_inf, R_fall, time_step = 1, prev_inf = None):
    '''
    

    Parameters
    ----------
    Ke : float
        Effective hydraulic conductivity from calc_Keff
    matricPot : float
        soil matric potential. From calcMatricPot
    delt_vmc : float
        change in volumetric moisture content across the wetting front.
        From calc_deltVMC
    cum_inf : numeric
        cumulative infiltration that has occured thus far (mm)
    R_fall : numeric
        amount of rainfall in timestep (mm)
    time_step : numeric
        length of time-step in hours

    Returns
    -------
    new_inf_rate:
    float
        infiltration rate during the time-step (mm)
    runoff:
        float
        amount of runoff during time-step (mm)

    From: SWAT technical docs, page 108
    '''
    
    
    
    new_inf_rate = 0
    
    term_1 = Ke * time_step
    term_2 = matricPot * delt_vmc
    inf_rate = term_1
    while True:
        try:
            new_inf_rate = term_1 + term_2 * log(((inf_rate + cum_inf) + term_2) / (cum_inf + term_2 ))
        except ValueError:
            print(new_inf_rate)
            print(inf_rate)
            print(term_1)
            print(term_2)
            raise
        if math.isclose(inf_rate, new_inf_rate, rel_tol = .01):
            
            break
        else:
            inf_rate = new_inf_rate
            #print(inf_rate)
    #print(new_inf_rate)
    if new_inf_rate >R_fall:
        new_inf_rate = R_fall
        runoff = 0
    else:
         runoff = R_fall - new_inf_rate
    return new_inf_rate, runoff


#%%
def greenAmpt_timestep2(Ke, matricPot, delt_vmc, cum_inf, R_fall, time_step = 1, prev_inf = None):
    '''
    

    Parameters
    ----------
    Ke : float
        Effective hydraulic conductivity from calc_Keff
    matricPot : float
        soil matric potential. From calcMatricPot
    delt_vmc : float
        change in volumetric moisture content across the wetting front.
        From calc_deltVMC
    cum_inf : numeric
        cumulative infiltration that has occured thus far (mm)
    R_fall : numeric
        amount of rainfall in timestep (mm)
    time_step : numeric
        length of time-step in hours

    Returns
    -------
    new_inf_rate:
    float
        infiltration rate during the time-step (mm)
    runoff:
        float
        amount of runoff during time-step (mm)

    From: SWAT technical docs, page 108
    '''
    
    
    inf_rate = 3
    new_inf_rate = 0
    
    term_1 = Ke * time_step
    term_2 = matricPot * delt_vmc
    while True:
        try:
            new_inf_rate = term_1 + term_2* log(1 + ((cum_inf+inf_rate)/term_2)) - cum_inf
            #new_inf_rate = term_1 + term_2 * log(((inf_rate + cum_inf) + term_2) / (cum_inf + term_2 ))
        except ValueError:
            print(new_inf_rate)
            print(inf_rate)
            print(term_1)
            print(term_2)
            raise
        if math.isclose(inf_rate, new_inf_rate, rel_tol = .01):
            
            break
        else:
            inf_rate = new_inf_rate
            #print(inf_rate)
    #print(new_inf_rate)
    if new_inf_rate >R_fall:
        new_inf_rate = R_fall
        runoff = 0
    else:
         runoff = R_fall - new_inf_rate
    return new_inf_rate, runoff

#%%
def greenAmpt_timestep3(Ke, matricPot, delt_vmc, cum_inf, R_fall, time_step = 1, prev_inf = None):
    '''
    

    Parameters
    ----------
    Ke : float
        Effective hydraulic conductivity from calc_Keff
    matricPot : float
        soil matric potential. From calcMatricPot
    delt_vmc : float
        change in volumetric moisture content across the wetting front.
        From calc_deltVMC
    cum_inf : numeric
        cumulative infiltration that has occured thus far (mm)
    R_fall : numeric
        amount of rainfall in timestep (mm)
    time_step : numeric
        length of time-step in hours

    Returns
    -------
    new_inf_rate:
    float
        infiltration rate during the time-step (mm)
    runoff:
        float
        amount of runoff during time-step (mm)

    From: SWAT technical docs, page 108
    '''
    
    
   
    
    term_1 = Ke * time_step
    term_2 = matricPot * delt_vmc
    
    new_inf_rate = term_1 * ((term_2+1) /(cum_inf +1e-6))
    
    if new_inf_rate >R_fall:
        new_inf_rate = R_fall
        runoff = 0
    else:
         runoff = R_fall - new_inf_rate
    return new_inf_rate, runoff



def greenAmpt_timestep4(Ke, matricPot, delt_vmc, cum_inf, R_fall, time_step = 1, prev_inf = None):
    '''
    

    Parameters
    ----------
    Ke : float
        Effective hydraulic conductivity from calc_Keff
    matricPot : float
        soil matric potential. From calcMatricPot
    delt_vmc : float
        change in volumetric moisture content across the wetting front.
        From calc_deltVMC
    cum_inf : numeric
        cumulative infiltration that has occured thus far (mm)
    R_fall : numeric
        amount of rainfall in timestep (mm)
    time_step : numeric
        length of time-step in hours

    Returns
    -------
    new_inf_rate:
    float
        infiltration rate during the time-step (mm)
    runoff:
        float
        amount of runoff during time-step (mm)

    From: SWAT technical docs, page 108
    
    Code translated directly from surq_greenampt.f in the swat source code
    '''
    
    term_1 = Ke * time_step
    term_2 = matricPot * delt_vmc
    if prev_inf:
        tst = prev_inf
    else:
        tst = term_1 
    counter = 0
    while True:
        try:
            f1 = cum_inf +  term_1  +   term_2 * log((tst + term_2)/(cum_inf + term_2))
        except ValueError:
            print(locals())
            raise
        
        
        if (abs(f1 - tst) <= 0.001):
            break
        else:
            tst = f1
        counter +=1
        if counter>10000:
            raise ValueError()
    new_inf_rate = f1 -cum_inf
    if new_inf_rate >R_fall:
        new_inf_rate = R_fall
        runoff = 0
    else:
         runoff = R_fall - new_inf_rate
    return new_inf_rate, runoff


def greenAmpt_timestep5(Ksat, matricPot, delt_vmc, cum_inf, R_fall, time_step = 1, prev_inf = None):
    '''
    

    Parameters
    ----------
    Ke : float
        Effective hydraulic conductivity from calc_Keff
    matricPot : float
        soil matric potential. From calcMatricPot
    delt_vmc : float
        change in volumetric moisture content across the wetting front.
        From calc_deltVMC
    cum_inf : numeric
        cumulative infiltration that has occured thus far (mm)
    R_fall : numeric
        amount of rainfall in timestep (mm)
    time_step : numeric
        length of time-step in hours

    Returns
    -------
    new_inf_rate:
    float
        infiltration rate during the time-step (mm)
    runoff:
        float
        amount of runoff during time-step (mm)

    From: 
        https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/infiltration-and-runoff-volume/green-and-ampt-loss-model
    '''
    
    term_1 = Ksat * time_step
    term_2 = matricPot * delt_vmc
    if prev_inf:
        tst = prev_inf
    else:
        tst = term_1 
    while True:
        f1 = term_1 * (1 + term_2)/(cum_inf + tst + .0001)
                      
        
        
        
        if (abs(f1 - tst) <= 0.001):
            break
        else:
            tst = f1
            #print(f1)
    
    f1
    if f1 >R_fall:
        f1 = R_fall
        runoff = 0
    else:
         runoff = cum_inf - f1
    return f1, runoff



def greenAmpt_timestep6(Ksat, matricPot, delt_vmc, cum_inf, R_fall, time_step = 1, prev_inf = None):
    '''
    

    Parameters
    ----------
    Ke : float
        Effective hydraulic conductivity from calc_Keff
    matricPot : float
        soil matric potential. From calcMatricPot
    delt_vmc : float
        change in volumetric moisture content across the wetting front.
        From calc_deltVMC
    cum_inf : numeric
        cumulative infiltration that has occured thus far (mm)
    R_fall : numeric
        amount of rainfall in timestep (mm)
    time_step : numeric
        length of time-step in hours

    Returns
    -------
    new_inf_rate:
    float
        infiltration rate during the time-step (mm)
    runoff:
        float
        amount of runoff during time-step (mm)

    From: 
        https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/infiltration-and-runoff-volume/green-and-ampt-loss-model
    '''
    
    term_1 = Ksat * time_step
    term_2 = matricPot * delt_vmc
    if prev_inf:
        tst = prev_inf
    else:
        tst = term_1 
    while True:
        f1 = term_1 * (1 + term_2)/(cum_inf + tst + .0001)
                      
        
        
        
        if (abs(f1 - tst) <= 0.001):
            break
        else:
            tst = f1
            #print(f1)
    
    f1
    if f1 >R_fall:
        f1 = R_fall
        runoff = 0
    else:
         runoff = cum_inf - f1
    return f1, runoff

#%%

def greenAmpt(Ke, matricPot, delt_vmc, cumulative_infl):
    '''
    

    Parameters
    ----------
    Ke : float
        hydraulic conductivity.
    matricPot : float
       wetting front matric potential (mm)
    delt_vmc : float
        change in voulumetric moisture content across wetting front
    cumulative_infl : float
        cumulative infiltration at time t (mm H2O)

    Returns
    -------
    float
        DESCRIPTION.

    
    
    
    
    Page 107
   
    '''
    Qsurf = Ke *(1 + (matricPot * delt_vmc/cumulative_infl))
    return Qsurf

def calcKeff(curve_num, Ksat):
    '''
    Calculate the effective hydraulic conductivity for a soil

    Parameters
    ----------
    curve_num : numeric
        SCS curve number for soil
    Ksat : numeric
        Saturated Hydraulic conductivity for soil (mm/hr)

    Returns
    -------
    float
    effective hydraulic conductivity in mm/hr

    
    
    Swat technical docs Page 108
    '''
    
    
    
    Ke = (56.82 * (Ksat)**.286) /(1+(.051* exp(.062*curve_num))) -2
    if Ke > .5 * Ksat:
        return .5 * Ksat
    
    return Ke

def calcMatricPot(porosity, perc_clay, perc_sand):
    '''Calculate matric potential of a soil
    

    Parameters
    ----------
    porosity : float
        total porosity of a soil, mm/mm.
        Value between 0 and 1
    perc_clay : numeric
        % clay content of the top layer of soil.
        Between 0 and 100
    perc_sand : numeric 
        % sand content of the top layer of soil.

    Returns
    -------
    float 
    matric potential of a soil (pascals?)

    '''
    matricPot = 10 * exp(
        6.5309 - 7.32561 * porosity + 0.001583 * perc_clay**2 + 3.809479 * porosity**2 + \
        .000344 * perc_sand* perc_clay - .049837 * perc_sand**2 * porosity**2 + \
        .001608 * perc_sand**2 * porosity**2 + .001602 * perc_clay**2 * porosity **2 -\
        .0000136 * perc_sand**2 * perc_clay - .003479 * perc_clay**2 * porosity - \
        .000799 * perc_sand **2 * porosity        
        )
    return matricPot
        
        
def calc_deltVMC(swc, field_capacity, pwp_capacity, porosity):
    '''
    Change in volumetric moisture content across the wetting front

    Parameters
    ----------
    swc : numeric
        current soil water content (mm)
    field_capacity : numeric
        water content at field capacity (mm)
    porosity : soil porosity (mm3/mm3)

    Returns
    -------
    delt_vmc
        Change in volumetric moisture content across the wetting front (mm/mm)

    '''
    
    if field_capacity< pwp_capacity or field_capacity<0:
        return np.nan
    if swc >= field_capacity:
        return calc_deltVMC(field_capacity*.999, field_capacity, pwp_capacity, porosity)
    else:
        swc = swc - pwp_capacity
        field_capacity = field_capacity- pwp_capacity
        delt_vmc = (1-(swc/field_capacity)) * (.95 * porosity)
    return delt_vmc

calc_deltVMC(9.5, 9.0, 9.8, .58)

def run_dayCN(Rday, curve_num, rainfall_i, area, musle_params):
    '''Run one day using the SCS Curve number method
    

    Parameters
    ----------
    Rday : numeric
        Rainfall for the day (mm)
    curve_num : numeric
        SCS curve numer for the soil
    rainfall_i : numeric
        maximum 30-minute rainfall intensity (mm/hr)    
    
    area : numeric
        area the field in hectares
    musle_params : dict
        a dictionary of MUSLE parameters that remain constant

    Returns
    -------
    Qsurf : numeric
        surface runoff from the field for the day (cm)
    sed_yield : numeric
        sediment yield from the field for the day (metric tons)

    '''
    Qsurf = SCN_eq(Rday, curve_num)
    Qpeak = rational_eq(Qsurf, Rday, rainfall_i)
    R_factor = math.pow(Qsurf * Qpeak * area, .56)
    
    sed_yield = musle(R_factor, area, **musle_params)
    
    return Qsurf, sed_yield


def musle_R(Qsurf, Rday, rainfall_i, area):
    Qpeak = rational_eq(Qsurf, Rday, rainfall_i)
    R_factor = math.pow(Qsurf * Qpeak * area, .56)
    return R_factor


def test_1():
        area = 10
        musle_params = dict(LS = 2, P = .8, C = .2, CFRG = .9)
        #tests / examples for using green-ampt equation for runoff
        for k, v in soil_data.items():
            curve_num = 70
            tex = k.lower()
            perc_clay, perc_silt, perc_sand = soil_tex_pyramid.get(tex)
            orgC = perc_clay/30
            K_factor = williamsK(perc_clay, perc_silt, perc_sand, orgC)
            Ksat = float(v['Saturated Hydrolic Conductivity']) *10
            porosity = float(v['Effective Porosity/Saturation'] ) 
            swc = 40
            fc_wat_cont = 150
            matricPot = calcMatricPot(porosity, perc_clay, perc_sand)
            delt_vmc = calc_deltVMC(swc, fc_wat_cont,  porosity)
            Ke = calcKeff(curve_num, Ksat)
            time_step = 1/24 
            cum_inf = 0
            runoffs = []
            
            R_fall = 5
            print(k)
            for i in range(20):
                inf, runoff = greenAmpt_timestep(Ke, matricPot, delt_vmc, cum_inf, R_fall, time_step)
                cum_inf += inf
                runoffs.append(runoff)
                if i%4 == 0:
                    pass
                    #print(inf, runoff)
            cum_runoff = sum(runoffs)
            max_runoff = max(runoffs)
            R_factor = np.power(cum_runoff * max_runoff * area, .56)
            sed_yield = musle(R_factor, K_factor, **musle_params) 
            print(cum_inf, cum_runoff, sed_yield)
    
    #%%
        results = []
        for drain_val in range(1,5):
            line = []
            for soil_str in range(1, 7):
                line.append(wishmeirK(30, 30, 20, drain_val, soil_str, 1))
            results.append(line)
        
        #%%
        results = []
        for drain_val in range(1,5):
            line = []
            for soil_str in range(1, 7):
                line.append(wishmeirK(30, 30, 20, drain_val, soil_str, 2))
            results.append(line)


#%%
def sim_storm(total_rfall, hrs, Ke, matric_pot, swc, field_capacity, pwp_capacity, porosity, time_step = 1, method = greenAmpt_timestep5):
    if np.any(np.isnan([Ke, matric_pot, swc, field_capacity, porosity])):
        return np.nan, np.nan
   
    elif field_capacity<pwp_capacity:
        return np.nan, np.nan
    
    elif any([x<0 for x in [Ke, matric_pot, swc, field_capacity, porosity, pwp_capacity]]):
        return np.nan, np.nan
    
    timestep_rfall = total_rfall / hrs *time_step
    cum_inf = 0
    inf = None
    infs = []
    
     
    
    for i in range(int(hrs/(time_step))):
        if swc<0:
            return np.nan, np.nan
        delt_vmc = calc_deltVMC(swc, field_capacity, pwp_capacity, porosity)
        try:
            inf, runoff = method(Ke, matric_pot, delt_vmc, cum_inf, 
                                         timestep_rfall, time_step, inf)
        except ValueError:
            print(f'{locals()}')
            raise
        cum_inf += inf
        swc+=inf
        infs.append(inf)
        
    return total_rfall-cum_inf, max(infs) * time_step/4



#to do: 
#link in soil_hydro_groups
#find a way to link residual water capacity and saturated water capacity 
# with field capacity and wilting point
    
#%%

rename_cols = {**{f'{frac}_tot_psa': f'perc_{frac}' for frac in 
                  ['clay', 'sand', 'silt']},
**{'oc': 'orgC', 'db_13b': 'bulk_dens', 'wf_0175': 'perc_coarse_frag',
   'sand_vf_psa': 'perc_vfs',
   
   }}


def simulate_rfalls(hrs_list, rates, data):
    runoffs = []
    Seffs = []
    CNes = []
    CNe1 = []
    CN_mins = []
    for hrs in hrs_list:
        for rate in rates:
            print(hrs, rates)
            total_rfall = hrs * rate
            #print(total_rfall)
            rainfall_i = total_rfall / (hrs*2)
            data['swc'] = (data['field_capacity'] - data['pwp_capacity']) *.4 + data['pwp_capacity']
        
            Qsurf, max_runoff = run_df(sim_storm)(data, dict(total_rfall = total_rfall, hrs = hrs,  time_step =1, method =greenAmpt_timestep4))
            #print(Qiurf)
            data['Qsurf'] = Qsurf
            data['max_runoff'] = max_runoff
            
            data['R_factor'] = run_df(musle_R)(data, dict(Rday=total_rfall, rainfall_i = rainfall_i,
                                                              area = area/10000))
               
            
            data['sed_yield'] = run_df( musle)(data, dict(C = C, P=P, LS = LS ))
            
            data['delt_vmc'] = run_df(calc_deltVMC)(data)
            data['ro_coef_GA'] = data['Qsurf'] / total_rfall
            data['CNe'] = run_df(estimate_CN)(data, dict(r_intens = rainfall_i *2, hrs =hrs, lamb = .2))
            data['rho'] = run_df(calc_rho)(data, dict(r_intens = rainfall_i*2))
            data['GA_t_scale'] = run_df(GA_t_scale)(data)
            data['Ksat_max'] = run_df(KsMax)(data, dict(i = rainfall_i*2, t = hrs))
            data['CN_min'] = run_df(CN_min)(data, dict(hrs = hrs, lamb =.2 ))
            runoffs.append(data['Qsurf'].tolist())
            
            CNes.append(data['CNe'].tolist())
            CN_mins.append(data['CN_min'].tolist())
    return runoffs, CNes, CN_mins

def fill_bulk_dens(perc_clay, perc_silt, perc_sand):
    '''Fill in Estimated Bulk Density for a soil based on its
    soil textural class.'''
    
#%%            
if __name__ == '__main__':
    results = []
    mul = .7
    dat_from_descrips = pd.read_csv(paths.soil_descrips, index_col= 0).dropna(subset = ['soil_str'])
    df = pd.read_csv(paths.soil_char_data).set_index(['samp_name'])
    df.rename(columns = rename_cols, inplace =True)
    #series_names = ['Vergennes', 'Hadley', 'Marlow', 'Limerick', 'Adams']
    #dat = df.loc[series_names]
    dat = df.loc['Walpole']
    #dat.index = dat.index.levels[0]
    dat = dat.merge(dat_from_descrips[[ 'soil_drain', 'soil_str']], 
                        left_index = True, right_index = True, how = 'inner' )
    
    
    
    dat.dropna(subset = ['orgC'], inplace= True)
    dat['orig_OrgC'] = dat['orgC']
    dat['orig_depth'] = dat['depth']
    dat['orig_bulk_dens'] = dat['bulk_dens']
    for i in range(15):
        dat['bulk_dens'] = dat['orig_bulk_dens'] * (1/mul)**.1
        dat['depth'] = dat['orig_depth'] / (dat['bulk_dens']/dat['orig_bulk_dens'])
        dat['orgC'] = dat['orig_OrgC'] * mul 
        dat['OM'] = dat['orgC'] *1.72
        res = rosetta3(dat['perc_clay'], dat['perc_sand'], dat['perc_silt'], dat['bulk_dens'])
        dat = dat.assign(residual_h2o = res[0], saturated_wc = res[1]* 100 , ret_param1 = res[2], 
                     ret_param2 = res[3], Ksat = res[4]) 
        dat['perm_code'] = dat['Ksat'].apply(soil_perm_class)
        dat['particle_dens'] = run_df(particle_density)(dat)
        dat['porosity'] = run_df(calc_porosity)(dat)
        dat['matric_pot'] = run_df(calcMatricPot)(dat)
        dat['will_K'] = run_df(williamsK)(dat )
        dat['wishm_K'] = run_df(wishmeirK)(dat)
        dat['field_capacity'] = run_df(fieldCapacity)(dat)
        dat['pwp_capacity'] = run_df(permWiltPoint)(dat)
        dat['pwp_capacity'] = dat['pwp_capacity'] *  dat['depth'] /10
        dat['field_capacity'] = dat['field_capacity'] * dat['depth'] /10
        dat['CFRG'] = run_df(calc_CFRG)(dat)
        
        #dat['Ke'] = dat['Ksat'] /2
        
        
        results.append( dat[['horizon_key', 'OM', 'bulk_dens', 'will_K', 'wishm_K', 
                             'field_capacity', 'pwp_capacity', 'Ksat', 'matric_pot', 'porosity', 
                             'CFRG']])
        mul +=.1
    
    df_sim = pd.concat([r[r['horizon_key']== 'A'] for r in results], axis = 0)
    df_sim['curve_num'] = 75
    df_sim['Ke'] = run_df(calcKeff)(df_sim)
    df_sim['K'] = df_sim['wishm_K']
    df_sim.set_index('OM', append = True, inplace = True)
    #df_sim['swc'] = df_sim['field_capacity'] * .75
    #df_sim['runoff'] = run_df(sim_storm)(df_sim, dict(total_rfall = 25.4 *2, hrs = 10, time_step = 1))
    
    
        
    
    
    
#%%






    
    field_shape = (88, 2200)
    slope = 5
    area = field_shape[0] * field_shape[1]
    LS = calc_LS(slope, field_shape[0])
    P = .8; C = .2
    
    
    
    
    
        
    
    hrs = 3


    
hrs_list = [1,2,3,4]; rates = [5,10,15, 20, 25]
runoffs, CNes, CN_mins = simulate_rfalls(hrs_list, rates, df_sim)

rdf = pd.DataFrame(runoffs).T
rdf.index = df_sim.index
CN_df = pd.DataFrame(CNes).T
CN_df.index = df_sim.index

results = []
series_names = df_sim.index.levels[0].unique()
for series_name in series_names:
    sub_r = rdf.loc[series_name]
    sub_CN = CN_df.loc[series_name]
    keep_cols = sub_r.min(axis = 0) > 0 
    
    out = sub_CN.max(axis = 0) - sub_CN.min(axis =0)
    out.iloc[~keep_cols] = np.nan
    results.append(out.to_list()  )

impacts_df = pd.DataFrame(results).set_index(series_names)

impacts_df = impacts_df[(~impacts_df.isna()).sum(axis=1) >0]
impacts_df.to_csv('results/CNimpact.csv')
CN_df.to_csv('results/CNs.csv')
rdf.to_csv('results/runoffs.csv')
df_sim.to_csv('results/dfsim.csv')

#df_sim['Seff'], _ = run_df(effective_S)(df_sim, dict(Rday = total_rfall))
#df_sim['norm_duration'] = run_df(calc_norm_duration)(df_sim, dict(lamb = .2))
#df_sim['CN_est'] = run_df(CN_est)(df_sim, dict(t = hrs, lamb =.2))


#%%
import matplotlib.pyplot as plt
      
fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('Soil OM')
ax1.set_ylabel('Simulated Runoff (mm)', color=color)
ax1.plot(df_sim['OM'], df_sim['Qsurf'] , color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim([0, max(df_sim['Qsurf'])*1.2])
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Erosion, T/ha', color=color)  # we already handled the x-label with ax1
ax2.plot(df_sim['OM'], df_sim['sed_yield'], color=color)
ax2.set_ylim([0, max(df_sim['sed_yield'])*1.2])
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()
plt.savefig('figures/Vergennes.png')
plt.show()
