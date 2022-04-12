# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 09:37:24 2022

@author: benja
"""


import math
import numpy as np
from math import exp, log, nan
from rosetta import rosetta, SoilData
from numbers import Number
from functools import partial

def rosetta3(perc_sand, perc_silt, perc_clay, bulk_dens = None, field_capacity_wc = None, pwp_capacity = None):
    '''Run the rosetta 3 model to retrieve predicted soil parameters.

    Parameters
    ----------
    perc_clay : numeric
        % clay content of the soil
    perc_sand : TYPE
        %  sand content of the soil
    perc_silt : TYPE
        %  silt content of the soil
    bulk_dens : numeric
        bulk density of soil (kg/m^3)

    Returns
    -------
    Array of residual water content (fraction), 
    saturated water content (fractional), 
    two water retention parameters and 
    the saturated hydraulic conductivity (mm/hr).

    '''
    
    vars_ = [v for v in locals().values() if v is not None]
    
    if isinstance(perc_clay, Number):
        for i, v in enumerate(vars_):
            vars_[i] = [v]
    
        single = True
    
    
    
    
    else:
        single = False
    try:
        data = SoilData.from_array(np.array(vars_).T)
    except:
        print(vars_)
        raise
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



def VG_swr(phi,theta_r, theta_s, npar, phi_entry, phi_entry_is_inv = True ): 
    '''
    Implement the van Genuchten soil water retention model to predict
    volumetric soil water content at a pressure level si. 
    
    For details on the model itself: 
    https://hwbdocuments.env.nm.gov/Los%20Alamos%20National%20Labs/TA%2054/11569.pdf
    
    
    Modified from https://github.com/tpclement/PySWR
     
    
    Parameters
    ----------
    si : numeric
        suction pressure to estimate soil water content at (kPa)
    teta_r : numeric
        residual water content (cm^3/cm^3).
    teta_s : TYPE
        saturated water content (cm^3/cm^3).
    npar : numeric
        pore-size distribution (dimensionless). Can be calculated from the rosetta model
    

    phi_entry: air entry at suction or inverse air entry at suction.
    inverse of air entry at suction Can be calculated from the rosetta model
    
    phi_entry_is_inv: Flag for whether phi_entry is the inverse or actual value.
    default is True.
    
    Returns
    -------
    teta : Volumetric Soil water content, value between 0 and 1

    '''
    
    #phi  = phi * 10.19716
    if not phi_entry_is_inv:
        alpha = 1/phi_entry
    
    m = 1-(1/npar)
    theta = theta_r + ((theta_s-theta_r)/(1+(alpha*phi)**npar)**m)
    return theta
#%%
fc_VG = partial(VG_swr, 33)
pwp_VG = partial(VG_swr, 1500)

def rosetta_solver(perc_sand, perc_silt, perc_clay, bulk_dens, OM):
    '''
    Combine the rosetta model, the van Genuchten soil water retention model
    and the balland pedotransfer functions to estimate field capacity, 
    permanent wilting point, saturated hydraulic conductivity and
    saturated water content.

    Parameters
    ----------
    perc_sand : fraction sand content of soil. between 0 and 1
        
    perc_silt : fraction silt content of soil. between 0 and 1
        
    perc_clay : fraction clay content of soil. between 0 and 1
        
    bulk_dens : bulk density of the soil in g/cm**3
        
    OM : organic matter % of soil. 

    Returns
    -------
    fc_avg : numeric
        volumetric water content at field capacity. between 0 and 1
    pwp_avg : numeric
        volumetric water content at the permanent wilting point. 
        between 0 and 1
    ksat : saturated hydraulic conductivity in mm/hr.

    '''
    
    if  isinstance(perc_sand, Number):
        prep_func = lambda x: x
    else:
        prep_func = np.vectorize
    fieldCapacity = prep_func(fieldCapacity_bal)
    permWiltPoint = prep_func(permWiltPoint_bal)
    fcVG = prep_func(fc_VG)
    pwpVG = prep_func(pwp_VG)
    

    
    o = rosetta3(perc_sand, perc_silt, perc_clay, bulk_dens)
    theta_s = o[1]
    theta_r = o[0]
    npar = o[3]
    alpha = o[2]
    ksat = o[4]
    
    for i in range(5):
        fc_bal = fieldCapacity(theta_s, bulk_dens, perc_clay, perc_sand, OM) 
        pwp_bal = permWiltPoint(fc_bal, perc_clay, perc_sand, OM, bulk_dens)
        #print(fc_bal)
        o = rosetta3(perc_sand, perc_silt, perc_clay, bulk_dens, fc_bal, pwp_bal)
        theta_s = o[1]
        theta_r = o[0]
        npar = o[3]
        alpha = o[2]
    fc_vg = fcVG(theta_r, theta_s, npar, alpha= alpha)
    pwp_vg = pwpVG(theta_r, theta_s, npar, alpha = alpha)
    fc_avg = fc_bal
    pwp_avg = pwp_bal
    #fc_avg = (fc_bal + fc_vg)/2
    #pwp_avg = (pwp_bal + pwp_vg)/2
    #o = rosetta3(perc_sand, perc_silt, perc_clay, bulk_dens, fc_avg, pwp_avg)
    #theta_s = o[1]
    return   fc_avg, pwp_avg, ksat, theta_s

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


def saturatedWC_byweight(bulk_dens, sat_wc):
    return sat_wc / bulk_dens
    
def fieldCapacity(perc_sand, OM, saturated_wc, bulk_dens):
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
    p_sand = perc_sand/100
    OM = OM/100
    saturated_wc_by_weight = saturated_wc/bulk_dens
    field_capacity_wc = 1- exp((-.588* (1-p_sand) -(1.73 * OM))/saturated_wc_by_weight)
    return field_capacity_wc * bulk_dens


def fieldCapacity_bal(saturated_wc, bulk_dens, perc_clay, perc_sand, OM):
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
    bulk_dens:

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    
    
    OM = OM/100
    p_clay = perc_clay/100
    p_sand = perc_sand/100
    sat_wc_wt = saturated_wc /bulk_dens
    fc = sat_wc_wt * (.565 + .426 * p_clay**.5) * exp((.103 * p_sand - .785 *(OM))/sat_wc_wt* -1)
    field_capacity = fc * bulk_dens
    return field_capacity
    

def permWiltPoint_bal(field_capacity, perc_clay, perc_sand, OM, bulk_dens):
    '''Equation to calculate water content of soil at permenant wilting point in %.
    From doi:10.1016/j.ecolmodel.2008.07.009, page  7

    Parameters
    ----------
    perc_sand : TYPE
        DESCRIPTION.
    OM : TYPE
        DESCRIPTION.
    saturated_wc : TYPE
        DESCRIPTION.
    bulk_dens: numeric
    bulk density in cm^3 /cm^3

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    p_clay = perc_clay/100
    OM = OM/100
    fc_wt = field_capacity /bulk_dens
    pwpwt = fc_wt * (.17 + .662 * p_clay**.5) * exp((0 -(-1.4 * OM)) /fc_wt * -1 )
    pwp = pwpwt * bulk_dens
    return pwp
    

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
 W. J. Rawls, ; D. L. Brakensiek, ; K. E. Saxtonn, (1982). Estimation of Soil Water Properties. Transactions of the ASAE, 25(5), 13161320. doi:10.13031/2013.33720 
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
    Predicts the field capacity water content in cm^3/cm^3
    One of two regression equations from:
    W. J. Rawls, ; D. L. Brakensiek, ; K. E. Saxtonn, (1982). Estimation of Soil Water Properties. Transactions of the ASAE, 25(5), 13161320. doi:10.13031/2013.33720 

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
    From doi:10.1016/j.ecolmodel.2008.07.009, page  8 (equation 9)
    May want a better method

    Parameters
    ----------
    perc_clay : TYPE
        DESCRIPTION.
    OM : TYPE
        DESCRIPTION.
    field_capacity : numeric
        field capacity in proportion moisture

    Returns
    -------
    permenant wilting point water content
        in m**3/m**3

    '''
    field_capacity = field_capacity *100
    pwp_capacity =field_capacity * (1- exp((-.511 * perc_clay - .865 * OM)/field_capacity)) 
    
    return pwp_capacity/100

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


