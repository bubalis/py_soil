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

def rosetta3(perc_sand, perc_silt, perc_clay, bulk_dens = None, field_capacity_wc = None, pwp_wc = None):
    '''Run the rosetta 3 model to retrieve predicted soil parameters.

    Parameters
    ----------
   
    perc_sand : numeric
        percent sand particles in soil. Range: 0-100
        
    perc_silt : numeric
        percent silt particles in soil. Range: 0-100
    
    perc_clay : numeric
        percent clay particles in soil. Range: 0-100
        
    
    bulk_dens : numeric
        bulk density of soil (kg/m^3)

    field_capacity_wc: volumetric water content of soil at field capacity.
   
    pwp_wc: volumetric water content of soil at the permanent wilting point.
        
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



def kSat_Nemes_Rawls(perc_sand, perc_clay, bulk_dens, OM):
    '''
    

    Parameters
    ----------
    
    perc_sand : numeric
        percent sand particles in soil. Range: 0-100
        

    perc_clay : numeric
        percent clay particles in soil. Range: 0-100
        
    
    bulk_dens : numeric
        bulk density of soil (kg/m^3)
   
   
    OM : TYPE
        Soil Organic Matter content. Range: 0-100

    Returns
    -------
    Ksat : saturated hydraulic conductivity in mm/hr

    '''
    
    
    x1 = -3.663 + 0.046 * perc_sand  #
    
    x2 =  -0.887 + 0.083 * perc_clay #
    x3 =  -9.699 + 6.451 * bulk_dens #
    
    x4 = -0.807 + 1.263 * OM #
    
    z1 = -0.428 + 0.998*x1 + 0.651*x1**2 + 0.130*x1**3
    z2 = 0.506*x1 - 0.188*x2 - 0.327*x3 - 0.094*x4
    
    z3 = -0.268 + 0.885*z1 + 0.544*z1**2 - 0.682*z1**3 +\
    0.320*z2 - 0.134*z1*z2 + 1.119*z1**2*z2 + 0.050*z2**2 -\
    0.645*z1*z2**2 + 0.160*z2**3 + 0.126*x4 - 0.144*z1*x4 -\
    0.372*z1**2*x4 + 0.247*z2*x4 + 0.795*z1*z2*x4 - 0.344*z2**2*x4 +\
    .038*x4**2 -  0.071*z1*x4**2 + 0.020*z2*x4**2 - 0.015*x4**3
    
    z4 = 0.102 + 1.383*z3 + 0.302*z3**2 + 0.103*z3**3 +\
    0.331*x2 + 0.693*z3*x2 + 0.541*z3**2*x2 + 0.198*x2**2 +\
    0.429*z3 *x2**2 + 0.092*x2**3 + 0.060*x3 + 0.277*z3*x3 +\
    0.417*z3**2*x3 + 0.242*x2*x3 + 0.929*z3*x2*x3 +\
    0.319*x2**2*x3 + 0.026*x3**2 + 0.094*z3*x3**2 + 0.116*x2*x3**2
    
    Ksat = 10**(0.571 + 0.956*z4)
    return Ksat


def kSat_Gupta(perc_sand, perc_clay, bulk_dens):
    '''
    

    Parameters
    ----------
    
    perc_sand : numeric
        percent sand particles in soil. Range: 0-100
        

    perc_clay : numeric
        percent clay particles in soil. Range: 0-100
        
    
    bulk_dens : numeric
        bulk density of soil (kg/m^3)
   


    Returns
    -------
    Ksat : saturated hydraulic conductivity in mm/hr

    '''
    
    b0 = 2.17 
    b1 = 0.9387 
    b2 = -0.8026 
    b3 = 0.0037 
    b4 = -0.017 
    b5 = 0.000015 
    b6 = 0.0025 
    b7 = 0.00086
    b8 = -0.00025 
    b9 = 0.000073 
    logKsat =  b0 + b1 * bulk_dens + b2 * bulk_dens**2 + b3 * perc_clay\
        + b4 * bulk_dens * perc_clay + b5 * perc_clay**2 + b6 * perc_sand\
        + b7 * bulk_dens * perc_sand + b8 * perc_clay * perc_sand + b9 * perc_sand**2
    return np.pow(10, logKsat)



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
    perc_sand : fraction sand content of soil. between 0 and 100
        
    perc_silt : fraction silt content of soil. between 0 and 100
        
    perc_clay : fraction clay content of soil. between 0 and 100
        
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





def saturatedWC_byweight(bulk_dens, sat_wc):
    return sat_wc / bulk_dens
    


def fieldCapacity_bal(saturated_wc, bulk_dens, perc_clay, perc_sand, OM):
    '''Equation to calculate volumetric water content of soil at field capacity
    From doi:10.1016/j.ecolmodel.2008.07.009, page  7

    Parameters
    ----------
    saturated_wc : volumetric water content of soil at saturation.
    bulk_dens : bulk density of the soil in g/cm**3 
    
    perc_clay : fraction clay content of soil. between 0 and 100
    perc_sand : fraction sand content of soil. between 0 and 100
    
    
    OM : organic matter % of soil. 
    
   
    

    Returns
    -------
    field_capacity_wc : volumetric water content at field capacity

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
    perc_sand : fraction sand content of soil. between 0 and 100
    
    OM : organic matter % of soil. 
    
    saturated_wc : volumetric water content at field capacity
    
    bulk_dens : bulk density of the soil in g/cm**3 

    Returns
    -------
    pwp_wc : volumetric water content at the permanent wilting point.

    '''
    
    
    p_clay = perc_clay/100
    OM = OM/100
    fc_wt = field_capacity /bulk_dens
    pwpwt = fc_wt * (.17 + .662 * p_clay**.5) * exp((0 -(-1.4 * OM)) /fc_wt * -1 )
    pwp = pwpwt * bulk_dens
    return pwp
    

def fieldCapacityRawls1(perc_sand,  OM, pwp_wc):
    '''
    Calculate the water content % at field capacity using the regression equation reported by Rawls (1982)

    Parameters
    ----------
    perc_sand : fraction sand content of soil. between 0 and 100
    
    OM : organic matter % of soil. 
    
    pwp_wc : volumetric water content at the permanent wilting point.

    Returns
    -------
    field_capacity_wc : volumetric water content at field capacity

    '''
    field_capacity_wc = .2391 -.0019 * perc_sand +.021 * OM + pwp_wc * .72 
    return field_capacity_wc


def permwiltPointRawls(perc_clay, OM):
    '''
    Caculate the pwp water content using the regression equation reported by Rawls (1982)
 W. J. Rawls, ; D. L. Brakensiek, ; K. E. Saxtonn, (1982). Estimation of Soil Water Properties. Transactions of the ASAE, 25(5), 13161320. doi:10.13031/2013.33720 
    Parameters
    ----------
    perc_clay : fraction clay content of soil. between 0 and 100
    
    OM : organic matter % of soil. 

    Returns
    -------
    pwp_wc : volumetric water content at the permanent wilting point.
    
    '''
    pwp_wc = .0260  + perc_clay * .005 + .0158 * OM 
    return pwp_wc



def fieldCapacityRawls2(perc_sand, perc_clay, OM):
    '''
    Predicts the field capacity water content in cm^3/cm^3
    One of two regression equations from:
    W. J. Rawls, ; D. L. Brakensiek, ; K. E. Saxtonn, (1982). Estimation of Soil Water Properties. Transactions of the ASAE, 25(5), 13161320. doi:10.13031/2013.33720 

    Parameters
    ----------
    perc_sand : fraction sand content of soil. between 0 and 100
    
    perc_clay : fraction clay content of soil. between 0 and 100

    
    OM : organic matter % of soil. 

    Returns
    -------
    None.

    '''
    field_capacity_wc = .2576 - perc_sand * .002 + perc_clay * .0036 + OM * .0299
    return field_capacity_wc






def calc_matric_pot(porosity, perc_clay, perc_sand):
    '''Calculate matric potential of a soil's wetting front.
    

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
    matric potential of a soil in pascals

equation from     Rawls and Brakensiek (1985)
Found on SWAT tech docs page 109

'''

    
    matricPot = 10 * exp(
        6.5309 - (7.32561 * porosity) + (0.001583 * perc_clay**2) + (3.809479 * porosity**2) + \
        (.000344 * perc_sand* perc_clay) - (.049837 * perc_sand * porosity) + \
        (.001608 * perc_sand**2 * porosity**2) + (.001602 * perc_clay**2 * porosity **2) -\
        (.0000136 * perc_sand**2 * perc_clay) - .003479 * perc_clay**2 * porosity - \
        (.000799 * perc_sand **2 * porosity)        
        )
    
    
    return matricPot
        
        
def delt_VMC(swc, field_capacity_wc, pwp_wc, porosity):
    '''
    Change in volumetric moisture content across the wetting front

    Parameters
    ----------
    swc : numeric
        current soil water content (mm)
    field_capacity_wc : numeric
        water content at field capacity (mm)
    
    pwp_wc : water content at the permanent wilting point (m)
        
    porosity : soil porosity (mm3/mm3)

    Returns
    -------
    delt_vmc
        Change in volumetric moisture content across the wetting front (mm/mm)

    '''
    
    if field_capacity_wc< pwp_wc or field_capacity_wc<0:
        return np.nan
    if swc >= field_capacity_wc:
        return delt_VMC(field_capacity_wc*.999, field_capacity_wc, pwp_wc, porosity)
    else:
        swc = swc - pwp_wc
        field_capacity_wc = field_capacity_wc- pwp_wc
        delt_vmc = (1-(swc/field_capacity_wc)) * (.95 * porosity)
    return delt_vmc
