#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 17:40:31 2021

@author: bdube
"""

#from libc.stdio cimport 
from libc.math cimport log, abs, fabs, round
from cpython cimport array
import cython

import numpy as np
cimport numpy as np

cdef double dubsum(double[:] arr):
    total =  0
    for item in arr:
        total += item
    return total

cdef (double, double) timestep1(double Ke, double matricPot, double delt_vmc, double cum_inf, 
                double R_fall, double time_step = 1, double prev_inf = 0):
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
    
    
    cdef double new_inf_rate, term_1, term_2, inf_rate, runoff
    new_inf_rate = 0
    
    term_1 = Ke * time_step
    term_2 = matricPot * delt_vmc
    inf_rate = term_1
    while True:
        new_inf_rate = term_1 + term_2 * log(((inf_rate + cum_inf) + term_2) / (cum_inf + term_2 ))
        if fabs(inf_rate - new_inf_rate) <.001:
            
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
cdef (double, double) timestep2(double Ke, double matricPot, double delt_vmc, double cum_inf, 
                double R_fall, double time_step = 1, double prev_inf = 0):
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
    
    cdef double new_inf_rate, term_1, term_2, inf_rate, runoff
    
    inf_rate = 3
    new_inf_rate = 0
    
    term_1 = Ke * time_step
    term_2 = matricPot * delt_vmc
    while True:
        
        new_inf_rate = term_1 + term_2* log(1 + ((cum_inf+inf_rate)/term_2))
       
        
        if fabs(inf_rate- new_inf_rate) < .001:
            
            break
        else:
            inf_rate = new_inf_rate
            
    if new_inf_rate >R_fall:
        new_inf_rate = R_fall
        runoff = 0
    else:
         runoff = R_fall - new_inf_rate
    return new_inf_rate, runoff

#%%
cdef (double, double) timestep3(double Ke, double matricPot, double delt_vmc, double cum_inf, 
                double R_fall, double time_step = 1, double prev_inf = 0):
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
    
    
    cdef double new_inf_rate, term_1, term_2, inf_rate, runoff
    
    term_1 = Ke * time_step
    term_2 = matricPot * delt_vmc
    
    new_inf_rate = term_1 * ((term_2+1) /(cum_inf +1e-6))
    
    if new_inf_rate >R_fall:
        new_inf_rate = R_fall
        runoff = 0
    else:
         runoff = R_fall - new_inf_rate
    return new_inf_rate, runoff


@cython.cdivision(True)
cdef (double, double) timestep4(double Ke, double matricPot, double delt_vmc, double cum_inf, 
                double R_fall, double time_step = 1, double prev_inf = 0):
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
    
    
    cdef double new_inf_rate, term_1, term_2, inf_rate, runoff, tst, f1
    cdef int counter
    
    term_1 = Ke * time_step 
    term_2 = matricPot * delt_vmc
    
    if prev_inf:
        tst = prev_inf
    else:
        tst = term_1 
    
    counter = 0
    while True:
        
        f1 = cum_inf +  term_1  +   term_2 * log((tst + term_2)/(cum_inf + term_2))
        
         
        if fabs(f1 - tst) < 0.001:
            break
        else:
            tst = f1
        counter = counter+ 1
        #if counter>10000:
        #    raise ValueError()
    new_inf_rate = f1 - cum_inf
    if new_inf_rate >R_fall:
        new_inf_rate = R_fall
        runoff = 0
    else:
         runoff = R_fall - new_inf_rate
    return new_inf_rate, runoff

cdef int wfLayer( double swcs[3], double fcs[3]):
    cdef int i
    for i in range(0,4):
        if swcs[i]< fcs[i]:
            return i
    return i



cdef (double, double) timestep5(
    double Ke, double matricPot, double delt_vmc, double cum_inf, 
                    double R_fall, double time_step = 1, double prev_inf = 0):
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
    cdef double term_1, term_2, inf_rate, runoff, tst, f1
    #cdef int counter
    
    
    term_1 = Ke * time_step
    term_2 = matricPot * delt_vmc
    if prev_inf:
        tst = prev_inf
    else:
        tst = term_1 
    while True:
        f1 = term_1 * (1 + term_2)/(cum_inf + tst + .00001)
                      
        
        
        
        if (abs(f1 - tst) <= 0.001):
            break
        else:
            tst = f1
            
    

    if f1 >R_fall:
        f1 = R_fall
        runoff = 0
    else:
         runoff = R_fall - f1
    return f1, runoff



@cython.cdivision(True)
cdef double multi_layer_Keff(double cum_inf, int wf_layer, 
                            double delt_vmcs[3],
                            double depths[3],
                            double  Kes[3]):
    
    cdef double term_1, num, denom
    if wf_layer == 0 :
        return Kes[0]
    elif wf_layer == 1:
        term1 = (cum_inf - (depths[0] * delt_vmcs[0]))/delt_vmcs[1]
        
        num = depths[0] + term1
        denom = (depths[0]/Kes[0]) + (term1/Kes[1])
        return num/denom
    else:
        term1 = (cum_inf - (depths[0] * delt_vmcs[0]) + (depths[1] * delt_vmcs[1]))/delt_vmcs[2]
        num = depths[0] + depths[1] + term1
        denom = depths[0]/ Kes[0] + depths[1]/Kes[1] + term1 / Kes[2]
        return num/denom


@cython.cdivision(True)
cdef (double, double) GA_storm_multilayer(
    double[:] rainfalls,
    int n_steps, 
    double layer_water_cont[4],
    double layer_fc[4],
    double layer_pwp[4],
    double layer_porosity[4],
    double layer_matric_p[4],
    double layer_KE[4],
    double layer_depth[4], 
    double time_step):
    '''
    Simullate a storm with a multilayer green-ampt model.

    Parameters
    ----------
    total_rfall : numeric
        rainfall quantity (mm)
    hrs : int
        number of hours.
    layer1 : dict or pandas series containing information about the top soil series
        contains .
    layer2 : dict or pandas series containing information about the top soil series
        DESCRIPTION.
    layer3 : dict or pandas series containing information about the top soil series
        DESCRIPTION.
    time_step : TYPE, optional
        DESCRIPTION. The default is 1/60.
    method : TYPE, optional
        DESCRIPTION. The default is timestep4.

    Returns
    -------
    TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.

    '''
    #layers = layers.set_index('horizon_key').sort_index()
    
    #if np.any(layers.loc['A', ['Ke', 'matric_pot', 'swc', 'field_capacity', 'porosity']].isna()):
    #    return np.nan, np.nan
   
    #if layers.loc['A', 'field_capacity']<layers.loc['A', 'pwp_capacity']:
    #    return np.nan, np.nan
    
    #layer_names = layers.index.tolist()
    
    cdef int c, wf_layer, old_wf_layer, new_wf_layer
    
     
    
    
    cdef double cum_inf 
    cdef double inf 
    
    cdef double[3] delt_vmcs 
    cdef double[:] runoffs 
    #cdef double[n_steps] infs
    


    cum_inf = 0
    inf = 0
    runoffs = np.zeros(n_steps)
    
    

    
    
    #wf_layers = []
    old_wf_layer = 0
    for i in range(0,4):
        delt_vmcs[i] = calc_deltVMC(layer_water_cont[i], 
                                    layer_fc[i], layer_pwp[i], layer_porosity[i])
    
    for c in range(len(rainfalls)):
        rfall = rainfalls[c]
       
        
        wf_layer = wfLayer(layer_water_cont, layer_fc)
        
        
        #if wf_layer!= old_wf_layer:
        #   delt_vmcs[old_wf_layer] = calc_deltVMC(layer_water_cont[old_wf_layer], layer_fc[old_wf_layer], 
        #                                          layer_pwp[old_wf_layer], layer_porosity[old_wf_layer])
            
        #method for calculating multilayer GA- equation with multiple layers
        # from https://www.gsshawiki.com/Infiltration:Multi-layer_Green_and_Ampt
        
        
        if wf_layer == 3: 
            if layer_water_cont[wf_layer] >= layer_fc[wf_layer]:
                Ke = .01        
            else:
                Ke = layer_KE[3]
            MatricPot = layer_matric_p[3]
            delt_vmc = .001
        
        else:
            
            matricPot = layer_matric_p[wf_layer]
            delt_vmc = delt_vmcs[wf_layer]
            Ke =  multi_layer_Keff(cum_inf, wf_layer, 
                                     delt_vmcs,
                                     layer_depth,
                                    layer_KE)
        
        inf, runoff = timestep4(Ke, matricPot , delt_vmc, cum_inf, 
                                     rfall, time_step, inf)
        
        
        
        cum_inf = cum_inf + inf
        
        layer_water_cont[wf_layer] = layer_water_cont[wf_layer] + inf
        
        #    layer_fc[wf_layer] = calc_deltVMC(layer_water_cont[wf_layer], layer_fc[wf_layer], 
        #                                           layer_pwp[wf_layer], 
        #                                           layer_porosity[wf_layer])
        
        
        old_wf_layer = wf_layer
        runoffs[c] = runoff
        
    runoff = dubsum(rainfalls) - cum_inf
    if runoff<0:
        runoff = 0
    
    
    max_runoff = np.max(runoffs) / (4*time_step)

    
    return runoff,  max_runoff
    
cdef double calc_deltVMC(double swc, double field_capacity, double pwp_capacity, double porosity):
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
    cdef double delt_vmc
    if field_capacity< pwp_capacity or field_capacity<0:
        return 0
    if swc >= field_capacity:
        return calc_deltVMC(field_capacity*.999, field_capacity, pwp_capacity, porosity)
    else:
        swc = swc - pwp_capacity
        field_capacity = field_capacity- pwp_capacity
        delt_vmc = (1-(swc/field_capacity)) * (.95 * porosity)
    return delt_vmc






    
def ts1(Ke, matricPot, delt_vmc, cum_inf, R_fall, time_step, prev_inf):
    return timestep1( Ke, matricPot, delt_vmc, cum_inf, R_fall, time_step, prev_inf )

def ts2(Ke, matricPot, delt_vmc, cum_inf, R_fall, time_step, prev_inf):
    return timestep2( Ke, matricPot, delt_vmc, cum_inf, R_fall, time_step, prev_inf )

def ts3(Ke, matricPot, delt_vmc, cum_inf, R_fall, time_step, prev_inf):
    return timestep3( Ke, matricPot, delt_vmc, cum_inf, R_fall, time_step, prev_inf )

def ts4(Ke, matricPot, delt_vmc, cum_inf, R_fall, time_step, prev_inf):
    return timestep4( Ke, matricPot, delt_vmc, cum_inf, R_fall, time_step, prev_inf )

def ts5(Ke, matricPot, delt_vmc, cum_inf, R_fall, time_step, prev_inf):
    return timestep5( Ke, matricPot, delt_vmc, cum_inf, R_fall, time_step, prev_inf )


 
    

def GA_storm_multilayer_df(layers, total_rfall,  hrs, time_step, rainfalls = []):
    '''API for wrapping the multi-layer green-ampt simulation.'''
    
    cdef array.array swc = array.array('d', layers['swc'].to_list())
    cdef array.array fc = array.array('d',  layers['field_capacity'].to_list())
    cdef array.array pwp = array.array('d', layers['pwp_capacity'].to_list())
    cdef array.array porosity = array.array('d', layers['porosity'].to_list())
    cdef array.array matric_p = array.array('d', layers['matric_pot'].to_list())
    cdef array.array Ke = array.array('d', layers['Ke'].to_list())
    cdef array.array depth = array.array('d', layers['depth'].to_list())
    
    
    if hrs and not rainfalls:
        n_steps =int(round(hrs /time_step))
        tstep_rfall = total_rfall/n_steps
        rainfalls = np.repeat(tstep_rfall, n_steps)
    else:
        n_steps = len(rainfalls)
        
    cdef array.array rainfall_array = array.array('d', rainfalls)
    
    try:
        return GA_storm_multilayer(rainfalls = rainfall_array, 
                               n_steps = n_steps, 
                                  layer_water_cont = swc.data.as_doubles, 
                                  layer_fc = fc.data.as_doubles,  
                                  layer_pwp = pwp.data.as_doubles, 
                                  layer_porosity = porosity.data.as_doubles,
                                  layer_matric_p = matric_p.data.as_doubles, 
                                  layer_KE = Ke.data.as_doubles, 
                                  layer_depth = depth.data.as_doubles,
                                  time_step = time_step)
    except:
        print(locals())
        raise



def test_ml_Keff(layers, double cum_inf, int wf_layer,):
    cdef array.array delt_vmcs = array.array('d', layers['delt_vmc'].to_list())
    cdef array.array depths  = array.array('d',  layers['depth'].to_list())
    cdef array.array Kes  = array.array('d', layers['Ke'].to_list())


    return multi_layer_Keff(cum_inf, wf_layer, delt_vmcs.data.as_doubles, depths.data.as_doubles, Kes.data.as_doubles)

def test_timestep4(Ke, matricPot , delt_vmc, cum_inf, 
                             timestep_rfall, time_step, inf):
    return timestep4(Ke, matricPot , delt_vmc, cum_inf, 
                                 timestep_rfall, time_step, inf)