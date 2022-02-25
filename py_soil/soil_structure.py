# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 10:14:07 2022

@author: benja
"""


from numpy import exp, log, nan, product


def calc_porosity(bulk_dens, particle_dens):
    '''Calculate porosity cm3/cm3 based on particle density and bulk density.
    
    Parameters
    ----------
    bulk_dens : numeric
        bulk density of the soil Mg/m**3
        
    particle_dens : numeric
         particle density of the soil Mg/m**3
    
    '''
    
    
    porosity =  1 - (bulk_dens/particle_dens)
    if porosity<0:
        return nan
    return porosity



def soil_perm_class(Ksat):
    '''
    Soil Permeability class for wishmeir method of the K factor.
    page 255 SWAT technical docs.

    Parameters
    ----------
    Ksat : numeric
        Saturated Hydraulic Conductivity in mm/hr

    Returns
    -------
    int
       a soil permeability class code. 

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







def struct_stability_index(perc_silt, perc_clay, orgC):
    '''
    Calculate the structural stability index in percents for the soil
    
    From Reynolds, et al 2009
    doi:10.1016/j.geoderma.2009.06.009

    Parameters
    ----------
    
    perc_silt : numeric
        percent silt particles in soil. Range: 0-100
    perc_clay : numeric
        percent clay particles in soil. Range: 0-100
        
    orgC : numeric
        % organic carbon of soil. Range: 0-100
        
    Returns
    -------
    SI : Soil Structural Stability Index
    
    Values >9% indicate stable soil structure
    Values 7-9% are low risk for degradation.
    Values 5-7% are high risk.
    Values below 5% are structurally degratded. 
    '''
    
    SI = 1.724 * orgC / (perc_silt + perc_clay)
    
    return SI




def particle_density(perc_clay, OM):
    '''Predict soil Particle Density based on clay fraction 
    and soil organic matter.
    
    equation from https://doi.org/10.1016/j.geoderma.2016.10.020
    
    Parameters
    ----------
    perc_clay : numeric
        percent clay particles in soil. Range: 0-100
    
    OM : numeric
        % Soil Organic Matter of Soil. Range: 0-100
    
    Returns
    -------
    float 
    particle density in Mg/m**3
    '''
    
    particle_dens = 2.652 + .216 * perc_clay/100 - 2.237 * OM/100
    return particle_dens






def williamsK(perc_clay, perc_silt, perc_sand,  orgC):
    '''Williams equation for calculating the USLE K factor.
    Developed for the EPIC / APEX model.
    
    SWAT technical documentation, page 255

    Parameters
    ----------
    perc_clay : numeric
        percent clay particles in soil. Range: 0-100
    perc_silt : numeric
        percent silt particles in soil. Range: 0-100
    perc_sand : numeric
        percent sand particles in soil. Range: 0-100
    
    orgC: numeric 
        percent organic carbon content of soil. Range: 0-100
    
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
    k_factor_williams = product([csand, cl_si, org_c_f, hisand])
    return k_factor_williams


def wishmK_continuous(perc_silt, perc_vfs, perc_clay, Ksat, soil_str, OM):
    '''Calculate the soil erodibility (K) factor using the wishmeir
    equation, with continuous, rather than values for the soil permeability code.
    
    perc_silt : numeric
        percent silt particles in soil. Range: 0-100
    
    perc_vfs : numeric
        percent very-fine sand particles in soil. Range: 0-100
    
    perc_clay : numeric
        percent clay particles in soil. Range: 0-100
   
    Ksat : Saturated Hydraulic Conductivity (mm/hr)
    
    soil_str : Soil structure code. range (1-4)
    
    OM : numeric
        % Soil Organic Matter of Soil. Range: 0-100
    
    '''
    
    perm_code = max([min([6 - log(Ksat), 6.5]), .9])
    
    return wishmeirK(perc_silt, perc_vfs, perc_clay, perm_code, soil_str, OM)
    


def wishmeirK(perc_silt, perc_vfs, perc_clay, perm_code, soil_str, OM):
    '''
    Wishmeir equation for calculating soil erodibility factor 
    (USLE K).
    For more information, see 
    the SWAT technical documentation page 254.
    

    Parameters
    ----------
    perc_silt : numeric
         percent silt particles in soil. Range: 0-100
    vfs : numeric
        percent very fine sand particles in soil. Range: 0-100
    perc_clay : numeric
        percent clay particles in soil. Range: 0-100
    perm_code : int
        soil permeability class, range:  1-6
    soil_str : int
        soil structure class: range 1-4
        Values correspond to: 
            1: Very Fine Granular
            2: Fine Granular
            3: Medium or Coarse Granular
            4: Blocky, platy, prism-like or massive.
        
    OM : numeric
        percent organic matter in soil. Range: 0-100
        

    Returns
    -------
    float
    USLE Soil erodibility factor K, between 0 and 1
    '''
    
    M = (perc_silt + perc_vfs) * (100-perc_clay)
    k_factor_wishmeir = ( .00021 * M**1.14 * (12-OM) + \
                         3.25*soil_str + 2.5*(perm_code-3))/100
    
    return max(k_factor_wishmeir, 0.01)




