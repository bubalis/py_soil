#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 09:03:56 2021

@author: bdube
"""

import paths
import pandas as pd
import json
import numpy as np


soil_tex_pyramid = { #clay, silt, sand
    'clay': (70, 20, 10),
    'silty clay' : (45, 50, 5),
    'silty clay loam': (35, 55, 10),
    'silt loam': (15, 65, 20), 
    'silt': (7, 85, 8), 
    'loam': (20, 40, 40), 
    'clay loam': (35, 35, 30),
    'sandy clay': (45, 5, 50),
    'sandy clay loam': (25, 10, 60), 
    'sandy loam': (10, 20, 60),
    'loamy sand' : (5, 10, 85),
    'sand': (4, 4, 92),
    'muck': (3,3,3)
    }

soil_data = {'Sand': {'Total Porosity/Saturation': '   0.437',
  'Effective Porosity/Saturation': '0.417',
  'Field Capacity Saturation': '0.091',
  'Wilting Point Saturation': '0.033',
  'Residual Saturation': '0.02',
  'Bubbling Pressure (mean)': '7.26',
  'Porse size distribution (mean)': '0.694',
  'Saturated Hydrolic Conductivity': '23.56',
  'Wetting Front Suction Head': '4.95'},
 'Loamy sand': {'Total Porosity/Saturation': '0.437',
  'Effective Porosity/Saturation': '0.401',
  'Field Capacity Saturation': '0.125',
  'Wilting Point Saturation': '0.055',
  'Residual Saturation': '0.035',
  'Bubbling Pressure (mean)': '8.69',
  'Porse size distribution (mean)': '0.553',
  'Saturated Hydrolic Conductivity': '5.98',
  'Wetting Front Suction Head': '6.13'},
 'Sandy loam': {'Total Porosity/Saturation': '0.453',
  'Effective Porosity/Saturation': '0.412',
  'Field Capacity Saturation': '0.207',
  'Wilting Point Saturation': '0.095',
  'Residual Saturation': '0.041',
  'Bubbling Pressure (mean)': '14.66',
  'Porse size distribution (mean)': '0.378',
  'Saturated Hydrolic Conductivity': '2.18',
  'Wetting Front Suction Head': '11.01'},
 'Loam': {'Total Porosity/Saturation': '0.463',
  'Effective Porosity/Saturation': '0.434',
  'Field Capacity Saturation': '0.27',
  'Wilting Point Saturation': '0.117',
  'Residual Saturation': '0.027',
  'Bubbling Pressure (mean)': '11.15',
  'Porse size distribution (mean)': '0.252',
  'Saturated Hydrolic Conductivity': '1.32',
  'Wetting Front Suction Head': '8.89'},
 'Silt loam': {'Total Porosity/Saturation': '0.501',
  'Effective Porosity/Saturation': '0.486',
  'Field Capacity Saturation': '0.33',
  'Wilting Point Saturation': '0.133',
  'Residual Saturation': '0.015',
  'Bubbling Pressure (mean)': '20.79',
  'Porse size distribution (mean)': '0.234',
  'Saturated Hydrolic Conductivity': '0.68',
  'Wetting Front Suction Head': '16.68'},
 'Sandy clay loam': {'Total Porosity/Saturation': '0.398',
  'Effective Porosity/Saturation': '0.330',
  'Field Capacity Saturation': '0.255',
  'Wilting Point Saturation': '0.148',
  'Residual Saturation': '0.068',
  'Bubbling Pressure (mean)': '28.08',
  'Porse size distribution (mean)': '0.319',
  'Saturated Hydrolic Conductivity': '0.30',
  'Wetting Front Suction Head': '21.85'},
 'Clay loam': {'Total Porosity/Saturation': '0.464',
  'Effective Porosity/Saturation': '0.390',
  'Field Capacity Saturation': '0.318',
  'Wilting Point Saturation': '0.197',
  'Residual Saturation': '0.075',
  'Bubbling Pressure (mean)': '25.89',
  'Porse size distribution (mean)': '0.242',
  'Saturated Hydrolic Conductivity': '0.20',
  'Wetting Front Suction Head': '20.88'},
 'Silty clay loam': {'Total Porosity/Saturation': '0.471',
  'Effective Porosity/Saturation': '0.432',
  'Field Capacity Saturation': '0.366',
  'Wilting Point Saturation': '0.208',
  'Residual Saturation': '0.040',
  'Bubbling Pressure (mean)': '32.56',
  'Porse size distribution (mean)': '0.177',
  'Saturated Hydrolic Conductivity': '0.20',
  'Wetting Front Suction Head': '27.30'},
 'Sandy clay': {'Total Porosity/Saturation': '0.430',
  'Effective Porosity/Saturation': '0.321',
  'Field Capacity Saturation': '0.339',
  'Wilting Point Saturation': '0.239',
  'Residual Saturation': '0.109',
  'Bubbling Pressure (mean)': '29.17',
  'Porse size distribution (mean)': '0.223',
  'Saturated Hydrolic Conductivity': '0.12',
  'Wetting Front Suction Head': '23.90'},
 'Silty clay': {'Total Porosity/Saturation': '0.479',
  'Effective Porosity/Saturation': '0.423',
  'Field Capacity Saturation': '0.387',
  'Wilting Point Saturation': '0.250',
  'Residual Saturation': '0.056',
  'Bubbling Pressure (mean)': '34.19',
  'Porse size distribution (mean)': '0.150',
  'Saturated Hydrolic Conductivity': '0.10',
  'Wetting Front Suction Head': '29.22'},
 'Clay': {'Total Porosity/Saturation': '0.475',
  'Effective Porosity/Saturation': '0.385',
  'Field Capacity Saturation': '0.396',
  'Wilting Point Saturation': '0.272',
  'Residual Saturation': '0.090',
  'Bubbling Pressure (mean)': '37.30',
  'Porse size distribution (mean)': '0.165',
  'Saturated Hydrolic Conductivity': '0.06',
  'Wetting Front Suction Head': '31.63'}}


HUC8_codes = {'Winooski' : '04150403',
              'Otter' : '04150402',
              'White' : '01080105',
              'Black/Otaqueechee': '01080106',
              'Lamoille': '04150405',
              "Missisquoi" : "04150407", 
              'Memphramagog and Canada' : "04150500", 
              'Batten Kill': "02020003",
              'West River': "01080107",
              }


CDL_codes_rev = dict(
    hay = [37, 58, 59, 61], 
    temp_hay = [36],
                 row_crops = [1,4,5,11,12, 43, 44, 47, 54, 216, 221, 222, 228, 229, ],
                 tree_crops = [68, 69, 70, 71, 242],
                 developed_open = [121],
                 developed_low = [122],
                 developed_medium = [123],
                 developed_high = [ 124],
                 fallow = [131],
                 forest = [63, 141, 142, 143, 152],
                 pasture = [176],
                 wetland = [190, 195], 
                 shrubs = [64],
                 small_grains = [6, 21,22, 23, 24, 25,26, 27, 28, 29, 30, 39,
                                 
                                 ],
                 open_water = [111]
                 )
lu_codes ={}
for k,v in CDL_codes_rev.items():
    for i in v:
        lu_codes[i] = k


CN_df = pd.read_csv(paths.CN)
CN_df = CN_df.groupby('crop cover')[['A', 'B', 'C', 'D']].mean()

CN_df.index = [i.lower().replace(' ',  '_') for i in CN_df.index]
CN_df['water'] = CN_df['D']
CN_df['not rated'] = CN_df['C']

drained_lcs =['row_crops', 'temp_hay', 'small_grains']
for st in ['A/D', 'B/D', 'C/D']:
    CN_df[st] = CN_df[st.split('/')[1]]
    CN_df.loc[drained_lcs, st] = CN_df.loc[drained_lcs, st.split('/')[0]]

    
with open(paths.hgrp_encoding) as f:
    hgrp_encoding = json.loads(f.read())

CN_dict = {}
for name, row in CN_df.iterrows():
    
    CN_dict[name] =  {hgrp_encoding[c]: row[c] for c in CN_df.columns} 

CN_dict = {code: CN_dict[name] for code, name in lu_codes.items()}



ag_codes = []
for s in ['pasture', 'small_grains', 'tree_crops', 'hay', 'temp_hay', 'row_crops']:
    ag_codes += CDL_codes_rev[s]


def CN_get(lc_name, hydro_group):
    return CN_dict.get(lc_name, {}).get(hydro_group, np.nan)