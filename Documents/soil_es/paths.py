#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 16:44:28 2021

@author: bdube
"""
import os
from os.path import join as pj


stations = os.path.join('data', 'weather_stations.csv')
_dir = os.getcwd()
county = '/mnt/c/Users/benja/VT_P_index/model/Source_data/VT_Data_-_County_Boundaries-shp/VT_Data_-_County_Boundaries.shp'

irene_georeffed = 'data/rainfall/geo_reffed_irene_rainfall.tif'
irene_rfall = pj(_dir, 'data', 'irene_rfall.tif')
lc_rast = os.path.join(_dir, 'data',
 'flood',
 'land_use_land_cover_NASS_CDL_vt_3951542_01',
 'land_use_land_cover',
 'nass_vt',
 'NASS_VT',
 'cdl_30m_r_vt_2020_utm18.tif')

waterbody = '/mnt/c/Users/benja/VT_P_index/model/Source_data/NHD_H_Vermont_State_Shape/Shape/NHDArea.shp'
wshed = os.path.join(_dir, 'data', 'flood', 'VT_Subwatershed_Boundaries_-_HUC12-shp')
soil_shp = os.path.join(_dir,'data/GeologicSoils_SO/Geologic_SO_poly.shp')
hydrogroup = 'data/hydrogroup.tif'
rainfall_raster = 'data/rainfall/ClimatePrecip_PRECIPA7100/ClimatePrecip_PRECIPA7100.txt'
hgrp_encoding = os.path.join(_dir, 'data', 'hydrogroup_mapping.txt')
irene_runoff = os.path.join(_dir,'intermed_data', 'irene_runoff.tif')

slope = os.path.join(_dir, 'intermed_data', 'slope.tif')
dem = os.path.join(_dir, 'data', '10M_dem.tif')

os.path.join(_dir, 'data', 'flood', 'VT_Subwatershed_Boundaries_-_HUC12-shp')

CN = os.path.join(_dir, 'data/curve_nums.csv')
CN_rast = pj(_dir, 'intermed_data', 'CNs.tif')

soil_descrips = os.path.join(_dir, 'data', 'soil_data_from_descrips.csv')

soil_char_data = os.path.join(_dir, 'data', 'soil_data_summary.csv')

soil_char_dict = os.path.join('data', 'data_dictionary.csv')

soil_char_db = '/mnt/c/Users/benja/Documents/soil_es/NCSS_Soil_Characterization_Database/NCSS_Lab_Data_Mart_09142018/NCSS_Lab_Data_Mart_09142018.mdb'

slope_dir = os.path.join(_dir,'intermed_data/slope_10M')


wwtf_tmdl = os.path.join(_dir, 'data', 'wwtf_tmdl.csv')
wwtf_costs = os.path.join(_dir, 'data', 'VT_wwtf_costs.csv')
