#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 15:56:54 2021

@author: bdube
"""
import rasterio as rio
import numpy as np
import raster_utils
import geopandas as gpd
import paths
from rasterio import fill
#%%




pixel_vals = {'>7': (230, 76, 0),
'>6': (255, 170, 0),
'>5': (255, 255, 115),
'>4': (73, 179, 20),
'>3': (113, 204, 67),
'>2': (159, 230, 124),
'<2': (159, 230, 124)}

rfall_vals = {'>7': 8.0,
'>6': 6.5,
'>5': 5.5,
'>4': 4.5,
'>3': 3.5,
'>2': 2.5,
'<2': 1.5}

r_fall_vals = {k: v*25.4 for k, v in rfall_vals.items()}


with rio.open(paths.irene_georeffed) as r:
    band1 = r.read(1)
    band2 = r.read(2)
    band3 = r.read(3)
    out_array = np.zeros(band3.shape) 
    meta = r.meta
    
for key, vals in pixel_vals.items():
    mask = np.where((np.isclose(band1, vals[0]) * np.isclose(band2, vals[1]) \
                     * np.isclose(band3, vals[2])) \
                     ==1)
    out_array[mask] = rfall_vals[key]

gdf = gpd.read_file(paths.county)

meta.update({'count' : 1, 
             'nodata': 0, 
             'dtype' : np.float32})

out_array = fill.fillnodata(out_array, np.where(out_array==0, 0, 1 ))


with rio.open(paths.irene_rfall, 'w+', **meta) as dst:
    dst.write(np.array([out_array]).astype(np.float32))
    




raster_utils.mask_raster(paths.irene_rfall, gdf.geometry)
raster_utils.snap_raster_to_other(paths.irene_rfall, paths.lc_rast, 
                                  out_path = paths.irene_rfall,
                                  )

