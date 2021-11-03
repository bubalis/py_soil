#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 08:39:23 2021

@author: bdube
"""

import raster_utils
import os
import geopandas as gpd
import rasterio as rio
import numpy as np

from runoff_MUSLE import SCN_eq, slope_adj_CN
import paths
import shutil
from rasterio.merge import merge
from whitebox import WhiteboxTools
from data_utils import CN_get, CDL_codes_rev, ag_codes
import matplotlib.pyplot as plt

def cdf(x, plot=True, *args, **kwargs):
    x, y = sorted(x), np.arange(len(x)) / len(x)
    return plt.plot(x, y, *args, **kwargs) if plot else (x, y)

# =============================================================================
# def interpolate_nan(array, 
#                     method = 'cubic'):
#     '''
#     Interpolate nan values of a 2d numpy array.
# 
#     Parameters
#     ----------
#     array : 2d numpy array
#         DESCRIPTION.
#     method : string, optional
#         method for 
#         DESCRIPTION. The default is 'cubic'.
# 
#     Returns
#     -------
#     TYPE
#         DESCRIPTION.
#         
#         
#     https://stackoverflow.com/questions/37662180/interpolate-missing-values-2d-python
# 
#     '''
# 
# 
#     x = np.arange(0, array.shape[1])
#     y = np.arange(0, array.shape[0])
#     
#     
#     
#     
#     #mask invalid values
#     array = np.ma.masked_invalid(array)
#     xx, yy = np.meshgrid(x, y)
#     #get only the valid values
#     x1 = xx[~array.mask]
#     y1 = yy[~array.mask]
#     
#     newarr = array[~array.mask]
#     
#     return interpolate.griddata(points = (x1, y1), values = newarr.ravel(),
#                               xi = (xx, yy),
#                                  method= method)
# 
# def interpolate_missing_pixels(
#         image: np.ndarray,
#         mask: np.ndarray,
#         method: str = 'nearest',
#         fill_value: int = 0
# ):
#     """
#     :param image: a 2D image
#     :param mask: a 2D boolean image, True indicates missing values
#     :param method: interpolation method, one of
#         'nearest', 'linear', 'cubic'.
#     :param fill_value: which value to use for filling up data outside the
#         convex hull of known pixel values.
#         Default is 0, Has no effect for 'nearest'.
#     :return: the image with missing values interpolated
#     """
#     from scipy import interpolate
# 
#     h, w = image.shape[:2]
#     xx, yy = np.meshgrid(np.arange(w), np.arange(h))
# 
#     known_x = xx[~mask]
#     known_y = yy[~mask]
#     known_v = image[~mask]
#     missing_x = xx[mask]
#     missing_y = yy[mask]
# 
#     interp_values = interpolate.griddata(
#         (known_x, known_y), known_v, (missing_x, missing_y),
#         method=method, fill_value=fill_value
#     )
# 
#     interp_image = image.copy()
#     interp_image[missing_y, missing_x] = interp_values
# 
#     return interp_image
# 
# =============================================================================

def make_slopes():
    from whitebox import WhiteboxTools
    wshed = gpd.read_file(paths.wshed )
    county_bounds = gpd.read_file(paths.county)

    wshed = wshed[wshed.geometry.intersects(county_bounds.dissolve().geometry.iloc[0])]
    wshed = gpd.clip(wshed, county_bounds)
    wshed['HUC8'] = wshed['HUC12'].apply(lambda x: x[0:8])
    dem_paths = []
    wbt = WhiteboxTools()
    with rio.open(paths.dem) as d:
        dem_crs = d.crs
    for i, row in wshed.dissolve('HUC8').iterrows():
        
        out_path = os.path.join(paths._dir, 'data', '10M_dem_pieces', str(i)+'.tif')
        if not os.path.exists('scratch'):
            os.makedirs('scratch')
        scratch_path = os.path.join(os.getcwd(), 'scratch', 'scratch.shp')
        gpd.GeoDataFrame([row], crs = wshed.crs, geometry = [row.geometry.buffer(20)]).to_file(scratch_path)
        wbt.clip_raster_to_polygon(paths.dem, scratch_path, out_path)
        dem_paths.append(out_path)
    shutil.rmtree('scratch')
    os.chdir(paths._dir)
   
    if not os.path.exists(paths.slope_dir):
        os.makedirs(paths.slope_dir)
    
    for path in dem_paths:
        
        slope_path = os.path.join(paths.slope_dir, os.path.split(path)[1])

        with rio.open(path) as src:
            data = src.read(1) * .3048
            meta = src.meta
            data = data.astype(meta['dtype'])
            data = np.array([data])
        with rio.open(slope_path,  'w+', **meta) as dst:
            dst.write(data)
        del data 
        
        wbt = WhiteboxTools()
        wbt.slope(dem = os.path.join(os.getcwd(), slope_path), 
              output = os.path.join(os.getcwd(), slope_path), units = 'degrees')

        with rio.open(slope_path, 'r+') as r:
            data = r.read(1)
            data = np.where(data>80, 0, data)
            r.write(data, 1)
            
        raster_utils.snap_raster_to_other(slope_path, paths.lc_rast, slope_path, src_crs = dem_crs)
        
    src_to_merge = [rio.open(os.path.join(paths.slope_dir, src)) for src in os.listdir(paths.slope_dir)]
    src_to_merge = [src for src in src_to_merge if src]
        
    rast, out_transform = merge(src_to_merge, method = 'max')
    out_meta = src_to_merge[0].meta.copy()
    out_meta.update({
                     "height": rast.shape[1],
                    "width": rast.shape[2],
                     "transform": out_transform,
                    }
                   )

    for src in src_to_merge:
        src.close()
    
    
    #na_mask = np.where((rast[0] >85) | (rast[0]<=0), True, False)
   
    
    #rast[0] = interpolate_missing_pixels(rast[0], mask = na_mask, 
    #                                     method = 'cubic',
    #                                     fill_value = 0)
    
    #rast[0] = interpolate_nan(rast[0])
   
    with rio.open(paths.slope, 'w+', **out_meta ) as dst:
        
        dst.write(rast.astype(out_meta['dtype']))
    os.chdir(paths._dir)
    
   
#%%
if __name__ == '__main__':
    make_slopes()
    if not os.path.exists(paths.slope):
        make_slopes()
    
    
    
    #%%
    with rio.open(paths.lc_rast) as r:
        lc_data = r.read(1)
        
    with rio.open(paths.irene_rfall) as r:
        rain_data = r.read(1)
        runoff_meta = r.meta
    
    with rio.open(paths.hydrogroup) as r:
        hgrp_data = r.read(1)
        
    with rio.open(paths.slope) as r:
        slope_data = r.read(1) /100
    #to do: adjust CNs to 
    
    rain_data = rain_data * 25.4 #inches to mm
    
    ag_cells = np.isin(lc_data, ag_codes)
    
    
    to_calc = np.where(slope_data != -99.99, True, False)
    slope_data = np.where(slope_data>.8, .08, slope_data)
    slope_data = np.where(slope_data<.033, .033, slope_data)
    
    CN_rast = np.vectorize(CN_get)(lc_data, hgrp_data)
    cdf(CN_rast[ag_cells])
    
    #adjust curve numbers for slope 
    CN_rast[to_calc] = np.vectorize(slope_adj_CN)( 
        CN_rast[to_calc], slope_data[to_calc])
    
    cdf(CN_rast[ag_cells], color = 'r')
    plt.show()
                         
    runoff_data = np.vectorize(SCN_eq)(rain_data, CN_rast)
    #runoff_data = np.where(np.isnan(runoff_data), -9999, runoff_data)
    
    runoff_meta.update({'dtype': np.float32, 'nodata' : -9999})
    #%%
    with rio.open(paths.CN_rast, 'w+', **runoff_meta) as dst:
        dst.write(np.array([CN_rast]).astype(np.float32))
    #%%
    with rio.open(paths.irene_runoff, 'w+', **runoff_meta) as dst:
        dst.write(np.array([runoff_data]).astype(np.float32))


    runoff_data = np.where(runoff_data ==  -9999, np.nan, runoff_data)
        
    
    water_cells = lc_data == CDL_codes_rev['open_water']
    print(runoff_data[ag_cells].sum()/runoff_data[~water_cells].sum())
    
    rd = np.where(rain_data ==0, np.nan, rain_data)
    
    print('Avg Rainfall on Ag land:')
    print(rain_data[ag_cells].mean())
    print('Avg Rainfall on non ag land:')
    print(np.nanmean(rd[~ag_cells]))
    
