#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 14:46:51 2021

@author: bdube
"""

import rasterio as rio
from rasterio.mask import mask
import os
import ogr
import gdal
import geopandas as gpd
import json
from rasterio.warp import calculate_default_transform, reproject, Resampling
import random
import string
from shapely.geometry import box
import numpy as np
import paths
from whitebox import WhiteboxTools
import shutil


def resize_raster(filepath, factor, out_path):
    '''Decrease the size of a raster by a scaling factor.'''
    with rio.open(filepath) as dataset:

    # resample data to target shape
        data = dataset.read(
            out_shape=(
                dataset.count,
                int(dataset.height / factor),
                int(dataset.width / factor)
            ),
            resampling = Resampling.bilinear
        )
    
        # scale image transform
        out_transform = dataset.transform * dataset.transform.scale(
            (dataset.width / data.shape[-1]),
            (dataset.height / data.shape[-2])
        )
        out_meta = dataset.meta.copy()
        out_meta.update({"transform": out_transform,
                     'width': data.shape[-1],
                     'height': data.shape[-2],
                     "driver": "GTiff",}
                   )
    with rio.open(out_path, 'w+', **out_meta) as dest:
        dest.write(data)

def rand_fp(length, tail = '.tif', head = ''):
    '''Generate a random filepath of length (length).
    tail: the end of the filepath, i.e. file type
    head: directory or other.
    '''
    return head + ''.join(
        random.choice(string.ascii_uppercase + string.digits) 
        for _ in range(length)) + tail


def Feature_to_Raster(input_shp, output_tiff, snap_raster,
                      field_name = False, NoData_value=-9999, data_type = gdal.GDT_Float32):
    """
    Converts a shapefile into a raster
    From:
    https://www.programcreek.com/python/example/101827/gdal.RasterizeLayer
    #6
    
    args: 
    input_shp: filepath for the shapefile
    output_tiff: filepath for output
    snap_raster: filepath to snap_raster
    field_name: optional. If passed, determines which data field in the shapefile
                to rasterize
    NoData_value: Value for missing data. Defaults to -9999
    data_type: type of data for the output file.
    """

    # Input
    try:
        inp_driver = ogr.GetDriverByName('ESRI Shapefile')
        
        inp_source = inp_driver.Open(input_shp, 0)
        inp_lyr = inp_source.GetLayer()
        inp_srs = inp_lyr.GetSpatialRef()
    except:
        globals().update(locals())
        raise

    # Extent
    with rio.open(snap_raster) as example:
        cellsize = example.transform[0]
        x_min = example.transform[2]
        y_max = example.transform[5]
        y_ncells, x_ncells = example.shape
        

    # Output
    out_driver = gdal.GetDriverByName('GTiff')
    if os.path.exists(output_tiff):
        out_driver.Delete(output_tiff)
    out_source = out_driver.Create(output_tiff, x_ncells, y_ncells,
                                   1, data_type)

    out_source.SetGeoTransform((x_min, cellsize, 0, y_max, 0, -cellsize))
    out_source.SetProjection(inp_srs.ExportToWkt())
    out_lyr = out_source.GetRasterBand(1)
    out_lyr.SetNoDataValue(NoData_value)

    # Rasterize
    if field_name:
        gdal.RasterizeLayer(out_source, [1], inp_lyr,
                            options=["ATTRIBUTE={0}".format(field_name)])
    else:
        gdal.RasterizeLayer(out_source, [1], inp_lyr, burn_values=[1])

    # Save and/or close the data sources
    inp_source = None
    out_source = None

    return output_tiff 





#%%
def mask_raster(raster_path, shapes, out_path = None, crs = None):
    '''Mask a raster with the geometries given in shapes.
    Save to out_path. If out_path is not specified, save to original path.'''
    
    if not crs:
        crs = shapes.crs
    
    with rio.open(raster_path) as src:
        rast_crs = src.crs
    
    if rast_crs != crs:
        print('Raster CRS does not match mask, reprojecting')
        reproj_raster(raster_path, raster_path, crs)
   
    with rio.open(raster_path) as src:
        out_image, out_transform = mask(src, shapes, crop = True)
        out_meta = src.meta
    
    out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_transform})
    
    
    if not out_path:
        out_path = raster_path
        
    with rio.open(out_path, "w+", **out_meta) as dest:
        dest.write(out_image)

def reproj_raster(file_path, out_path, crs):
    with rio.open(file_path) as src:
        transform, width, height = calculate_default_transform(
            src.crs, crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': crs,
            'transform': transform,
            'width': width,
            'height': height
        })
    
        with rio.open(out_path, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source = rio.band(src, i),
                    destination = rio.band(dst, i),
                    src_transform = src.transform,
                    src_crs = src.crs,
                    dst_transform = transform,
                    dst_crs = crs,
                    resampling = Resampling.nearest)


def encode_as_ints(gdf, column):
    '''Create encoding for categorical values to use in encoding rasters.
    args: gdf- a geodataframe
    column- str: a column name of that geodataframe.
    
    returns: the gdf with the categorical value reassigned as a number.
    
        mapping- a dictionary of the form {category_name: number}
        '''
    
    mapping  =  {val:i for i, val in enumerate(gdf[column].unique())}
    gdf[column]  =  gdf[column].apply(lambda x: mapping.get(x))
    return gdf, mapping



def categorical_feature_to_raster(input_shp, output_tiff, snap_raster,
                      column, NoData_value=-9999, 
                      data_type = gdal.GDT_Float32):
    
    output_tiff = os.path.join(os.getcwd(), output_tiff)
    snap_raster = os.path.join(os.getcwd(), snap_raster)
    with rio.open(snap_raster) as s:
        crs = s.crs
    
    gdf = gpd.read_file(input_shp)
    gdf, mapping = encode_as_ints(gdf, column)
    gdf = gdf[['geometry', column]].to_crs(crs)
    scratch =  rand_fp(5, tail = '')
    
    scratch_shp = os.path.join(os.getcwd(), 'data', scratch, scratch + '.shp')
    
    map_path = output_tiff.replace('.tif', '_mapping.txt')
    
    with open(map_path, 'w+') as f:
        f.write(json.dumps(mapping))
    os.makedirs(os.path.split(scratch_shp)[0])
    
    #%%
    gdf.dissolve(column)
    gdf.to_file(scratch_shp)
    
    wbt = WhiteboxTools()
    start_dir = os.getcwd()
    wbt.vector_polygons_to_raster(scratch_shp, output_tiff, field = column,
                                 base = snap_raster)
    #Feature_to_Raster(scratch_shp, output_tiff, snap_raster, 
    #                  column, NoData_value, data_type)
    os.chdir(start_dir)
    
    shutil.rmtree(os.path.split(scratch_shp)[0])
    
    
def snap_raster_to_other(raster_path, ex_path, out_path, src_crs= None):
    '''Resample a raster to match the size and extent of another raster.
    args : 
        raster_path - raster to alter.
        ex_path - raster to use as template
        out_path - path of output raster.
    '''

    #collect template data
    with rio.open(ex_path) as ex:
    
        dst_res = (ex.transform[0], ex.transform[4] * -1)
        dst_crs = ex.crs
        height = ex.height
        width = ex.width
        dst_transform = ex.transform
        bounds = ex.bounds
        crs = ex.crs

    if src_crs:
        with rio.open(raster_path, 'r') as src:
            data = src.read()
            meta = src.meta
       
    
        meta.update({'crs': src_crs})
        with rio.open(raster_path, 'w+', **meta) as dst:
            dst.write(data)
    

    src = rio.open(raster_path)
    
    #if extents dont match, mask the raster to bounding box of the example.
    if src.bounds!= bounds:
        src.close()
        mask_raster(raster_path, [box(*bounds)], out_path, crs)
        src = rio.open(out_path)
    
    # make the out array
    dst_data = np.zeros((1, height, width))
    
    try: 
        rio.warp.reproject(rio.band(src, 1), 
                                destination = dst_data,
                            dst_resolution = dst_res, 
                            dst_crs = dst_crs,
                            src_crs = src.crs,
                            dst_transform = dst_transform, )
        
        out_meta = src.meta.copy()
    finally:
        src.close()
    
    assert width == dst_data.shape[-1]
    assert height == dst_data.shape[-2]
    
    out_meta.update({'transform': dst_transform,
                     'width': dst_data.shape[-1],
                     'height': dst_data.shape[-2],
                     })
    
    with rio.open(out_path, 'w+', **out_meta) as dst:
        dst.write(dst_data.astype(out_meta['dtype']))
        
    assert rio.open(out_path).read(1).shape == rio.open(ex_path).read(1).shape
#%%

if __name__ == '__main__':
    #if not os.path.exists(paths.hydrogroup):
        categorical_feature_to_raster(paths.soil_shp, 'data/hydrogroup.tif', paths.lc_rast, 'HYDROGROUP', 
                  NoData_value = -1)
    
    
#%%

