#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 16:37:16 2021

@author: bdube
"""

#notes: exponential norm seems to work best for modelling maxes



import pandas as pd
import geopandas as gpd
import os
from scipy import stats
import numpy as np
import rasterio as rio
from runoff_MUSLE import SCN_eq, slope_adj_CN
import paths



#%%

    



def prep_data(df):
    df = df.rename(columns=
        {col: col[5:] for col in df.columns if 'prcp.' in col})
    df['Year'] = pd.to_datetime(df['date']).apply(lambda x: x.year)
    df['Month'] = pd.to_datetime(df['date']).apply(lambda x: x.month)
    df['prcp'] = df['prcp'].fillna(0)
    return df


def AEP1(df):
    df = df.sort_values('prcp', ascending = False).reset_index(drop =True)
    df['AEP'] = (df.index + 1) / df.shape[0] *365.24
    return df

def AEP2(df, return_periods = [2, 5, 10, 25, 100, 500, 1000]):
    '''
    Estimate Values for annual exceddence probabilities and return periods for precipitation events.

    Parameters
    ----------
    df : dataframe 
        A Dataframe of daily precipitation values for a weather station.

    return_periods : array-like
        return periods (in years) to calculate precipations levels for 
    
    Returns
    -------
    dataframe
        a dataframe of of 

    '''   

    dist = stats.exponnorm(*stats.exponnorm.fit(df.groupby('Year')['prcp'].max()))
    sim = dist.rvs(100000)
    AEPs = 1/np.array(return_periods)
    values = np.quantile(sim, AEPs)
    AEPs = 1- AEPs
    return_periods = 1/AEPs
    return pd.DataFrame({'AEP': AEPs, 'prcp': values, 'return_period': return_periods})



if __name__ == '__main__':
    stations = pd.read_csv(paths.stations)
#%%
    stations = gpd.GeoDataFrame(stations, geometry = gpd.points_from_xy(stations['longitude'], stations['latitude'],
                                             )).set_crs(epsg=4326)
    #%%
    county_bounds = gpd.read_file(paths.county)
    stations.to_crs(county_bounds.crs, inplace = True)
    ax = county_bounds.plot(color = 'b')
    stations.plot(ax=ax, color = 'k')
    AEP_data = []
    
    for station in stations['id'].unique():
        df = prep_data(pd.read_csv(os.path.join('data', 'rainfall', f'{station}.csv')))
        df = df.sort_values('prcp', ascending = False).reset_index(drop =True)
        new_data = AEP2(df)
        new_data['id'] = station
        AEP_data.append(new_data)
        
    AEP_df = pd.concat(AEP_data)
        
        #%%
    
    watersheds = gpd.read_file(
         paths.wshed).to_crs(
            county_bounds.crs)
            
    watersheds = watersheds[watersheds.geometry.intersects(county_bounds.dissolve().geometry.iloc[0])]
    
    #%%
    station_ids= stations['id']
    distances = np.array([watersheds.geometry.centroid.distance(geom) for geom in stations['geometry'] ])
    
    #%%
    closest = [station_ids[ np.argmin(distances[:,i])]
    for i in range(distances.shape[1])]
    
    #%%
    watersheds['nearest_station'] = closest 
    
   
    
    
    