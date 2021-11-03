#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 11:04:02 2021

@author: bdube
"""

import os
import paths
import geopandas as gpd
import rasterio
from rasterstats import zonal_stats
import numpy as np
import matplotlib.pyplot as plt
from data_utils import CDL_codes_rev, HUC8_codes


#%%


wshed = gpd.read_file(paths.wshed)
county_bounds = gpd.read_file(paths.county)

wshed = wshed[wshed.geometry.intersects(county_bounds.dissolve().geometry.iloc[0])]
wshed = gpd.clip(wshed, county_bounds)
wshed['HUC8'] = wshed['HUC12'].apply(lambda x: x[0:8])

with rasterio.open(paths.lc_rast) as r:
    crs = r.crs
    

def aggregate_ag(dic):
    ag_names = ['hay', 'row_crops', 'fallow', 'tree_crops', 'pasture']
    dic['ag_land_total'] = sum([dic.get(a, 0) for a in ag_names])
    
    return {k: v for k,v in dic.items() if k not in ag_names}


def unpack_zonal_stats(dic):
    out = {}
    for code, nums in CDL_codes_rev.items():
        out[code] = sum([dic.get(n, 0) for n in nums])
    return out
        
    
def perc_dict(dic):
    s = sum(dic.values())
    if s:
        return {k: v/s for k, v in dic.items()}
    else:
        return {}
    

#%%

#%%
res = zonal_stats(wshed.to_crs(crs), raster = paths.lc_rast, categorical = True)
res2 = [unpack_zonal_stats(d) for d in res]
res2 = [aggregate_ag(d) for d in res2]
res2 = [perc_dict(d) for d in res2]



wshed['% forest'] = [r.get('forest', np.nan) for r in res2]
wshed['% ag'] = [r.get('ag_land_total', np.nan) for r in res2]
wshed[wshed['% ag']>.1].plot('% ag', legend = True)
plt.show()
#%%
dic = {code : {'ToHUC': to_huc, 'n_downstream' : None} for code, to_huc in 
       zip(wshed['HUC12'].tolist(),  wshed['ToHUC'].tolist())}



def walk_HUCs(HUC12, dic):
    if dic[HUC12]['n_downstream']:
        return dic[HUC12]['n_downstream']
    elif dic[HUC12]['ToHUC'] not in dic.keys():
        return 0
    else:
        return walk_HUCs(dic[HUC12]['ToHUC'], dic) +1
    
for k in dic:
    dic[k]['n_downstream'] = walk_HUCs(k, dic)

out = {k: v['n_downstream'] for k, v in dic.items()}
wshed['n_downstream'] = wshed['HUC12'].apply(lambda x: out.get(x))



plt.scatter(wshed['n_downstream'], wshed['% ag'])
plt.show()
wshed.plot('n_downstream')

#%%
towns = gpd.read_file(os.path.join('data', 'flood', 'VT_Data_-_Town_Boundaries'))
town_aid = gpd.pd.read_csv(os.path.join('data', 'flood', 'vpr-irene-20130827-federal-map.csv'
                                        ))
towns['TOWNNAME'] = towns['TOWNNAME'].str.title()




flood_relevant = wshed[(wshed['n_downstream']>2) & (wshed['% ag']>=.2)]


wshed['HUC8'] = wshed['HUC12'].apply(lambda x: x[0:8])
#%%
h8 = wshed[['geometry', 'HUC8', '% ag', '% forest']].dissolve(by = 'HUC8', aggfunc = 'mean')

bodies = gpd.read_file(paths.waterbody).to_crs(h8.crs)

#%%
fig, axes = plt.subplots(2, figsize= (6, 10), )

h8.plot('% ag', legend = True, 
        legend_kwds = {'label':  'Proportion land in Ag'}, cmap = 'Reds',
             ax = axes[0])

axes[0].axis('off')



wshed[wshed['HUC8']=='04150402'].plot('% ag', legend = True, cmap = 'Reds',
        legend_kwds = {'label':  'Proportion land in Ag'}, ax = axes[1])
axes[1].axis('off')


for town in ['Middlebury', 'Rutland City', 'Vergennes']:
    point =  towns[towns['TOWNNAME'] == town].centroid
    point.plot(ax =axes[1], 
              color = 'k',
             marker = '*')
    plt.annotate(town.split(' ')[0], (point.x + 5000, point.y -2000))



bodies[bodies['Permanent_'] == '{F9B526FC-7DEA-4DDB-81BB-B01250E74011}'].geometry.boundary.plot(color = 'b', 
                                                                              ax= axes[1], linewidth = .25
                                                                              )
axes[0].set_title('Ag Land Cover by Major Watershed in VT')
axes[1].set_title('Inset of the Otter Creek')
plt.savefig('figures/ag_land.png')
plt.show()
#%%
fig, ax = plt.subplots(figsize = (10, 8))
wshed[wshed['HUC8']==HUC8_codes['Winooski']].plot('% ag', legend = True, cmap = 'Reds',
        legend_kwds = {'label':  'Proportion land in Ag'}, ax = ax)





for town in ['Montpelier', 
             'Waterbury', 
             'Essex',
             ]:
    point =  towns[towns['TOWNNAME'] == town].geometry.centroid
    point.plot(ax =ax, 
              color = 'k',
             marker = '*')
    x = point.x + 2000
    print(x)
    y = point.y -1000
    print(y)
    plt.annotate(town.replace('City', ''), (x, y))
    
bodies[bodies['Permanent_'] =='{EE97D214-7A59-4C44-80E9-BA74CD1CDEA2}'].geometry.boundary.plot(color = 'b', 
                                                                              ax= ax, linewidth = .25)

plt.savefig('figures/winooski.png')
#%%
towns = towns.merge(town_aid, left_on = 'TOWNNAME', right_on = 'town', how = 'left')


towns.plot('fema_total')
plt.show()
#%%



flood_df = gpd.pd.read_csv( '/mnt/c/Users/benja/Documents/soil_es/data/flood/FimaNfipClaims.csv')
flood_df = flood_df[flood_df['state']  =="VT"]

flood_df['CNTY'] = flood_df['countyCode']%1000 


cpi = gpd.pd.read_csv('data/flood/cpi.tsv', sep ='\t')
cpi.rename(columns = {c: c.strip() for c in cpi.columns}, inplace = True)
cpi['series_id'] = cpi['series_id'].str.strip()


cpi = cpi[(cpi['period'] == 'M01') & (cpi['series_id']== 'CWSR0000SA0' )] 

cpi['value'] = cpi['value'] / cpi[cpi['year'] == 2011]['value'].iloc[0]
cpi = cpi[['year', 'value']]

          
flood_df = flood_df.merge(cpi, left_on = 'yearOfLoss', right_on = 'year')
flood_df['amountPaidOnBuildingClaim'] = flood_df['amountPaidOnBuildingClaim'] / flood_df['value']
flood_df['amountPaidOnContentsClaim'] = flood_df['amountPaidOnContentsClaim']/ flood_df['value']

flood_df['CNTY'] = flood_df['countyCode']%1000 
irene_claims = flood_df[(flood_df['yearOfLoss'] == 2011) & (flood_df['CNTY']!=13)]






county_claims = irene_claims.groupby('CNTY')[['amountPaidOnBuildingClaim', 'amountPaidOnContentsClaim']].sum()


county_bounds = county_bounds.merge(county_claims, left_on = 'CNTY', right_index = True)

#county_bounds.plot('amountPaidOnContentsClaim', legend = True)

print('Irene % of 45-year damages:')
print(irene_claims['amountPaidOnBuildingClaim'].sum() / flood_df['amountPaidOnBuildingClaim'].sum())
irene_claims = gpd.GeoDataFrame(irene_claims, 
                                geometry = gpd.points_from_xy(y= irene_claims['latitude'] + np.random.random(irene_claims.shape[0])*.1-.05, 
                                x = irene_claims['longitude'] + np.random.random(irene_claims.shape[0])*.1-.05),
                                ).set_crs(epsg=4326).to_crs(wshed.crs)


#%%

(flood_df.groupby('yearOfLoss')['amountPaidOnBuildingClaim'].sum()/1000000).plot() 
plt.ylabel('Millions of $ (2011)')
plt.xlabel('Year')
plt.title('US Flood Insurance Payouts for Building Damages in VT') 
plt.savefig('figures/flood_by_year.png')
plt.show()

#%%