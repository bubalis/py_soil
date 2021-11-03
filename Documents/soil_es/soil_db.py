#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 09:22:52 2021

@author: bdube


Routine to assemble needed data from soil description files downloaded by downloader.py
"""

import os

import pandas as pd
import re
import paths
import os
import subprocess
import pandas as pd
import math
import re

import geopandas as gpd
import requests
import os
import re
import time
from bs4 import BeautifulSoup


#%%
def download_soil_data(names, directory = 'soil_descrips'):
    if not os.path.exists(directory):
        os.makedirs(directory)
    for name in names:
        fp = os.path.join('soil_descrips', f'{name}.txt')
        if not os.path.exists(fp):
            url = f'https://soilseries.sc.egov.usda.gov/OSD_Docs/{name[0]}/{name}.html'
            response = requests.get(url)
       
        
            with open(fp, 'w+') as f:
                f.write(BeautifulSoup(response.text).text)
        time.sleep(2)
    



def matcher(pattern, text, group = 0):
    match = re.search(pattern, text)
    if match:
        return match.group(group)
    else:
        return ''

def assign_struct_code(text):
    if not text:
        return None
    elif 'very fine granular' in text:
        return 1
    elif 'fine granular' in text:
        return 2
    elif  'granular' in text and ('coarse' in text or 'granular' in text):
        return 3
    elif any([b in text for b in ['blocky', 'platy', 'prism', 'massive']]):
        return 4
    else:
        print('Not Found')
        print(text)
        return None


def parse_soil_descrip(f):
    '''Parse a text containg'''
    
    name = f.split('.')[0]
    text = open(f).read()
    soil_tex = matcher(f'TYPICAL PEDON:[\s]+{name} ([\w\s]+),*', text, 1)
    soil_tex = re.split(' [io]n ', soil_tex)[0]
    
    if not soil_tex:
        print(name)

    soil_drain = matcher('(DRAINAGE[\w\s]+:) ([\w\s]+)', text, group =2)
    
    a_horizon = matcher('A[2p\s]*--(.+(?=\n))', text, 1)
    if a_horizon:
        try:
            soil_structure = [p.strip() for p in re.split('[;:]', a_horizon) if 'structure' in p.lower()][0]
        except:
            print(name)
            soil_structure = ''
        
    else:
        B_horizon = matcher('B[hs\s]*--(.+(?=\n))', text, 1)
        if B_horizon:
            try:
                soil_structure = [p.strip() for p in re.split('[;:]', B_horizon) if 'structure' in p.lower()][0]
            except:
                print(name)
                soil_structure = ''
        else:
            soil_structure = ''
    return dict(soil_tex = soil_tex, soil_drain = soil_drain, soil_structure = soil_structure)



def assemble_text_data(_dir = 'soil_descrips'):
    '''Parse all downloaded soil description data into a dataframe.
    Includes: drainage rate, soil texture classification, and soil structure code.
    
    args: _dir- directory location of the textfiles
    
    '''
    
    
    os.chdir(_dir)


    files =  os.listdir()
    

    data = {f.split('.')[0]: parse_soil_descrip(f) for f in files}
    
    
    df = pd.DataFrame(data).T
    
    df['soil_str'] = df['soil_structure'].apply(assign_struct_code)
    
    os.chdir('..')
    df.to_csv(paths.soil_descrips)



def load(dbase_path, tbl, path):
    cmd = ['mdb-export', dbase_path, tbl]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    df = pd.read_csv(proc.stdout)
    return df


def give_first_valid(series):
        for i in series.index:
            item = series[i]
            if isinstance(item, str) and item:
                return item
        return ''

def extract_from_DataBase(db_path, texture_names):
    


    bad_dbs = []
    colnames = {}
    tables = subprocess.check_output(['mdb-tables', db_path]).decode('latin-1')
    tables = tables.strip().split(" ")
    
    
    
    
    #texture_names = [f.split('.')[0] for f in os.listdir('soil_descrips')]
    
    
    
    taxon = load(db_path, 'NCSS_Pedon_Taxonomy')
    taxon = taxon[taxon['samp_name'].isin(texture_names)]
    
    layer = load(db_path, 'NCSS_Layer')
    
    #layer = layer[layer['layer_type'] == 'horizon']
    layer['horizon_key'] = layer[['hzn_desgn_old', 'hzn_desgn', 'hzn_master']].apply(give_first_valid, axis =1 )
    layer = layer[layer['horizon_key'].str.contains('A|B|C|D|E|O')]
    layer['horizon_key'] = layer['horizon_key'].apply(lambda x: re.search('A|B|C|D|E|O', x).group())
    layer= layer[['pedon_key', 'labsampnum', 'horizon_key', 'hzn_top', 'hzn_bot']]
    
    
    
    soil_tex = load(db_path, 'PSDA_and_Rock_Fragments')
    
    taxon = taxon.merge(layer, on = 'pedon_key')
    soil_tex_data = taxon.merge(soil_tex, on = 'labsampnum')
    
    
    columns = [ 'clay_tot_psa', 'silt_tot_psa',
           'sand_tot_psa', 'clay_f', 'co3_cly', 'silt_f_psa', 'silt_c_psa',
           'sand_vf_psa', 'sand_f_psa', 'sand_m_psa', 'sand_c_psa', 'sand_vc_psa',
           'wf_25', 'wf_520', 'wf_2075', 'wf_0175', 'wpg2']
    
    bulk_dens_cols = ['w3cld', 'w15l2', 'db_od', 'db_13b', 'aggstb']
    
    soil_tex_data = soil_tex_data[columns + ['labsampnum', 'samp_name', 'horizon_key', 'hzn_top', 'hzn_bot']]
    

    h2o = load(db_path, 'Bulk_Density_and_Moisture')[bulk_dens_cols + ['labsampnum']]
    soil_tex_data = soil_tex_data.merge(h2o, on = 'labsampnum')
    

    c = load(db_path, 'Carbon_and_Extractions')[['labsampnum', 'oc']]
    soil_tex_data = soil_tex_data.merge(c, on = 'labsampnum')
    
    
    soil_tex_data.to_csv(os.path.join('data', 'soils_data.csv'))
    soil_tex_data['depth'] = soil_tex_data['hzn_bot'] - soil_tex_data['hzn_top']
    
    cols = columns + bulk_dens_cols+ ['oc', 'depth']
    avg_soil_tex_dat = soil_tex_data.groupby(['samp_name', 'horizon_key'])[columns + bulk_dens_cols+ ['oc', 'depth']].mean()
    
    
    avg_soil_tex_dat.to_csv(paths.soil_char_data )
    
    cols = load(db_path, 'NCSS_Data_Dictionary_Data_Tier')
    
    data_dict = cols[(cols['column_name'].isin(columns+['w3cld', 'w15l2', 'oc']+bulk_dens_cols)) & 
                     (cols['data_tier_desc'].isin(['Bulk Density and Moisture', 
                                                   'Carbon and Extractions',
                                                   'PSDA and Rock Fragments']))]
    
    data_dict.to_csv(paths.soil_char_dict)
    
    
if __name__ == '__main__':
    gdf = gpd.read_file(paths.soil_shp)
    names = gdf['MUNAME'].apply(lambda x: re.split('[\s-]', x)[0]).unique()
    del gdf
    download_soil_data(names)
    assemble_text_data()
    extract_from_DataBase(paths.soil_char_db)
    
    

