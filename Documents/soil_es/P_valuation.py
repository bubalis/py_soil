#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 09:07:21 2021

@author: bdube
"""

import paths
import pandas as pd



def fix_costs_df():
    costs = pd.read_csv(paths.wwtf_costs)
    costs['Percent Hydraulic Capacity Used']= costs['Percent Hydraulic Capacity Used'].apply(lambda x: x.replace('%', '')).astype(int) / 100
    
    costs['Percent P'] = costs['Percent P'].apply(lambda x: x.replace('%', '')).astype(int) / 100
    
    for c in costs.columns[6:]:
        if c== 'Combined Scenario':
            continue
        costs[c] = costs[c].apply(lambda x: str(x).strip().split(' ')[0].replace('*', '')).astype(float)
    costs['Facility'] = costs['Facility'].str.replace('\xa0', ' ')
    costs['Facility'] = costs['Facility'].str.replace('IBM', 'Global Foundries'
                                                      ).str.replace('Wyeth \(PBM Nutritionals\)', 'Perrigo  Nutritionals'
                                                                    ).str.replace('\n', '')
    costs.to_csv(paths.wwtf_costs, index = False)


wwtf_alloc = pd.read_csv(paths.wwtf_tmdl)
costs = pd.read_csv(paths.wwtf_costs)


