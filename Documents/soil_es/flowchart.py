#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 11:58:08 2021

@author: bdube
"""

import pygraphviz as pgv

#%%
freshwater_outcomes = ['ALS Incidence', 'Tourism', 'Home Prices']
saltwater_outcomes = ['Fisheries / Benefits to CT/NYS', ]

nodes = {
         'Soil Health Improvements': {'Less Erosion': {'attr': {'color' : 'blue'}}, 'Less Runoff': None,
                      'Carbon Sequestration': None,
                      'Improved N Cycling': None
                      },
         'Improved N Cycling': {'Climate Mitigation': None,
                                'Lower N Loading': None,
                                'Lower Groundwater N' : None},
         
         'Less Runoff' : {'Lower P Loading': None, 'Lower N Loading': None,
                          'Flood Mitigation': None
                          },
         
         'Carbon Sequestration': {'Climate Mitigation': None},
         'Climate Mitigation': {'Social Cost of Carbon': None},
         'Less Erosion': {'Lower P Loading': None},
         
         'Lower P Loading': {'Improved Water Quality (Fresh)': None},
         'Lower N Loading': {'Improved Water Quality (Ocean)': None},
         'Lower Groundwater N': {'Mitigation/Treatment Costs': None},
         'Improved Water Quality (Fresh)': freshwater_outcomes,
         'Improved Water Quality (Ocean)': saltwater_outcomes,
         'Flood Mitigation': {'Avoided Damages'},
         
         }

g = pgv.AGraph(directed = True)



o = freshwater_outcomes + saltwater_outcomes + ['Avoided Damages' , "Social Cost of Carbon", 'Mitigation/Treatment Costs']
g.add_nodes_from(nodes.keys())

#soil_health =  g.add_subgraph(['Soil Health Improvements'], rank = "same")
#processes = g.add_subgraph(['Improved N Cycling', 'Less Runoff', 'Less Erosion', "Carbon Sequestration"], rank ='same')
#sub_outcomes = g.add_subgraph(['Improved Water Quality (Fresh)', 'Improved Water Quality (Ocean)'], rank = "same")


for n, v in nodes.items():
    for i in v:
        g.add_edge(n, i)

outcomes = g.add_subgraph(o, rank = 'same')
g.draw('thought_model.png', prog = 'dot') 




#%% Abstract_graph

abs_nodes = {'Practice': {'Soil Health Metric': None, 'Ecological Process': None},
             'Soil Health Metric': {'Ecological Process': None},
             'Ecological Process': {'Ecosystem Service': None},
             'Ecosystem Service': {'Valuation': None}
             }

g2 = pgv.AGraph(abs_nodes, directed = True)



g2.layout('dot')
g2.draw('figures/abstract.png')


#%%
ex_nodes = {'Cover Crop': {'Soil Organic Matter': None, 'Rainfall Infiltration': None},
             'Soil Organic Matter': {'Rainfall Infiltration': None},
             'Rainfall Infiltration': {'Flood Mitigation': None},
             'Flood Mitigation': {'$ Value of Damages Aovided': None}
             }

g3 = pgv.AGraph(ex_nodes, directed = True)



g3.layout('dot')
g3.draw('figures/example.png')
