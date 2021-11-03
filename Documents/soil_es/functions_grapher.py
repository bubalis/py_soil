#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 18:25:33 2021

@author: bdube
"""
import re

def split_and_strip(string, splitter = ','):
    return [p.strip() for p in string.split(splitter)]

def graphviz_code(text):
    funcnames = re.findall('(?<=def )[^\s^\\(]+', text)
    func_args = [split_and_strip(re.search(f'def {name}\s*\((.*)\)', text).group(1))
                 for name in funcnames]
    func_outputs = [split_and_strip(r) for r in re.findall('(?<=return ).+\n', text)]
    
    out_dict = {}
    for i, arglist in enumerate(func_args):
        for arg in arglist:
            if arg not in out_dict:
                out_dict[arg] = []
            out_dict[arg] += func_outputs[i]
    return out_dict        
        

path = '/mnt/c/Users/benja/Documents/soil_es/runoff_MUSLE.py'

text = open(path).read()
