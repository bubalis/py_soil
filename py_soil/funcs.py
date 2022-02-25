#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 13:09:59 2021

@author: bdube
"""



import math
import numpy as np
from math import exp, log, nan
from rosetta import rosetta, SoilData
from numbers import Number

from numpy import vectorize
from sympy.solvers import solve
from sympy import Symbol, Eq
from functools import partial



def run_df(func):
    '''
    Wrapper to automatically run a function row-wise on a dataframe with 
    column names which match the function argument names.

    Parameters
    ----------
    func : A function

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    
    
    def wrapper(df, standard_args = {}):
        n_args = func.__code__.co_argcount
        columns = func.__code__.co_varnames[:n_args]
        
        if standard_args:
            columns = [c for c in columns if c not in standard_args]
            f = partial(func, **standard_args)
            try:
                return vectorize(f)(**{col: df[col] for col in columns})
            except:
                print(columns)
                #print(standard_args)
                print(f)
                raise
                
        return vectorize(func)(*[df[col] for col in columns])
    
    return wrapper












