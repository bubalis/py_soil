#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 17:50:46 2021

@author: bdube
"""

from setuptools import setup
from Cython.Build import cythonize
import os
import Cython.Compiler.Options
import numpy
Cython.Compiler.Options.annotate = True


if __name__ == '__main__':
    #path = os.path.join(os.getcwd(), "dumb.pyx")
    path = os.path.join(os.getcwd(), "greenampt_cy.pyx")
    
    #assert os.path.exists(path)
    setup(
        ext_modules = cythonize(path, annotate = True),
        package_dir = {'runoff_erosion': ''},
        include_dirs=[numpy.get_include()]
    )
