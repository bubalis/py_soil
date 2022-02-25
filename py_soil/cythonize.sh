#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 17:53:17 2021

@author: bdube
"""

cd /home/bdube/miniconda3
conda init bash
conda activate soilES
cd /mnt/c/Users/benja/Documents/soil_es/runoff_erosion
python setup.py build_ext --inplace
