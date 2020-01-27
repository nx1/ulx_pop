#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 10:43:52 2020

@author: nk7g14
"""

#import curvefunc
import pandas as pd
import time
import glob


time0 = time.time()


with open('./new_curves/a7232574-74cb-4240-a5b1-2e0fbf5262ca.txt', 'r') as file:
    data = file.read().splitlines()
    simulation_info = data[-9:] 
    
# curve = load_curve_file('./new_curves/a7232574-74cb-4240-a5b1-2e0fbf5262ca.txt')
tt = time.time() - time0

print('time taken:', tt)
print('time for 1.2m:', 1.2E6*tt/60/60)