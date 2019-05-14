#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  8 10:51:53 2019

@author: nk7g14
"""

from multiprocessing import Pool
import numpy as np
import time
    
p = Pool(5)
    
def square(x):
    temp = [i**2 for i in range(10000-x)]
    print(sum(temp))
    return x*x

mylist = np.arange(0,100,1)
emptylist = np.empty(len(mylist))

'''
print('Performing non-multi method')
time0 = time.time()

for i in mylist:
    emptylist[i] = square(i)
    
print('time taken:', time.time() - time0)
'''

print('Performing multi method')
time1 = time.time()
mymap = p.map(square, mylist)
print('time taken:', time.time() - time1)





