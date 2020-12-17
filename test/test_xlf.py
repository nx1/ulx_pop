# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 14:16:57 2020

@author: norma
"""

import sys
sys.path.insert(0, '../src')

import pytest
import populations


import xlf


@pytest.fixture
def pop():
    df = populations.startrack_v2_mt_1_all(nrows=1000)
    pop = populations.Population(df) 
    return pop

def test_():
    pass

def test_main(pop):
    xlf.main(pop, N=100, N_mc=100)
#    import matplotlib.pyplot as plt
#    plt.show()
    
