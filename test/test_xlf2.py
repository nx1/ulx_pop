# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 14:16:57 2020

@author: norma
"""

import sys
import os
sys.path.insert(0, '../src')
os.chdir('../src')

import pytest
import populations
from constants import set_latex_font

import xlf2


@pytest.fixture
def pop():
    df = populations.startrack_v2_mt_1_all(nrows=10000)
    pop = populations.Population(df)
    return pop

def test_():
    pass

def test_main(pop):
    # xlf.main(pop, N=100, N_mc=10)
#    import matplotlib.pyplot as plt
#    plt.show()
    pass
