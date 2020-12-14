"""
This file contains some translated routines I've used in C
but i needed to make sure they work properly so I've typed them up in python.
"""

import pytest

def get_first_transient_cycle(is_ulx):
    for i in range(1,8):
        diff = is_ulx[i] - is_ulx[i-1]
        if diff != 0:
            return i
    return -1

def test_first_transient_cycle_index():
    case0 = [0,0,0,0,0,0,0,0]   
    case1 = [1,1,1,1,1,1,1,1]  
    case2 = [1,0,0,0,0,0,0,0]
    case3 = [0,0,0,0,0,0,0,1]
    case4 = [0,1,0,1,0,1,0,1]
    case5 = [0,0,0,1,0,0,0,0]

    assert get_first_transient_cycle(case0) == -1
    assert get_first_transient_cycle(case1) == -1
    assert get_first_transient_cycle(case2) == 1
    assert get_first_transient_cycle(case3) == 7
    assert get_first_transient_cycle(case4) == 1
    assert get_first_transient_cycle(case5) == 3


