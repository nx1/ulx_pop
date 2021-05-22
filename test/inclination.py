import numpy as np
import matplotlib.pyplot as plt


# Randomly generated inclinations in radians assuming uniform spacial distribution in two dimensions
def sample_inclination_sin():
    i = np.arccos(np.random.random())
    return i

def test_bulk_arccos():
    i = np.arccos(np.random.random(size=1000000)) * 180 / np.pi
    assert max(i) < 91
    assert min(i) >= 0
    assert max(i) < 91
    assert min(i) >= 0

def P_i(i):
    """i=inclination in deg, (integer values)"""
    i_upper = i+0.5
    i_lower = i-0.5
    i_upper_rad = i_upper*np.pi/180
    i_lower_rad = i_lower*np.pi/180
    P_i = -np.cos(i_upper_rad) --np.cos(i_lower_rad)
    return P_i



i = np.array([sample_inclination_sin() for n in range(100000)])
i_deg = i * 180 / np.pi


P_is = [P_i(i) for i in range(91)]

plt.figure()
plt.plot(range(91), P_is, label='P(i) calculated')
plt.hist(i_deg, bins=100, density=True, label='P(i) via Monte-carlo')
plt.xlabel('inclination')
plt.ylabel('P(i)')
plt.legend()
plt.show()

