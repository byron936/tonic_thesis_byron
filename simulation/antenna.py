# antenna.py
import numpy as np
from scipy.special import j1, jn
from config import Gmax, theta_3db

def antenna_gain(theta):
    # Main lobe, Bessel function order 1 and 3, formula [6]
    mu = 2.07123 * np.sin(theta) / np.sin(theta_3db)
    with np.errstate(divide='ignore', invalid='ignore'):
        pattern = (j1(mu) / 2 * mu)**2 + 36 * (jn(3, mu) / mu**3)**2
        gain = Gmax * (pattern**2)
        gain = np.where(mu == 0, Gmax, gain)
    return gain
