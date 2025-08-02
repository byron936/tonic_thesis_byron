# channel.py
import numpy as np
from scipy.constants import c, pi
from scipy.special import gamma, gammainc, factorial
from config import sat_height, carrier_freq, k_rician, md, ms

def free_space_path_loss(sat_pos, cell_centers):
    wavelength = c / carrier_freq
    d = np.sqrt((sat_pos[0] - cell_centers[:, 0])**2 +
                (sat_pos[1] - cell_centers[:, 1])**2 +
                sat_pos[2]**2)
    loss = (wavelength / (4 * pi * d)) ** 2
    return loss

def shadowed_rician_cdf(x, m, Omega, b, n_max=30):
    # Pre-compute constants
    two_b = 2 * b
    K = (two_b * m / (two_b * m + Omega))**m / two_b
    delta = (Omega / (two_b * m + Omega)) / two_b
    sum_terms = 0.0
    for n in range(n_max + 1):
        poch = gamma(m + n) / gamma(m)
        coeff = poch * (delta**n) * (two_b ** (1 + n)) / ((factorial(n))**2)
        term = coeff * gamma_lower(1 + n, x / two_b)
        sum_terms += term
    cdf = K * sum_terms
    return cdf

def gamma_lower(a, x):
    # Lower incomplete gamma: gamma(a, x) = gammainc(a, x) * gamma(a)
    return gammainc(a, x) * gamma(a)
