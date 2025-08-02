# config.py
import numpy as np

# System parameters
N = 1             # satellites
M = 106           # beams per satellite
K = 1058          # ground cells
U = 5_290_000     # user equipments
sat_height = 600e3  # 600 km in meters

# Hex grid parameters
cell_radius = 2296  # meters, adjust as needed
num_cols = 34       # Approximates K=1058 when used with num_rows
num_rows = 32

# Frequency (Hz)
carrier_freq = 20e9

# Shadowed-Rician parameters (example, refer to recommended channel settings)
k_rician = 2.4     # Rician K-factor
md = 1.5           # Nakagami parameter for dominant
ms = 1.5           # Nakagami parameter for scattered
Omega = 2.0  # Average LoS power (can be adjusted based on scenario)
b = 1.0      # Multipath average power (can be adjusted)
threshold = 1.0  # Channel gain threshold used for outage probability

# Antenna parameters (example TR 38.901, ref [6])
Gmax = 30         # dBi (antenna gain peak, example value)
theta_3db = np.deg2rad(7)   # 3 dB beamwidth in radians

# Channel parameters


