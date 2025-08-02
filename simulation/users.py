# users.py
import numpy as np
from config import K, U, cell_radius
from hexgrid import generate_hex_grid

def create_users():
    centers = generate_hex_grid()
    users_cell_idx = np.random.choice(K, U)
    sigma = cell_radius / 2.5
    offsets = np.random.normal(0, sigma, (U, 2))
    users_xy = centers[users_cell_idx] + offsets
    return users_xy, users_cell_idx
