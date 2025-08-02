# hexgrid.py
import numpy as np
from config import cell_radius, num_rows, num_cols

def generate_hex_grid():
    width = 2 * cell_radius
    height = np.sqrt(3) * cell_radius
    horiz_spacing = width * 0.75
    vert_spacing = height
    centers = []
    for r in range(num_rows):
        for c in range(num_cols):
            x = c * horiz_spacing
            y = r * vert_spacing + (c % 2) * (vert_spacing / 2)
            centers.append((x, y))
    return np.array(centers)
