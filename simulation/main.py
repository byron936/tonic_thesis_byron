# main.py

import numpy as np
import matplotlib.pyplot as plt

from config import N, M, K, U, sat_height, md, Omega, b, threshold
from hexgrid import generate_hex_grid
from users import create_users
from channel import free_space_path_loss, shadowed_rician_cdf
from antenna import antenna_gain
from ssb import assign_ssb_periodicity, ue_random_access_delay
from utils import compute_boresight_angle

# Generate hexagonal ground cell centers
cell_centers = generate_hex_grid()

# Define satellite position at zenith over grid center
region_side = cell_centers.max(axis=0)
sat_pos = np.array([region_side[0] / 2, region_side[1] / 2, sat_height])

# Assign users to cells
users_xy, users_cell_idx = create_users()

# Assign SSB periodicity to each cell
ssb_periodicity = assign_ssb_periodicity(K)

# Calculate free space path loss from satellite to each cell
pl_cells = free_space_path_loss(sat_pos, cell_centers)

# Example: Evaluate analytical channel gain CDF at threshold for each cell
cdf_per_cell = [shadowed_rician_cdf(threshold, md, Omega, b) for _ in range(K)]

# Pick a specific UE for inspection
example_user_idx = 0
user_xy = users_xy[example_user_idx]
cell_idx = users_cell_idx[example_user_idx]
cell_center = cell_centers[cell_idx]
pl = pl_cells[cell_idx]
cdf = cdf_per_cell[cell_idx]

# Compute antenna gain for this user
theta = compute_boresight_angle(sat_pos, user_xy, cell_center)
gain = antenna_gain(theta)

# Display results
print(f"User {example_user_idx} in cell {cell_idx}:")
print(f"  Cell center: {cell_center}")
print(f"  Path Loss: {pl:.2e}")
print(f"  Antenna Gain: {gain:.2e}")
print(f"  Channel Gain CDF at threshold {threshold}: {cdf:.5f}")

# Compute and print average random access delay
random_access_delays = ue_random_access_delay(users_cell_idx, ssb_periodicity, U)
print(f"Random access delay mean: {random_access_delays.mean():.2f} ms")

# === Visualization: Plot hexagonal grid ===
plt.figure(figsize=(8, 7))
plt.scatter(cell_centers[:, 0], cell_centers[:, 1], c='deepskyblue', s=30, edgecolors='k')
plt.scatter(user_xy[0], user_xy[1], color='red', label='Example User')
plt.scatter(cell_center[0], cell_center[1], color='orange', label='Serving Cell Center')
plt.title('Hexagonal Grid of Ground Cells')
plt.xlabel('x (meters)')
plt.ylabel('y (meters)')
plt.legend()
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()

# For non-interactive backend (headless environment): save plot to file
plt.savefig("hexgrid_cells.png")
print("Hex grid plot saved as 'hexgrid_cells.png'")
