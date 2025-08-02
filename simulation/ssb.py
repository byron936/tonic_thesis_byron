# ssb.py
import numpy as np

def assign_ssb_periodicity(K):
    # Each cell SSB periodicity: randomly assign as per spec
    return np.random.choice([20, 40, 80, 160], K)

def ue_random_access_delay(users_cell_idx, ssb_periodicity, N_UE):
    # Initial wait: uniform over periodicity
    initial_wait = np.random.uniform(0, ssb_periodicity[users_cell_idx], N_UE)
    return initial_wait
