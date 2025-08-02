# utils.py
import numpy as np

def compute_boresight_angle(sat_pos, user_xy, beam_center_xy):
    vec_beam = np.hstack([beam_center_xy, 0]) - sat_pos
    vec_user = np.hstack([user_xy, 0]) - sat_pos
    # Angle between vectors (dot product)
    cos_theta = np.dot(vec_beam, vec_user) / (np.linalg.norm(vec_beam) * np.linalg.norm(vec_user))
    theta = np.arccos(np.clip(cos_theta, -1., 1.))
    return theta
