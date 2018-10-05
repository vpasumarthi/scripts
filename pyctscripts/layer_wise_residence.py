#!/usr/bin/env python

import numpy as np


def read_occupancy(src_path, n_traj):
    """Reads the occupancy data from traj-level directories and return a 
       dictionary of trajectory-wise occupancy data, each as a numpy array
    :param src_path:
    :param n_traj:
    :return: occupancy_data:
    """
    occupancy_file_name = 'occupancy.npy'
    occupancy_data = {}
    for traj_index in range(n_traj):
        traj_dir_path = src_path / f'traj{traj_index+1}'
        occupancy_data[traj_index+1] = np.load(traj_dir_path / occupancy_file_name)
    return occupancy_data

def compute_layer_wise_residence(src_path, system_size, n_traj, gradient_ld,
                                 layer_length_ratio):
    """Returns the layer wise residence of charge carriers
    :param src_path:
    :param system_size:
    :param n_traj:
    :param gradient_ld:
    :param layer_length_ratio:
    :return:
    """
    occupancy_data = read_occupancy(src_path, n_traj)
    return None
