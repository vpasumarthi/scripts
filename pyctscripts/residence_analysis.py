#!/usr/bin/env python

import numpy as np


def traj_shell_wise_residence(src_path, traj_index):
    # TODO: parse yaml file for output data file names
    site_indices_data = np.load(f'{src_path}/traj{traj_index}/site_indices.npy')[()]
    occupancy = np.load(f'{src_path}/traj{traj_index}/occupancy.npy')[()]
    time = np.load(f'{src_path}/traj{traj_index}/time_data.npy')[()]
    time_step_data = np.diff(time)

    return None

def shell_wise_residence(src_path, n_traj, kBT, shell_wise_penalties):
    # TODO: parse yaml file to reduce number of inputs
    num_shells = len(shell_wise_penalties)
    relative_residence_data = np.zeros((n_traj, num_shells))
    for traj_index in range(n_traj):
        relative_residence_data[traj_index, :] = traj_shell_wise_residence(src_path, traj_index+1)
    return None
