#!/usr/bin/env python

import numpy as np


def traj_shell_wise_residence(src_path, traj_index):
    # TODO: parse yaml file for output data file names
    site_indices_data = np.load(f'{src_path}/traj{traj_index}/site_indices.npy')[()]
    occupancy = np.load(f'{src_path}/traj{traj_index}/occupancy.npy')[()]
    time = np.load(f'{src_path}/traj{traj_index}/time_data.npy')[()]
    time_step_data = np.diff(time)

    # TODO: change the following hard-code for 2-shell implementation to any number of shells
    first_shell_site_indices = site_indices_data[site_indices_data[:, 2] == 0][:, 0]
    second_shell_site_indices = site_indices_data[site_indices_data[:, 2] == 1][:, 0]
    third_shell_site_indices = site_indices_data[site_indices_data[:, 2] == 2][:, 0]
    bulk_site_indices = site_indices_data[site_indices_data[:, 2] > 2][:, 0]

    num_first_shell_sites = len(first_shell_site_indices)
    num_second_shell_sites = len(second_shell_site_indices)
    num_third_shell_sites = len(third_shell_site_indices)
    num_bulk_sites = len(bulk_site_indices)
    shell_wise_num_sites = [num_first_shell_sites, num_second_shell_sites, num_third_shell_sites, num_bulk_sites]

    return None

def shell_wise_residence(src_path, n_traj, kBT, shell_wise_penalties):
    # TODO: parse yaml file to reduce number of inputs
    num_shells = len(shell_wise_penalties)
    relative_residence_data = np.zeros((n_traj, num_shells))
    for traj_index in range(n_traj):
        relative_residence_data[traj_index, :] = traj_shell_wise_residence(src_path, traj_index+1)
    return None
