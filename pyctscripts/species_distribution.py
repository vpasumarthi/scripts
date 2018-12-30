#!/usr/bin/env python

import numpy as np


def get_species_distribution(src_path):
    occupancy = np.load(src_path / 'occupancy.npy')[()]
    time = np.load(src_path / 'time_data.npy')[()]
    unique_occupancies = np.unique(occupancy)
    bin_edges = np.hstack((unique_occupancies[0], unique_occupancies[:-1] + np.diff(unique_occupancies) / 2, unique_occupancies[-1]))
    hist = np.histogram(occupancy[:-1], bins=bin_edges, weights=np.diff(time)[:, None])[0]
    relative_residence_data = hist / np.sum(hist)
    site_indices = np.load(src_path / 'site_indices.npy')[()]
    num_shells = np.max(site_indices[:, 2])
    intersection_indices = np.in1d(site_indices[:, 0], unique_occupancies)
    intersection_site_indices = site_indices[intersection_indices]
    site_occupancy_data = np.delete(intersection_site_indices, 1, 1)
    site_occupancy_data = np.hstack((site_occupancy_data, relative_residence_data[:, None]))
    # (shell_index, mean_relative_residence, min_relative_residence, max_relative_residence)
    shell_wise_occupancy_data = np.zeros((num_shells+1, 4))
    shell_wise_occupancy_data[:, 0] = np.arange(num_shells+1)
    for shell_index in range(num_shells+1):
        shell_occupancy_data = site_occupancy_data[np.where(site_occupancy_data[:, 1]==shell_index)[0], 2]
        if len(shell_occupancy_data) > 0:
            shell_wise_occupancy_data[shell_index, 1] = np.sum(shell_occupancy_data)
            shell_wise_occupancy_data[shell_index, 2] = np.min(shell_occupancy_data)
            shell_wise_occupancy_data[shell_index, 3] = np.max(shell_occupancy_data)
    np.set_printoptions(precision=4, suppress=True)
    print('shell_wise_occupancy_data:')
    print(shell_wise_occupancy_data)
    return None
