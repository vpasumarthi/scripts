# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np


def get_unit_cell_indices(system_size, total_elements_per_unit_cell, n_traj,
                          occupancy_data):
    """Returns the unit cell indices of the element
    :param system_size:
    :param total_elements_per_unit_cell:
    :param system_element_index:
    :return:
    """
    unit_cell_index_data = {}
    for traj_index in range(n_traj):
        traj_occupancy_data = occupancy_data[traj_index+1]
        (num_states, num_species) = traj_occupancy_data.shape
        traj_unit_cell_indices = np.zeros((num_states, num_species, 3), int)
        unit_cell_element_indices = traj_occupancy_data % total_elements_per_unit_cell
        total_filled_unit_cells = ((traj_occupancy_data - unit_cell_element_indices)
                                   // total_elements_per_unit_cell)
        for index in range(3):
            traj_unit_cell_indices[:, :, index] = total_filled_unit_cells / system_size[index+1:].prod()
            total_filled_unit_cells -= traj_unit_cell_indices[:, :, index] * system_size[index+1:].prod()
        unit_cell_index_data[traj_index+1] = traj_unit_cell_indices
    return unit_cell_index_data

def occupancy_analysis(src_path, occupancy, time_data, site_indices):
    unique_occupancies = np.unique(occupancy)
    bin_edges = np.hstack((unique_occupancies[0], unique_occupancies[:-1] + np.diff(unique_occupancies) / 2, unique_occupancies[-1]))
    hist = np.histogram(occupancy, bins=bin_edges, weights=time_data)[0]
    relative_residence_data = hist / np.sum(hist)
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
    print()
    
    bulk_shell_wise_occupancy_data = np.copy(shell_wise_occupancy_data[:4])
    bulk_shell_wise_occupancy_data[3, 1] = np.sum(shell_wise_occupancy_data[3:, 1])
    bulk_shell_wise_occupancy_data[3, 2] = np.min(shell_wise_occupancy_data[3:, 2])
    bulk_shell_wise_occupancy_data[3, 3] = np.min(shell_wise_occupancy_data[3:, 3])
    print('bulk_shell_wise_occupancy_data: (assume shell index 3. to be 3+)')
    print(bulk_shell_wise_occupancy_data)
    print()
    return None

def get_species_distribution(src_path, system_size, total_elements_per_unit_cell,
                             n_traj, gradient_ld):
    occupancy = np.load(src_path / 'occupancy.npy')[()][:-1]
    time = np.load(src_path / 'time_data.npy')[()]
    time_data = np.diff(time)[:, None]
    site_indices = np.load(src_path / 'site_indices.npy')[()]
    print('full system:')
    occupancy_analysis(src_path, occupancy, time_data, site_indices)
    
    occupancy_data = {1: occupancy}
    unit_cell_index_data = get_unit_cell_indices(
            system_size, total_elements_per_unit_cell, n_traj, occupancy_data)
    gradient_unit_cell_index = unit_cell_index_data[1][:, 0, 0]
    site_indices_data = {1: site_indices[:, 0][:, None]}
    site_indices_unit_cell_index_data = get_unit_cell_indices(
            system_size, total_elements_per_unit_cell, n_traj, site_indices_data)
    gradient_site_indices_unit_cell_index = site_indices_unit_cell_index_data[1][:, 0, 0]
    
    first_segment_indices = np.where(gradient_unit_cell_index < 5)[0]
    first_segment_occupancy = occupancy[first_segment_indices]
    first_segment_time_data = time_data[first_segment_indices]
    gradient_unit_cell_index = unit_cell_index_data[1][:, 0, 0]
    first_segment_site_indices = site_indices[np.where(gradient_site_indices_unit_cell_index < 5)[0]]
    print()
    print('first segment:')
    occupancy_analysis(src_path, first_segment_occupancy, first_segment_time_data, first_segment_site_indices)

    second_segment_indices = np.where(gradient_unit_cell_index >= 5)[0]
    second_segment_occupancy = occupancy[second_segment_indices]
    second_segment_time_data = time_data[second_segment_indices]
    gradient_unit_cell_index = unit_cell_index_data[1][:, 0, 0]
    second_segment_site_indices = site_indices[np.where(gradient_site_indices_unit_cell_index >= 5)[0]]
    print()
    print('second segment:')
    occupancy_analysis(src_path, second_segment_occupancy, second_segment_time_data, second_segment_site_indices)
    return None
