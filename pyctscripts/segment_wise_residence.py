#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


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

def compute_segment_wise_residence(src_path, system_size, total_elements_per_unit_cell,
                                   n_traj, gradient_ld, segment_length_ratio,
                                   segmentwise_num_dopants, num_acceptor_sites_per_unit_cell):
    """Returns the segment wise residence of charge carriers
    :param src_path:
    :param system_size:
    :param n_traj:
    :param gradient_ld:
    :param segment_length_ratio:
    :return:
    """
    occupancy_data = read_occupancy(src_path, n_traj)
    unit_cell_index_data = get_unit_cell_indices(
            system_size, total_elements_per_unit_cell, n_traj, occupancy_data)
    num_elemental_segments = np.sum(segment_length_ratio)
    elemental_segment_system_size = np.copy(system_size)
    elemental_segment_system_size[gradient_ld] //= num_elemental_segments
    num_segments = len(segment_length_ratio)
    bin_edges = [0]
    segmentwise_doping_level = []
    segment_wise_residence = np.zeros((n_traj, num_segments), int)
    for segment_index in range(num_segments):
        segment_system_size = elemental_segment_system_size * segment_length_ratio[segment_index]
        bin_edges.append(bin_edges[-1] + segment_system_size[gradient_ld])
        segmentwise_num_acceptor_sites = segment_system_size.prod() * num_acceptor_sites_per_unit_cell
        segmentwise_doping_level.append(segmentwise_num_dopants[segment_index] / segmentwise_num_acceptor_sites * 100)
    for traj_index in range(n_traj):
        trajwise_unit_cell_index_data = unit_cell_index_data[traj_index+1]
        segment_wise_residence[traj_index] = np.histogram(trajwise_unit_cell_index_data[:, :, gradient_ld], bin_edges)[0]
    mean_segment_wise_residence = np.mean(segment_wise_residence, axis=0)
    std_segment_wise_residence = np.std(segment_wise_residence, axis=0)
    
    plt.switch_backend('Agg')
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    segment_index_list = range(1, num_segments+1)
    ax1.bar(segment_index_list, mean_segment_wise_residence, width=0.8,
            color='#0504aa', alpha=0.7)
    ax1.errorbar(segment_index_list, mean_segment_wise_residence,
                 yerr=std_segment_wise_residence, fmt='ko', capsize=3, mfc='none',
                 mec='none')
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(segment_index_list)
    ax2.set_xticklabels(segmentwise_doping_level)
    ax2.set_xlabel('Doping level (%)')
    ax1.set_xticks(segment_index_list)
    ax1.set_xticklabels(segment_index_list)
    ax1.set_xlabel('Segment Index')
    ax1.set_ylabel('Frequency')
    plt.tight_layout()
    plt.savefig(str(src_path / 'Segment-wise Residence.png'))

    segment_wise_relative_residence = segment_wise_residence / np.sum(segment_wise_residence, axis=1)[:, None]
    mean_segment_wise_relative_residence = np.mean(segment_wise_relative_residence, axis=0)
    std_segment_wise_relative_residence = np.std(segment_wise_relative_residence, axis=0)
    
    plt.switch_backend('Agg')
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    segment_index_list = range(1, num_segments+1)
    ax1.plot(segment_index_list, mean_segment_wise_relative_residence, 'o-',
             c='#0504aa', mfc='#0504aa', mec='black')
    ax1.errorbar(segment_index_list, mean_segment_wise_relative_residence,
                 yerr=std_segment_wise_relative_residence, fmt='o', capsize=3,
                 c='#0504aa', mfc='none', mec='none')
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(segment_index_list)
    ax2.set_xticklabels(segmentwise_doping_level)
    ax2.set_xlabel('Doping level (%)')
    ax1.set_xticks(segment_index_list)
    ax1.set_xticklabels(segment_index_list)
    ax1.set_xlabel('Segment Index')
    ax1.set_ylabel('Relative Frequency')
    plt.tight_layout()
    plt.savefig(str(src_path / 'Segment-wise Relative Residence.png'))
    return None
