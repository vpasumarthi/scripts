#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from PyCT import constants


def traj_analysis(dst_path, disp_prec):
    position_array = (np.loadtxt(dst_path.joinpath('unwrappedTraj.dat'))
                      / constants.ANG2BOHR)
    num_positions = len(position_array)
    desired_indices = [0]
    for stepIndex in range(1, num_positions):
        if not np.array_equal(
                        np.round(position_array[stepIndex, :], disp_prec),
                        np.round(position_array[stepIndex-1, :], disp_prec)):
            desired_indices.append(stepIndex)
    newposition_array = np.copy(position_array[desired_indices])

    disp_vec_array = np.diff(newposition_array, axis=0)
    disp_array = np.linalg.norm(disp_vec_array, axis=1)

    num_steps = len(newposition_array) - 1
    rattle_list = []
    num_rattles = 1
    for stepIndex in range(2, num_steps+1):
        if np.array_equal(
                    np.round(newposition_array[stepIndex, :], disp_prec),
                    np.round(newposition_array[stepIndex-2, :], disp_prec)):
            num_rattles += 1
        else:
            if num_rattles > 2:
                rattle_dist = np.linalg.norm(
                                        newposition_array[stepIndex-1, :]
                                        - newposition_array[stepIndex-2, :])
                escape_dist = np.linalg.norm(
                                        newposition_array[stepIndex, :]
                                        - newposition_array[stepIndex-1, :])
                rattle_list.append([np.round(rattle_dist, disp_prec),
                                    num_rattles,
                                    np.round(escape_dist, disp_prec)])
                num_rattles = 1
    rattle_array = np.asarray(rattle_list)
    exceptions = np.array([2.8541, 2.8600, 2.9958])

    polyhedral_rattle_list = []
    num_rattles = 0
    for rattle_index in range(len(rattle_array)):
        num_rattles += rattle_array[rattle_index, 1]
        escape_dist = rattle_array[rattle_index, 2]
        if escape_dist not in exceptions:
            polyhedral_rattle_list.append([num_rattles, escape_dist])
            num_rattles = 0
    polyhedral_rattle_array = np.asarray(polyhedral_rattle_list)

    # report
    report_file_name = 'traj_analysis.log'
    num_kmc_steps = disp_array.shape[0]
    total_rattle_steps = int(np.sum(polyhedral_rattle_array[:, 0]))
    [uni_escape_dist, escape_counts] = np.unique(polyhedral_rattle_array[:, 1],
                                                 return_counts=True)
    escape_dist_list = list(uni_escape_dist)
    with open(report_file_name, 'w') as report_file:
        report_file.write(f'Total number of kmc steps in simulation: '
                          f'{num_kmc_steps}\n')
        report_file.write(f'Cumulative number of kmc steps in rattling: '
                          f'{total_rattle_steps}\n')
        report_file.write(
                        f'Average number of rattles per rattle event:'
                        f'{np.mean(polyhedral_rattle_array[:, 0]):{7}.{5}}\n')
        report_file.write(
                    f'List of escape distances: '
                    f'{", ".join(str(dist) for dist in escape_dist_list)}\n')

    # analysis on choice among available processes
    # round displacements to given precision
    disp_array = np.round(disp_array, disp_prec)

    # collect to bins
    [unique, counts] = np.unique(disp_array, return_counts=True)

    plt.switch_backend('Agg')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    proc_indices = np.arange(len(unique))
    xtick_items = ['%1.4f' % item for item in unique]
    plt.bar(proc_indices, counts, align='center', alpha=0.5, edgecolor='black')
    plt.xticks(proc_indices, xtick_items, rotation='vertical')

    for i, v in enumerate(counts):
        ax.text(i - 0.2, v + 100, str(v), color='green', rotation='vertical',
                fontweight='bold')
    add_rectangle = 1
    if add_rectangle:
        ax.add_patch(patches.Rectangle((13.5, -10), 4, 500, fill=False,
                                       color='red'))
    ax.set_xlabel('Hopping Distance')
    ax.set_ylabel('Counts')
    ax.set_title('Histogram of processes')
    filename = 'process_histogram'
    figure_name = filename + '.png'
    figure_path = dst_path / figure_name
    plt.tight_layout()
    plt.savefig(str(figure_path))

    # analysis on escape distances
    fig = plt.figure()
    ax = fig.add_subplot(111)
    escape_dist_indices = np.arange(len(uni_escape_dist))
    xtick_items = ['%1.4f' % item for item in uni_escape_dist]
    plt.bar(escape_dist_indices, escape_counts, align='center', alpha=0.5,
            edgecolor='black')
    plt.xticks(escape_dist_indices, xtick_items, rotation='vertical')

    for i, v in enumerate(escape_counts):
        ax.text(i - 0.2, v, str(v), color='green', rotation='vertical',
                fontweight='bold')
    add_rectangle = 1
    if add_rectangle:
        ax.add_patch(patches.Rectangle((7.5, -10), 2, 25, fill=False,
                                       color='red'))
    ax.set_xlabel('Escape Distance')
    ax.set_ylabel('Counts')
    ax.set_title('Histogram of Escape Distances')
    filename = 'escape_distance_histogram'
    figure_name = filename + '.png'
    figure_path = dst_path / figure_name
    plt.tight_layout()
    plt.savefig(str(figure_path))
    return None
