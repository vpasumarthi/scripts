#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from PyCT import constants

def plot_process_analysis(disp_array_prec, max_hop_dist, bar_color, annotate,
                          dst_path):
    # analysis on choice among available processes
    [unique_hop_dist, counts_hops] = np.unique(disp_array_prec,
                                               return_counts=True)
    plt.switch_backend('Agg')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    hop_proc_indices = np.where((0 < unique_hop_dist) & (unique_hop_dist <= max_hop_dist))[0]
    xtick_items = ['%1.4f' % item for item in unique_hop_dist[hop_proc_indices]]
    plt.bar(hop_proc_indices, counts_hops[hop_proc_indices], align='center', alpha=0.5,
            edgecolor='black', color=bar_color)
    plt.xticks(hop_proc_indices, xtick_items, rotation='vertical')

    if annotate:
        for i, v in enumerate(counts_hops[hop_proc_indices]):
            ax.text(i + 0.8, v + 100, str(v), color='green', rotation='vertical',
                    fontweight='bold')

    ax.set_xlabel('Hopping Distance')
    ax.set_ylabel('Frequency')
    ax.set_yscale('log')
    ax.set_title('Histogram of processes')
    filename = 'process_histogram'
    figure_name = filename + '.png'
    figure_path = dst_path / figure_name
    plt.tight_layout()
    plt.savefig(str(figure_path))
    return (counts_hops, hop_proc_indices)

def plot_escape_dist_analysis(uni_escape_dist, escape_proc_indices,
                              escape_counts, bar_color, annotate, dst_path):
    # analysis on escape distances
    fig = plt.figure()
    ax = fig.add_subplot(111)
    escape_dist_indices = np.arange(len(uni_escape_dist[escape_proc_indices]))
    xtick_items = ['%1.4f' % item for item in uni_escape_dist[escape_proc_indices]]
    plt.bar(escape_dist_indices, escape_counts[escape_proc_indices], align='center', alpha=0.5,
            edgecolor='black', color=bar_color)
    plt.xticks(escape_dist_indices, xtick_items, rotation='vertical')

    if annotate:
        for i, v in enumerate(escape_counts[escape_proc_indices]):
            ax.text(i - 0.2, v, str(v), color='green', rotation='vertical',
                    fontweight='bold')
    ax.set_xlabel('Escape Distance')
    ax.set_ylabel('Frequency')
    ax.set_yscale('log')
    ax.set_title('Histogram of Escape Distances')
    filename = 'escape_distance_histogram'
    figure_name = filename + '.png'
    figure_path = dst_path / figure_name
    plt.tight_layout()
    plt.savefig(str(figure_path))
    return None

def plot_mobility_analysis(mobility_dist_array, max_hop_dist, bar_color,
                           annotate, dst_path):
    # analysis on hopping distance contributing to mobility
    [unique_mobil_hop_dist, counts_mobil_hops] = np.unique(mobility_dist_array,
                                                           return_counts=True)
    mobil_proc_indices = np.where((0 < unique_mobil_hop_dist) & (unique_mobil_hop_dist <= max_hop_dist))[0]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    mobil_dist_indices = np.arange(len(unique_mobil_hop_dist[mobil_proc_indices]))
    xtick_items = ['%1.4f' % item for item in unique_mobil_hop_dist[mobil_proc_indices]]
    plt.bar(mobil_dist_indices, counts_mobil_hops[mobil_proc_indices], align='center', alpha=0.5,
            edgecolor='black', color=bar_color)
    plt.xticks(mobil_dist_indices, xtick_items, rotation='vertical')

    if annotate:
        for i, v in enumerate(counts_mobil_hops[mobil_proc_indices]):
            ax.text(i - 0.2, v, str(v), color='green', rotation='vertical',
                    fontweight='bold')
    ax.set_xlabel('Hop Distance')
    ax.set_ylabel('Frequency')
    ax.set_yscale('log')
    ax.set_title('Histogram of Hop Distances contributing to mobility')
    filename = 'mobil_hop_distance_histogram'
    figure_name = filename + '.png'
    figure_path = dst_path / figure_name
    plt.tight_layout()
    plt.savefig(str(figure_path))
    return None

def plot_rattle_analysis(rattle_dist_array, bar_color, annotate, dst_path):
    # analysis on hopping distance contributing to rattling
    [unique_rattle_hop_dist, counts_rattle_hops] = np.unique(
                                        rattle_dist_array, return_counts=True)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rattle_dist_indices = np.arange(len(unique_rattle_hop_dist))
    xtick_items = ['%1.4f' % item for item in unique_rattle_hop_dist]
    plt.bar(rattle_dist_indices, counts_rattle_hops, align='center', alpha=0.5,
            edgecolor='black', color=bar_color)
    plt.xticks(rattle_dist_indices, xtick_items, rotation='vertical')

    if annotate:
        for i, v in enumerate(counts_rattle_hops):
            ax.text(i - 0.2, v, str(v), color='green', rotation='vertical',
                    fontweight='bold')
    ax.set_xlabel('Hop Distance')
    ax.set_ylabel('Frequency')
    ax.set_yscale('log')
    ax.set_title('Histogram of Hop Distances contributing to rattling')
    filename = 'rattle_hop_distance_histogram'
    figure_name = filename + '.png'
    figure_path = dst_path / figure_name
    plt.tight_layout()
    plt.savefig(str(figure_path))
    return None

def traj_analysis(dst_path, intra_poly_dist_list, max_hop_dist, disp_prec,
                  annotate, bar_color):
    #NOTE: currently works with unwrapped_traj.dat which has positions at every
    # step written to it using 'write_every_step' branch of PyCT
    position_array = np.load(dst_path / 'traj1/unwrapped_traj.npy') / constants.ANG2BOHR

    disp_vec_array = np.diff(position_array, axis=0)
    disp_array = np.linalg.norm(disp_vec_array, axis=1)
    # round displacements to given precision
    disp_array_prec = np.round(disp_array, disp_prec)
    
    num_steps = position_array.shape[0] - 1
    num_rattles = 0
    rattle_dist_list = []
    rattle_event_list = []
    mobility_dist_list = []
    for step_index in range(num_steps):
        hop_dist = disp_array_prec[step_index]
        if hop_dist in intra_poly_dist_list:
            num_rattles += 1
            if num_rattles == 2:
                hop_dist_old = disp_array_prec[step_index - 1]
                rattle_dist_list.append(hop_dist_old)
                rattle_dist_list.append(hop_dist)
            elif num_rattles > 1:
                rattle_dist_list.append(hop_dist)
        else:
            if num_rattles == 1:
                mobility_dist_list.append(disp_array_prec[step_index - 1])
            elif num_rattles > 1:
                escape_dist = hop_dist
                rattle_event_list.append([num_rattles, escape_dist])
            mobility_dist_list.append(hop_dist)
            num_rattles = 0

    rattle_event_array = np.asarray(rattle_event_list)
    rattle_dist_array = np.asarray(rattle_dist_list)
    mobility_dist_array = np.asarray(mobility_dist_list)

    (counts_hops, hop_proc_indices) = plot_process_analysis(
                                                disp_array_prec, max_hop_dist,
                                                bar_color, annotate, dst_path)

    # report
    report_file_name = 'traj_analysis.log'
    num_kmc_steps = sum(counts_hops[hop_proc_indices])
    total_rattle_steps = int(np.sum(rattle_event_array[:, 0]))
    [uni_escape_dist, escape_counts] = np.unique(rattle_event_array[:, 1],
                                                 return_counts=True)
    escape_proc_indices = np.where((0 < uni_escape_dist) & (uni_escape_dist <= max_hop_dist))[0]
    escape_dist_list = list(uni_escape_dist[escape_proc_indices])
    with open(report_file_name, 'w') as report_file:
        report_file.write(f'Total number of kmc steps in simulation: '
                          f'{num_kmc_steps}\n')
        report_file.write(f'Cumulative number of kmc steps in rattling: '
                          f'{total_rattle_steps}\n')
        report_file.write(
                        f'Average number of rattles per rattle event:'
                        f'{np.mean(rattle_event_array[:, 0]):{7}.{5}}\n')
        report_file.write(
                    f'List of escape distances: '
                    f'{", ".join(str(dist) for dist in escape_dist_list)}\n')

    plot_escape_dist_analysis(uni_escape_dist, escape_proc_indices,
                              escape_counts, bar_color, annotate, dst_path)
    plot_mobility_analysis(mobility_dist_array, max_hop_dist, bar_color,
                           annotate, dst_path)
    plot_rattle_analysis(rattle_dist_array, bar_color, annotate, dst_path)
    return None
