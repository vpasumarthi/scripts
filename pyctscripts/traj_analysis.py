#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import yaml

from PyCT import constants

def generate_report(hop_dist_count_array, hop_proc_indices,
                    rattle_event_array_dict, max_hop_dist):
    report_file_name = f'traj_analysis.log'
    n_traj = len(rattle_event_array_dict)
    num_kmc_steps_array = np.zeros(n_traj, int)
    total_rattle_steps_array = np.zeros(n_traj, int)
    average_rattles_per_event_array = np.zeros(n_traj)
    escape_dist_list_array = np.empty(n_traj, object)
    traj_wise_escape_count_array = np.empty(n_traj, object)
    escape_proc_indices_array = np.empty(n_traj, object)
    for traj_index in range(n_traj):
        num_kmc_steps_array[traj_index] = sum(hop_dist_count_array[traj_index, :][hop_proc_indices])
        if len(rattle_event_array_dict[traj_index+1]):
            total_rattle_steps_array[traj_index] = int(np.sum(rattle_event_array_dict[traj_index+1][:, 0]))
            average_rattles_per_event_array[traj_index] = np.mean(rattle_event_array_dict[traj_index+1][:, 0])
            [uni_escape_dist, escape_counts] = np.unique(rattle_event_array_dict[traj_index+1][:, 1],
                                                         return_counts=True)
            escape_proc_indices = np.where((0 < uni_escape_dist) & (uni_escape_dist <= max_hop_dist))[0]
            escape_dist_list_array[traj_index] = np.copy(uni_escape_dist[escape_proc_indices])
            traj_wise_escape_count_array[traj_index] = np.copy(escape_counts)
            escape_proc_indices_array[traj_index] = np.copy(escape_proc_indices)
        else:
            total_rattle_steps_array[traj_index] = 0
            average_rattles_per_event_array[traj_index] = 0
            escape_dist_list_array[traj_index] = np.array([])
            traj_wise_escape_count_array[traj_index] = np.array([])
            escape_proc_indices_array[traj_index] = np.array([])
        if traj_index == 0:
            cumulative_escape_dist_list = np.copy(escape_dist_list_array[traj_index])
        else:
            cumulative_escape_dist_list = np.append(cumulative_escape_dist_list, escape_dist_list_array[traj_index])
    unique_escape_dist_array = np.unique(cumulative_escape_dist_list.round(4))
    with open(report_file_name, 'w') as report_file:
        report_file.write(f'Total number of kmc steps in simulation: {num_kmc_steps_array.mean():4.3e} +/- {num_kmc_steps_array.std() / np.sqrt(n_traj):4.3e}\n')
        report_file.write(f'Cumulative number of kmc steps in rattling: {total_rattle_steps_array.mean():4.3e} +/- {total_rattle_steps_array.std() / np.sqrt(n_traj):4.3e}\n')
        report_file.write(f'Average number of rattles per rattle event: {average_rattles_per_event_array.mean():4.3f} +/- {average_rattles_per_event_array.std() / np.sqrt(n_traj):4.3f}\n')
        report_file.write(f'List of escape distances: {", ".join(str(dist) for dist in unique_escape_dist_array)}\n')
    return (escape_dist_list_array, escape_proc_indices_array, traj_wise_escape_count_array)

def plot_process_analysis(disp_array_prec_dict, max_hop_dist, xlabel_choice,
                          dist_to_barrier_height_dict, bar_color, annotate,
                          dst_path, plot_style):
    n_traj = len(disp_array_prec_dict)
    for traj_index in range(n_traj):
        if traj_index == 0:
            cumulative_hop_dist_array = np.copy(disp_array_prec_dict[traj_index+1])
        else:
            cumulative_hop_dist_array = np.append(cumulative_hop_dist_array, disp_array_prec_dict[traj_index+1])
    cumulative_unique_hop_dist = np.unique(cumulative_hop_dist_array)
    num_unique_hop_dist = len(cumulative_unique_hop_dist)
    hop_dist_count_array = np.zeros((n_traj, num_unique_hop_dist), int)

    # analysis on choice among available processes
    for traj_index in range(n_traj):
        [unique_hop_dist, hop_count] = np.unique(disp_array_prec_dict[traj_index+1],
                                                        return_counts=True)
        for hop_dist_index, hop_dist in enumerate(unique_hop_dist):
            dest_index = np.where(cumulative_unique_hop_dist == hop_dist)[0][0]
            hop_dist_count_array[traj_index, dest_index] = hop_count[hop_dist_index]

    mean_hop_count = np.mean(hop_dist_count_array, axis=0)
    sem_hop_count = np.std(hop_dist_count_array, axis=0) / np.sqrt(n_traj)

    plt.switch_backend('Agg')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    hop_proc_indices = np.where((0 < cumulative_unique_hop_dist) & (cumulative_unique_hop_dist <= max_hop_dist))[0]
    if xlabel_choice == 'hop_dist':
        xtick_items = ['%1.4f' % item for item in cumulative_unique_hop_dist[hop_proc_indices]]
    elif xlabel_choice == 'activation_energy':
        xtick_items = ['%1.4f' % dist_to_barrier_height_dict[item] for item in cumulative_unique_hop_dist[hop_proc_indices]]
    plt.errorbar(hop_proc_indices, mean_hop_count[hop_proc_indices],
                 yerr=sem_hop_count[hop_proc_indices], fmt='o', capsize=3,
                 color=bar_color, mfc='none', mec='none')
    plt.bar(hop_proc_indices, mean_hop_count[hop_proc_indices], align='center', alpha=1,
            edgecolor='black', color=bar_color)
    plt.xticks(hop_proc_indices, xtick_items, rotation='vertical')

    if annotate:
        for i, v in enumerate(mean_hop_count[hop_proc_indices]):
            ax.text(i + 0.8, v + 100, str(v), color='green', rotation='vertical',
                    fontweight='bold')

    ax.set_xlabel('Hopping Distance')
    ax.set_ylabel('Frequency')
    ax.set_yscale(plot_style)
    ax.set_title('Histogram of processes')
    filename = f'process_histogram_{plot_style}'
    figure_name = filename + '.png'
    figure_path = dst_path / figure_name
    plt.tight_layout()
    plt.savefig(str(figure_path), dpi=600)
    return (hop_dist_count_array, hop_proc_indices)

def plot_escape_dist_analysis(escape_dist_list_array, escape_proc_indices_array,
                              traj_wise_escape_count_array, xlabel_choice,
                              dist_to_barrier_height_dict, bar_color, annotate,
                              dst_path, plot_style):
    # analysis on escape distances
    n_traj = len(escape_dist_list_array)
    for traj_index in range(n_traj):
        if traj_index == 0:
            cumulative_escape_dist_list = np.copy(escape_dist_list_array[traj_index])
        else:
            cumulative_escape_dist_list = np.append(cumulative_escape_dist_list, escape_dist_list_array[traj_index])
    unique_escape_dist = np.unique(cumulative_escape_dist_list)
    num_unique_escape_dist = len(unique_escape_dist)
    escape_dist_escape_count_array = np.zeros((n_traj, num_unique_escape_dist), int)
    for traj_index in range(n_traj):
        traj_escape_dist_count = traj_wise_escape_count_array[traj_index][escape_proc_indices_array[traj_index]]
        for source_index, escape_dist in enumerate(escape_dist_list_array[traj_index]):
            dest_index = np.where(unique_escape_dist == escape_dist)[0][0]
            escape_dist_escape_count_array[traj_index, dest_index] = traj_escape_dist_count[source_index]
    mean_escape_count_array = np.mean(escape_dist_escape_count_array, axis=0)
    sem_escape_count_array = np.std(escape_dist_escape_count_array, axis=0) / np.sqrt(n_traj)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    escape_dist_indices = np.arange(num_unique_escape_dist)
    if xlabel_choice == 'hop_dist':
        xtick_items = ['%1.4f' % item for item in unique_escape_dist]
    elif xlabel_choice == 'activation_energy':
        xtick_items = ['%1.4f' % dist_to_barrier_height_dict[item] for item in unique_escape_dist]
    plt.errorbar(escape_dist_indices, mean_escape_count_array,
                 yerr=sem_escape_count_array, fmt='o', capsize=3,
                 color=bar_color, mfc='none', mec='none')
    plt.bar(escape_dist_indices, mean_escape_count_array, align='center', alpha=1,
            edgecolor='black', color=bar_color)
    plt.xticks(escape_dist_indices, xtick_items, rotation='vertical')

    if annotate:
        for i, v in enumerate(mean_escape_count_array):
            ax.text(i - 0.2, v, str(v), color='green', rotation='vertical',
                    fontweight='bold')
    ax.set_xlabel('Escape Distance')
    ax.set_ylabel('Frequency')
    ax.set_yscale(plot_style)
    ax.set_title('Histogram of Escape Distances')
    filename = f'escape_distance_histogram_{plot_style}'
    figure_name = filename + '.png'
    figure_path = dst_path / figure_name
    plt.tight_layout()
    plt.savefig(str(figure_path), dpi=600)
    return None

def plot_mobility_analysis(mobility_dist_array_dict, max_hop_dist,
                           xlabel_choice, dist_to_barrier_height_dict,
                           bar_color, annotate, dst_path, plot_style):
    n_traj = len(mobility_dist_array_dict)
    # analysis on hopping distance contributing to mobility
    for traj_index in range(n_traj):
        if traj_index == 0:
            cumulative_mobil_dist_array = np.copy(mobility_dist_array_dict[traj_index+1])
        else:
            cumulative_mobil_dist_array = np.append(cumulative_mobil_dist_array, mobility_dist_array_dict[traj_index+1])
    cumulative_unique_mobil_hop_dist = np.unique(cumulative_mobil_dist_array)
    num_unique_mobil_hop_dist = len(cumulative_unique_mobil_hop_dist)
    mobil_dist_hop_count_array = np.zeros((n_traj, num_unique_mobil_hop_dist), int)
    for traj_index in range(n_traj):
        [unique_mobil_hop_dist, counts_mobil_hops] = np.unique(mobility_dist_array_dict[traj_index+1],
                                                               return_counts=True)
        mobil_proc_indices = np.where((0 < unique_mobil_hop_dist) & (unique_mobil_hop_dist <= max_hop_dist))[0]
        for source_index, mobil_hop_dist in enumerate(unique_mobil_hop_dist[mobil_proc_indices]):
            dest_index = np.where(cumulative_unique_mobil_hop_dist == mobil_hop_dist)[0][0]
            mobil_dist_hop_count_array[traj_index, dest_index] = counts_mobil_hops[mobil_proc_indices][source_index]
    mean_mobil_hop_count_array = np.mean(mobil_dist_hop_count_array, axis=0)
    sem_mobil_hop_count_array = np.std(mobil_dist_hop_count_array, axis=0) / np.sqrt(n_traj)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    mobil_dist_indices = np.arange(num_unique_mobil_hop_dist)
    if xlabel_choice == 'hop_dist':
        xtick_items = ['%1.4f' % item for item in cumulative_unique_mobil_hop_dist]
    elif xlabel_choice == 'activation_energy':
        xtick_items = ['%1.4f' % dist_to_barrier_height_dict[item] for item in cumulative_unique_mobil_hop_dist]
    plt.errorbar(mobil_dist_indices, mean_mobil_hop_count_array,
                 yerr=sem_mobil_hop_count_array, fmt='o', capsize=3,
                 color=bar_color, mfc='none', mec='none')
    plt.bar(mobil_dist_indices, mean_mobil_hop_count_array, align='center', alpha=1,
            edgecolor='black', color=bar_color)
    plt.xticks(mobil_dist_indices, xtick_items, rotation='vertical')

    if annotate:
        for i, v in enumerate(counts_mobil_hops[mobil_proc_indices]):
            ax.text(i - 0.2, v, str(v), color='green', rotation='vertical',
                    fontweight='bold')
    ax.set_xlabel('Hop Distance')
    ax.set_ylabel('Frequency')
    ax.set_yscale(plot_style)
    ax.set_title('Histogram of Hop Distances contributing to mobility')
    filename = f'mobil_hop_distance_histogram_{plot_style}'
    figure_name = filename + '.png'
    figure_path = dst_path / figure_name
    plt.tight_layout()
    plt.savefig(str(figure_path), dpi=600)
    return None

def plot_rattle_analysis(rattle_dist_array_dict, xlabel_choice,
                         dist_to_barrier_height_dict, bar_color, annotate,
                         dst_path, plot_style):
    n_traj = len(rattle_dist_array_dict)
    # analysis on hopping distance contributing to rattling
    for traj_index in range(n_traj):
        if traj_index == 0:
            cumulative_rattle_hop_dist_array = np.copy(rattle_dist_array_dict[traj_index+1])
        else:
            cumulative_rattle_hop_dist_array = np.append(cumulative_rattle_hop_dist_array, rattle_dist_array_dict[traj_index+1])
    cumulative_unique_rattle_hop_dist = np.unique(cumulative_rattle_hop_dist_array)
    num_unique_rattle_hop_dist = len(cumulative_unique_rattle_hop_dist)
    rattle_dist_hop_count_array = np.zeros((n_traj, num_unique_rattle_hop_dist), int)
    for traj_index in range(n_traj):
        [unique_rattle_hop_dist, counts_rattle_hops] = np.unique(
                    rattle_dist_array_dict[traj_index+1], return_counts=True)
        for source_index, rattle_hop_dist in enumerate(unique_rattle_hop_dist):
            dest_index = np.where(cumulative_unique_rattle_hop_dist == rattle_hop_dist)[0][0]
            rattle_dist_hop_count_array[traj_index, dest_index] = counts_rattle_hops[source_index]
    mean_rattle_dist_hop_count = np.mean(rattle_dist_hop_count_array, axis=0)
    sem_rattle_dist_hop_count = np.std(rattle_dist_hop_count_array, axis=0) / np.sqrt(n_traj)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rattle_dist_indices = np.arange(num_unique_rattle_hop_dist)
    if xlabel_choice == 'hop_dist':
        xtick_items = ['%1.4f' % item for item in cumulative_unique_rattle_hop_dist]
    elif xlabel_choice == 'activation_energy':
        xtick_items = ['%1.4f' % dist_to_barrier_height_dict[item] for item in cumulative_unique_rattle_hop_dist]
    plt.errorbar(rattle_dist_indices, mean_rattle_dist_hop_count,
                 yerr=sem_rattle_dist_hop_count, fmt='o', capsize=3,
                 color=bar_color, mfc='none', mec='none')
    plt.bar(rattle_dist_indices, mean_rattle_dist_hop_count, align='center', alpha=1,
            edgecolor='black', color=bar_color)
    plt.xticks(rattle_dist_indices, xtick_items, rotation='vertical')

    if annotate:
        for i, v in enumerate(counts_rattle_hops):
            ax.text(i - 0.2, v, str(v), color='green', rotation='vertical',
                    fontweight='bold')
    ax.set_xlabel('Hop Distance')
    ax.set_ylabel('Frequency')
    ax.set_yscale(plot_style)
    ax.set_title('Histogram of Hop Distances contributing to rattling')
    filename = f'rattle_hop_distance_histogram_{plot_style}'
    figure_name = filename + '.png'
    figure_path = dst_path / figure_name
    plt.tight_layout()
    plt.savefig(str(figure_path), dpi=600)
    return None

def traj_analysis(dst_path, rattle_distance_pool, rattle_definition,
                  max_hop_dist, disp_prec, xlabel_choice,
                  dist_to_barrier_height_dict, annotate, bar_color, plot_style):
    #NOTE: currently works with unwrapped_traj.dat which has positions at every
    # step written to it using 'write_every_step' branch of PyCT

    # Load simulation parameters
    sim_param_file_name = 'simulation_parameters.yml'
    sim_param_file_path = dst_path / sim_param_file_name
    with open(sim_param_file_path, 'r') as stream:
        try:
            sim_params = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    n_traj = int(sim_params['n_traj'])
    disp_array_prec_dict = {}
    rattle_dist_array_dict = {}
    rattle_event_array_dict = {}
    mobility_dist_array_dict = {}
    for traj_index in range(n_traj):
        position_array = np.load(f'{dst_path}/traj{traj_index+1}/unwrapped_traj.npy') / constants.ANG2BOHR
    
        disp_vec_array = np.diff(position_array, axis=0)
        disp_array = np.linalg.norm(disp_vec_array, axis=1)
        # round displacements to given precision
        disp_array_prec = np.round(disp_array, disp_prec)
        
        num_steps = position_array.shape[0] - 1
        num_rattles = 0
        rattle_dist_list = []
        rattle_event_list = []
        mobility_dist_list = []
        if rattle_definition == 'inclusive':
            for step_index in range(num_steps):
                hop_dist = disp_array_prec[step_index]
                if hop_dist in rattle_distance_pool:
                    num_rattles += 1
                    if num_rattles == 2:
                        hop_dist_old = disp_array_prec[step_index - 1]
                        rattle_dist_list.append(hop_dist_old)
                        rattle_dist_list.append(hop_dist)
                    elif num_rattles > 2:
                        rattle_dist_list.append(hop_dist)
                else:
                    if num_rattles == 1:
                        mobility_dist_list.append(disp_array_prec[step_index - 1])
                    elif num_rattles > 1:
                        escape_dist = hop_dist
                        rattle_event_list.append([num_rattles, escape_dist])
                    mobility_dist_list.append(hop_dist)
                    num_rattles = 0
        elif rattle_definition == 'exclusive':
            hop_dist_old = disp_array_prec[0]
            num_rattles = 1
            for step_index in range(1, num_steps):
                hop_dist_new = disp_array_prec[step_index]
                if hop_dist_new == hop_dist_old:
                    num_rattles += 1
                else:
                    if num_rattles > 1:
                        rattle_dist_list.extend([hop_dist_old] * num_rattles)
                        escape_dist = hop_dist_new
                        rattle_event_list.append([num_rattles, escape_dist])
                        num_rattles = 1
                    mobility_dist_list.append(hop_dist_old)
                hop_dist_old = hop_dist_new

        rattle_dist_array = np.asarray(rattle_dist_list)
        rattle_event_array = np.asarray(rattle_event_list)
        mobility_dist_array = np.asarray(mobility_dist_list)
        rattle_dist_array_dict[traj_index+1] = rattle_dist_array
        rattle_event_array_dict[traj_index+1] = rattle_event_array
        mobility_dist_array_dict[traj_index+1] = mobility_dist_array
        disp_array_prec_dict[traj_index+1] = disp_array_prec

    (hop_dist_count_array, hop_proc_indices) = plot_process_analysis(
                            disp_array_prec_dict, max_hop_dist, xlabel_choice,
                            dist_to_barrier_height_dict, bar_color, annotate,
                            dst_path, plot_style)
    (escape_dist_list_array, escape_proc_indices_array,
     traj_wise_escape_count_array) = generate_report(hop_dist_count_array,
                                                     hop_proc_indices,
                                                     rattle_event_array_dict,
                                                     max_hop_dist)
    len_escape_dist_list_array = [len(traj_escape_dist_list_array) for traj_escape_dist_list_array in escape_dist_list_array]
    if np.sum(len_escape_dist_list_array):
        plot_escape_dist_analysis(escape_dist_list_array, escape_proc_indices_array,
                                  traj_wise_escape_count_array, xlabel_choice,
                                  dist_to_barrier_height_dict, bar_color,
                                  annotate, dst_path, plot_style)
    plot_mobility_analysis(mobility_dist_array_dict, max_hop_dist, xlabel_choice,
                           dist_to_barrier_height_dict, bar_color, annotate,
                           dst_path, plot_style)
    plot_rattle_analysis(rattle_dist_array_dict, xlabel_choice,
                         dist_to_barrier_height_dict, bar_color, annotate,
                         dst_path, plot_style)
    return None
