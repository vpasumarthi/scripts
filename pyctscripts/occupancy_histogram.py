#!/usr/bin/env python

import re

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


class Occupancy(object):
    """Class definition to generate occupancy histogram files"""

    def __init__(self, src_path, color_list):
        # Load occupancy parameters
        self.color_list = color_list
        self.num_colors = len(self.color_list)
        self.src_path = src_path
        return None

    def generate_occupancy_histogram(self, shell_wise, site_wise, n_traj):
        for traj_index in range(n_traj):
            traj_number = traj_index + 1
            (num_shells,
             probe_indices,
             site_population_list) = self.read_trajectory_data(traj_number)
        if shell_wise:
            self.generate_shell_wise_occupancy(num_shells,
                                               site_population_list)
        if site_wise:
            self.generate_site_wise_occpancy(num_shells,
                                             probe_indices,
                                             site_population_list)
        return None

    def read_trajectory_data(self, traj_number):
        site_indices_dir_name = 'site_indices_data'
        site_indices_file_name = f'site_indices_{traj_number}.csv'
        site_indices_file_path = self.src_path / site_indices_dir_name / site_indices_file_name
        occupancy_dir_name = 'occupancy_data'
        occupancy_file_name = f'occupancy_{traj_number}.dat'
        occupancy_file_path = self.src_path / occupancy_dir_name / occupancy_file_name
        shell_indices_dict = {}
        with site_indices_file_path.open('r') as site_indices_file:
            for line in site_indices_file:
                split_elements = re.split(',|\n', line)
                shell_index = int(split_elements[3])
                site_index = int(split_elements[0])
                if shell_index in shell_indices_dict:
                    shell_indices_dict[shell_index].append(site_index)
                else:
                    shell_indices_dict[shell_index] = [site_index]
                
        num_shells = len(shell_indices_dict) - 1
        probe_indices = []
        site_population_list = []
        for shell_index in range(num_shells+1):
            if shell_index == num_shells + 1:
                probe_indices.append(sorted(
                    [item
                     for sublist in shell_indices_dict[num_shells+1:]
                     for item in sublist]))
            else:
                probe_indices.append(sorted(
                    [item for item in shell_indices_dict[shell_index]]))
            site_population_list.append(
                                    [0] * len(probe_indices[shell_index]))

        with occupancy_file_path.open('r') as occupancy_file:
            for line in occupancy_file:
                site_index = int(line.split('\n')[0])
                for shell_index in range(num_shells+1):
                    if site_index in probe_indices[shell_index]:
                        list_index = probe_indices[shell_index].index(site_index)
                        site_population_list[shell_index][list_index] += 1
                        break
        return (num_shells, probe_indices, site_population_list)

    def generate_site_wise_occpancy(self, num_shells, probe_indices,
                                    site_population_list):
        plt.switch_backend('Agg')
        fig = plt.figure()
        plt.title('Site-wise occupancy')
        plt.axis('off')
        num_sub_plots = num_shells + 1
        num_cols = 2
        num_rows = num_sub_plots // num_cols + num_sub_plots % num_cols
        num_data = 0
        gs = gridspec.GridSpec(num_rows, num_cols, hspace=0.4)
        for shell_index in range(num_sub_plots):
            row_index = shell_index // num_cols
            col_index = shell_index % num_cols
            if shell_index == num_sub_plots - 1 and num_sub_plots % num_cols:
                ax = plt.subplot(gs[row_index, :])
            else:
                ax = plt.subplot(gs[row_index, col_index])
            length = len(probe_indices[shell_index])
            ax.bar(range(num_data, num_data+length),
                   site_population_list[shell_index],
                   color=self.color_list[shell_index % self.num_colors])
            ax.set_ylabel(f'{shell_index}')
            ax.set_xticks([])
            ax.set_yticks([])
            num_data += length
        figure_name = 'site-wise_occupancy.png'
        figure_path = self.src_path / figure_name
        plt.savefig(str(figure_path))
        return None

    def generate_shell_wise_occupancy(self, num_shells, site_population_list):
        plt.switch_backend('Agg')
        fig = plt.figure()
        plt.title('Shell-wise occupancy')
        ax = fig.add_subplot(111)
        for shell_index in range(num_shells+1):
            mean_value = int(np.mean(site_population_list[shell_index]))
            ax.bar(shell_index, mean_value,
                   color=self.color_list[shell_index % self.num_colors])
            ax.text(shell_index, 1.01 * mean_value, str(mean_value),
                    color='black', horizontalalignment='center')
        ax.set_xlabel('Shell Number')
        ax.set_ylabel('Average shell occupancy')
        xticks_list = [str(index) for index in range(num_shells+1)]
        plt.xticks(range(num_shells+1), xticks_list)
        figure_name = 'shell-wise_occupancy.png'
        figure_path = self.src_path / figure_name
        plt.tight_layout()
        plt.savefig(str(figure_path))
        return None
