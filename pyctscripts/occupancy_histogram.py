#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


class Occupancy(object):
    """Class definition to generate occupancy histogram files"""

    def __init__(self, num_probe_shells, color_list, dst_path,
                 shell_index_file_path, occupancy_file_path):
        # Load occupancy parameters
        self.num_probe_shells = num_probe_shells
        self.color_list = color_list
        self.dst_path = dst_path
        self.shell_index_file_path = shell_index_file_path
        self.occupancy_file_path = occupancy_file_path
        self.shell_indices_list = []
        with self.shell_index_file_path.open('r') as shell_index_file:
            for line in shell_index_file:
                self.shell_indices_list.append(
                    [int(index) for index in line.split('\n')[0].split(', ')])
        self.probe_indices = []
        self.site_population_list = []
        for shell_index in range(num_probe_shells+2):
            if shell_index == num_probe_shells + 1:
                self.probe_indices.append(sorted(
                    [item
                     for sublist in self.shell_indices_list[num_probe_shells+1:]
                     for item in sublist]))
            else:
                self.probe_indices.append(sorted(
                    [item for item in self.shell_indices_list[shell_index]]))
            self.site_population_list.append(
                                    [0] * len(self.probe_indices[shell_index]))

        with self.occupancy_file_path.open('r') as occupancy_file:
            for line in occupancy_file:
                site_index = int(line.split('\n')[0])
                for shell_index in range(num_probe_shells+2):
                    if site_index in self.probe_indices[shell_index]:
                        list_index = self.probe_indices[shell_index].index(site_index)
                        self.site_population_list[shell_index][list_index] += 1
                        break


    def generate_site_wise_occpancy(self):
        plt.switch_backend('Agg')
        fig = plt.figure()
        plt.title('Site-wise occupancy')
        plt.axis('off')
        num_sub_plots = self.num_probe_shells + 2
        num_data = 0
        for shell_index in range(num_sub_plots):
            ax = fig.add_subplot(num_sub_plots // 2, 2, shell_index+1)
            length = len(self.probe_indices[shell_index])
            ax.bar(range(num_data, num_data+length),
                   self.site_population_list[shell_index],
                   color=self.color_list[shell_index])
            num_data += length
            ax.set_xlabel(f'Shell {shell_index+1}')
            ax.set_ylabel('Site occupancy')
        ax = fig.add_subplot(111)
        ax.set_xlabel(f'Site Number')
        ax.set_ylabel('Site occupancy')
        xticks_list = [str(index) for index in range(self.num_probe_shells+2)]
        xticks_list[-1] += '+'
        plt.xticks(range(self.num_probe_shells+2), xticks_list)
        figure_name = 'site-wise_occupancy.png'
        figure_path = self.dst_path / figure_name
        plt.tight_layout()
        plt.savefig(str(figure_path))
        return None

    def generate_shell_wise_occupancy(self):
        plt.switch_backend('Agg')
        fig = plt.figure()
        plt.title('Shell-wise occupancy')
        ax = fig.add_subplot(111)
        for shell_index in range(self.num_probe_shells+2):
            mean_value = int(np.mean(self.site_population_list[shell_index]))
            ax.bar(shell_index, mean_value, color=self.color_list[shell_index])
            ax.text(shell_index, 1.01 * mean_value, str(mean_value),
                    color='black', horizontalalignment='center')
        ax.set_xlabel(f'Shell Number')
        ax.set_ylabel('Average shell occupancy')
        xticks_list = [str(index) for index in range(self.num_probe_shells+2)]
        xticks_list[-1] += '+'
        plt.xticks(range(self.num_probe_shells+2), xticks_list)
        figure_name = 'shell-wise_occupancy.png'
        figure_path = self.dst_path / figure_name
        plt.tight_layout()
        plt.savefig(str(figure_path))
        return None
