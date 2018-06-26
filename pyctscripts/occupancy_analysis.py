#!/usr/bin/env python

import re

import numpy as np
from scipy.stats import linregress
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

    def generate_occupancy_histogram(self, shell_wise, res_time, site_wise,
                                     barrier_shell_index, n_traj):
        for traj_index in range(n_traj):
            traj_number = traj_index + 1
            (num_shells, probe_indices, site_population_list,
             traj_res_time_pool) = self.read_trajectory_data(traj_number,
                                                             res_time,
                                                             barrier_shell_index)
            if site_wise:
                self.generate_site_wise_occpancy(num_shells,
                                                 probe_indices,
                                                 site_population_list,
                                                 traj_number)
            if shell_wise:
                if traj_number == 1:
                    avg_site_population_list = [
                                [site_population
                                 for site_population in shell_population]
                                for shell_population in site_population_list]
                else:
                    avg_site_population_list = [
                        [site_population + avg_site_population_list[shell_index][site_index]
                         for site_index, site_population in enumerate(shell_population)]
                        for shell_index, shell_population in enumerate(site_population_list)]
                if traj_number == n_traj:
                    avg_site_population_list = [
                        [site_population / n_traj
                         for site_population in shell_population]
                        for shell_population in avg_site_population_list]
                    self.generate_shell_wise_occupancy(num_shells,
                                                       avg_site_population_list,
                                                       n_traj)
            if res_time:
                if traj_number == 1:
                    res_time_pool = []
                res_time_pool.append(traj_res_time_pool)
                if traj_number == n_traj:
                    self.generate_res_time_distribution(res_time_pool, n_traj)
        return None

    def read_trajectory_data(self, traj_number, res_time, barrier_shell_index):
        site_indices_dir_name = 'site_indices_data'
        site_indices_file_name = f'site_indices_{traj_number}.csv'
        site_indices_file_path = self.src_path / site_indices_dir_name / site_indices_file_name
        shell_indices_dict = {}
        site_indices_dict = {}
        with site_indices_file_path.open('r') as site_indices_file:
            for line in site_indices_file:
                split_elements = re.split(',|\n', line)
                shell_index = int(split_elements[3])
                site_index = int(split_elements[0])
                site_indices_dict[site_index] = shell_index
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

        traj_res_time_pool = []
        res_count = 0
        occupancy_dir_name = 'occupancy_data'
        occupancy_file_name = f'occupancy_{traj_number}.dat'
        occupancy_file_path = self.src_path / occupancy_dir_name / occupancy_file_name
        with occupancy_file_path.open('r') as occupancy_file:
            for line in occupancy_file:
                site_index = int(line.split('\n')[0])
                if res_time:
                    shell_index = site_indices_dict[site_index]
                    if shell_index <= barrier_shell_index:
                        res_count += 1
                    else:
                        if res_count != 0:
                            traj_res_time_pool.append(res_count)
                            res_count = 0
                for shell_index in range(num_shells+1):
                    if site_index in probe_indices[shell_index]:
                        list_index = probe_indices[shell_index].index(site_index)
                        site_population_list[shell_index][list_index] += 1
                        break
        return (num_shells, probe_indices, site_population_list,
                traj_res_time_pool)

    def read_occupancy_data(self, traj_number):
        occupancy_dir_name = 'occupancy_data'
        occupancy_file_name = f'occupancy_{traj_number}.dat'
        occupancy_file_path = self.src_path / occupancy_dir_name / occupancy_file_name
        num_kmc_steps = 0
        with occupancy_file_path.open('r') as occupancy_file:
            for index, line in enumerate(occupancy_file):
                if line.strip():
                    if index == 0:
                        num_carriers = len(line.strip().split())
                    num_kmc_steps += 1
        occupancy_data = np.zeros((num_kmc_steps, num_carriers), int)
        with occupancy_file_path.open('r') as occupancy_file:
            for index, line in enumerate(occupancy_file):
                str_site_indices = line.strip().split()
                occupancy_data[index, :] = [int(str_site_index) for str_site_index in str_site_indices]
        return occupancy_data

    def generate_site_wise_occpancy(self, num_shells, probe_indices,
                                    site_population_list, traj_number):
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
        figure_name = f'site-wise_occupancy_{traj_number}.png'
        figure_path = self.src_path / figure_name
        plt.savefig(str(figure_path))
        return None

    def generate_shell_wise_occupancy(self, num_shells, site_population_list,
                                      n_traj):
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
        figure_name = f'avg_shell-wise_occupancy_{n_traj}.png'
        figure_path = self.src_path / figure_name
        plt.tight_layout()
        plt.savefig(str(figure_path))
        return None

    def generate_res_time_distribution(self, res_time_pool, n_traj):
        plt.switch_backend('Agg')
        fig = plt.figure()
        plt.title('Residence length distribution in the doping region')
        ax1 = fig.add_subplot(111)
        cumulative_res_time_pool = [res_time
                                    for traj_res_time_pool in res_time_pool
                                    for res_time in traj_res_time_pool]
        min_res_time = min(cumulative_res_time_pool)
        max_res_time = max(cumulative_res_time_pool)
        num_bins = max_res_time - min_res_time + 1
        bin_edges = np.arange(min_res_time-1, max_res_time+1) + 0.5
        cumulative_pop, _ = np.histogram(cumulative_res_time_pool, bin_edges)
        bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        traj_freq = []
        mean_residence_length_pool = []
        for traj_res_time_pool in res_time_pool:
            i_traj_freq, _ = np.histogram(traj_res_time_pool, bin_edges)
            i_total_freq = sum(i_traj_freq)
            i_prob_freq = i_traj_freq / i_total_freq
            i_mean_residence_length = np.inner(bin_centers, i_prob_freq)
            mean_residence_length_pool.append(i_mean_residence_length)
            traj_freq.append(i_traj_freq)
        traj_freq = np.asarray(traj_freq)
        mean_freq = np.mean(traj_freq, axis=0)
        std_pop = np.std(traj_freq, axis=0)
        plt.plot(bin_centers, mean_freq, 'o-', color='black',
                mfc=self.color_list[1], mec='black')
        plt.errorbar(bin_centers, mean_freq, yerr=std_pop, fmt='o', capsize=3,
                    color='black', mfc='none', mec='none')

        plt.yscale('linear')
        xmin, xmax = plt.xlim()
        ymin, ymax = plt.ylim()
        mean_residence_length = np.mean(mean_residence_length_pool)
        std_mean_residence_length = np.std(mean_residence_length_pool)
        plt.plot([mean_residence_length] * 2, [ymin, ymax], '--',
                 color=self.color_list[2])
        plt.ylim([ymin, ymax])
        floor_index = np.where(bin_centers==np.floor(mean_residence_length))[0]
        err_y_pos = (mean_freq[0] + mean_freq[floor_index]) / 4
        ax1.errorbar(mean_residence_length, err_y_pos,
                    xerr=std_mean_residence_length, fmt='o', capsize=3,
                    color='black', markerfacecolor='none',
                    markeredgecolor='none')
        ax1.text(mean_residence_length + std_mean_residence_length + 0.5,
                0.98 * err_y_pos,
                f'Mean residence length ~ {mean_residence_length:.3f} steps',
                verticalalignment='center')

        ax2 = fig.add_axes([0.5, 0.5, 0.35, 0.35], yscale='log')
        limited_indices = np.where(mean_freq > 1)[0]
        ax2.plot(bin_centers[limited_indices], mean_freq[limited_indices],
                 'o-', color='black', mfc=self.color_list[1], mec='black')
        slope, intercept, r_value, _, _ = \
            linregress(bin_centers[limited_indices], np.log10(mean_freq[limited_indices]))
        ax2.plot(bin_centers[limited_indices],
                 10**(intercept + slope * bin_centers[limited_indices]), '--',
                 color=self.color_list[2])
        sign = '-' if slope < 0 else '+'
        ax2.text(mean_residence_length + std_mean_residence_length + 2,
                 1.2*err_y_pos,
                 'log$_{10}$(y) = ' + f'{intercept:.3f} {sign} {abs(slope):.3f}x',
                 verticalalignment='center')
        ax2.text(mean_residence_length + 12, 0.4*err_y_pos,
                 'r$^2$ = ' + f'{r_value**2:.3f}', verticalalignment='center')

        ax1.set_xlabel('Residence length in number of kmc steps')
        ax1.set_ylabel('Frequency')
        figure_name = f'Residence length distribution_{n_traj}traj.png'
        figure_path = self.src_path / figure_name
        plt.savefig(str(figure_path))
        res_time_data = np.hstack((bin_centers[:, None], mean_freq[:, None],
                                   std_pop[:, None]))
        res_time_data_file_name = f'Residence length distribution_{n_traj}traj.dat'
        np.savetxt(res_time_data_file_name, res_time_data)
        return None

    def get_cell_indices(self, system_size, system_element_index,
                         num_elements_per_unit_cell):
        cell_indices = np.zeros(3, dtype=int)
        unit_cell_element_index = system_element_index % num_elements_per_unit_cell
        n_filled_unit_cells = ((system_element_index - unit_cell_element_index)
                               / num_elements_per_unit_cell)
        for index in range(3):
            cell_indices[index] = (n_filled_unit_cells
                                      / system_size[index+1:].prod())
            n_filled_unit_cells -= (cell_indices[index]
                                    * system_size[index+1:].prod())
        return cell_indices

    def stepwise_res_time(self, system_size, ld, step_length_ratio,
                          num_elements_per_unit_cell, n_traj):
        sum_step_length_ratio = sum(step_length_ratio)
        num_steps = len(step_length_ratio)
        step_system_size_master_list = []
        step_limits = [0]
        old_index = 0
        for step_index in range(num_steps):
            step_system_size = np.copy(system_size)
            step_system_size[ld] *= step_length_ratio[step_index] / sum_step_length_ratio
            step_system_size_master_list.append(step_system_size)
            step_limits.append(old_index + step_system_size[ld])
            old_index += step_system_size[ld]
        step_limits = np.asarray(step_limits)
        step_res_count = np.zeros((n_traj, num_steps), int)
        for traj_number in range(1, n_traj+1):
            traj_step_res_count = np.zeros(num_steps, int)
            occupancy_data = self.read_occupancy_data(traj_number)
            old_kmc_stepwise_step_res_count = np.zeros(num_steps, int)
            if num_steps > 2:
                up_transition_record = np.zeros(num_steps, int)
                down_transition_record = np.zeros(num_steps, int)
            else:
                transition_record = np.zeros(1, int)
            for kmc_step_index, kmc_stepwise_occupancy_data in enumerate(occupancy_data):
                new_kmc_stepwise_step_res_count = np.zeros(num_steps, int)
                for site_index in kmc_stepwise_occupancy_data:
                    cell_indices = self.get_cell_indices(system_size, site_index,
                                                         num_elements_per_unit_cell)
                    site_step_index = sum(step_limits < cell_indices[ld]) - 1
                    new_kmc_stepwise_step_res_count[site_step_index] += 1
                if kmc_step_index != 0 and not np.array_equal(old_kmc_stepwise_step_res_count, new_kmc_stepwise_step_res_count):
                    if num_steps > 2:
                        return
                    else:
                        transition_record += 1
                old_kmc_stepwise_step_res_count = np.copy(new_kmc_stepwise_step_res_count) 
                traj_step_res_count += new_kmc_stepwise_step_res_count
            step_res_count[traj_number-1] = traj_step_res_count
        mean_step_res_count = np.mean(step_res_count, axis=0)
        std_step_res_count = np.std(step_res_count, axis=0)
        stat_width = 10
        stat_decimals = 3
        print(f'Stepwise mean occupancy of electrons is: [' + "".join(f'{val:{stat_width}.{stat_decimals}f}' for val in mean_step_res_count) + ']')
        print(f'Stepwise standard deviation in occupancy of electrons is: [' + "".join(f'{val:{stat_width}.{stat_decimals}f}' for val in std_step_res_count) + ']')
        step_res_time_data_file_name = 'stepwise_residence_occupany.txt'
        step_res_time_data_file_path = self.src_path / step_res_time_data_file_name
        np.savetxt(step_res_time_data_file_path, step_res_count)
        return None
