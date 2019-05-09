#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

class Residence():
    def __init__(self):
        return None

    def traj_shell_wise_residence(self, src_path, traj_index):
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
    
        for index, site_index in enumerate(first_shell_site_indices):
            if index == 0:
                first_shell_occupancy_data = np.where(occupancy[:-1] == site_index)[0]
            else:
                first_shell_occupancy_data = np.hstack((first_shell_occupancy_data, np.where(occupancy[:-1] == site_index)[0]))
    
        for index, site_index in enumerate(second_shell_site_indices):
            if index == 0:
                second_shell_occupancy_data = np.where(occupancy[:-1] == site_index)[0]
            else:
                second_shell_occupancy_data = np.hstack((second_shell_occupancy_data, np.where(occupancy[:-1] == site_index)[0]))
    
        for index, site_index in enumerate(third_shell_site_indices):
            if index == 0:
                third_shell_occupancy_data = np.where(occupancy[:-1] == site_index)[0]
            else:
                third_shell_occupancy_data = np.hstack((third_shell_occupancy_data, np.where(occupancy[:-1] == site_index)[0]))
    
        for index, site_index in enumerate(bulk_site_indices):
            if index == 0:
                bulk_occupancy_data = np.where(occupancy[:-1] == site_index)[0]
            else:
                bulk_occupancy_data = np.hstack((bulk_occupancy_data, np.where(occupancy[:-1] == site_index)[0]))
    
        first_site_time = np.sum(time_step_data[first_shell_occupancy_data])
        second_site_time = np.sum(time_step_data[second_shell_occupancy_data])
        third_site_time = np.sum(time_step_data[third_shell_occupancy_data])
        bulk_site_time = np.sum(time_step_data[bulk_occupancy_data])
    
        site_times = np.array([first_site_time, second_site_time, third_site_time, bulk_site_time])
        traj_relative_residence_data = site_times / np.sum(site_times)
        return (traj_relative_residence_data, shell_wise_num_sites)
    
    def shell_wise_residence(self, src_path, n_traj, kBT, shell_wise_penalties):
        # TODO: parse yaml file to reduce number of inputs
        num_shells = len(shell_wise_penalties)
        shell_wise_pop_factors = np.exp(-shell_wise_penalties / kBT)
        relative_residence_data = np.zeros((n_traj, num_shells))
        for traj_index in range(n_traj):
            relative_residence_data[traj_index, :] = self.traj_shell_wise_residence(src_path, traj_index+1)[0]
        shell_wise_num_sites = self.traj_shell_wise_residence(src_path, traj_index+1)[1]
    
        # TODO: Change abs to exact
        abs_relative_residence = np.multiply(shell_wise_num_sites, shell_wise_pop_factors) / np.dot(shell_wise_num_sites, shell_wise_pop_factors)
        mean_relative_residence_data = np.mean(relative_residence_data, axis=0)
        sem_relative_residence_data = np.std(relative_residence_data, axis=0) / np.sqrt(n_traj)
    
        np.save(src_path / 'abs_relative_residence.npy', abs_relative_residence)
        np.save(src_path / 'mean_relative_residence_data.npy', mean_relative_residence_data)
        np.save(src_path / 'sem_relative_residence_data.npy', sem_relative_residence_data)
        return None
    
    def plot_shell_wise_residence(self, src_path, shell_wise_penalties):
        abs_relative_residence = np.load(src_path / 'abs_relative_residence.npy')
        mean_relative_residence_data = np.load(src_path / 'mean_relative_residence_data.npy')
        sem_relative_residence_data = np.load(src_path / 'sem_relative_residence_data.npy')
    
        plt.switch_backend('Agg')
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
        # TODO: Avoid all hard-coding
        num_shells = len(shell_wise_penalties)
        shell_index_list = np.arange(num_shells)
        ax.plot(shell_index_list, mean_relative_residence_data, 'o-',
                 c='#0504aa', mfc='#0504aa', mec='black')
        ax.errorbar(shell_index_list, mean_relative_residence_data,
                     yerr=sem_relative_residence_data, fmt='o', capsize=3,
                     c='#0504aa', mfc='none', mec='none')
        for shell_index in shell_index_list:
            ax.plot([shell_index - 0.1, shell_index + 0.1], [abs_relative_residence[shell_index]] * 2,
                     '-', c='#d62728')
        ax.set_xlabel('Shell Index')
        ax.set_ylabel('Relative Residence')
        ax.set_title('W10: [0.6596 eV, -0.0168 eV, -0.0154 eV, 0.0000 eV]')
        plt.tight_layout()
        plt.savefig(str(src_path / 'Relative Residence_Shell_wise.png'))
        return None
