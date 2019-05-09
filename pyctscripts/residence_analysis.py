#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import yaml

from PyCT import constants


class Residence(object):
    def __init__(self, src_path, temp):
        # Load simulation parameters
        sim_param_file_name = 'simulation_parameters.yml'
        sim_param_file_path = src_path / sim_param_file_name
        with open(sim_param_file_path, 'r') as stream:
            try:
                self.sim_params = yaml.load(stream)
            except yaml.YAMLError as exc:
                print(exc)

        self.kBT = constants.KB / constants.EV2J * temp  # kBT in eV

        # doping parameters
        doping_params = self.sim_params['doping']
        self.num_dopants = doping_params['num_dopants']
        self.relative_energies = []
        self.num_shells = []
        self.substitution_element_type_list = []
        substitution_element_type_count = {}
        self.dopant_element_type_list = []
        for element_map in doping_params['doping_element_map']:
            substitution_element_type, dopant_element_type = element_map.split(':')
            self.substitution_element_type_list.append(substitution_element_type)
            self.dopant_element_type_list.append(dopant_element_type)
            if substitution_element_type in substitution_element_type_count:
                substitution_element_type_count[substitution_element_type] += 1
            else:
                substitution_element_type_count[substitution_element_type] = 1
            count_index = substitution_element_type_count[substitution_element_type] - 1
            self.relative_energies.append(self.sim_params['relative_energies']['doping'][substitution_element_type][count_index])
            self.num_shells.append(len(self.relative_energies[-1]) - 1)
        return None

    def traj_shell_wise_residence(self, src_path, traj_index):
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
    
    def shell_wise_residence(self, src_path, n_traj):
        bulk_site_relative_energies = 0.000
        for map_index, relative_energies in enumerate(self.relative_energies):
            if self.num_dopants[map_index]:
                # append relative energies for bulk sites as 0
                relative_energies.append(bulk_site_relative_energies)

                shell_wise_pop_factors = np.exp(- np.asarray(relative_energies) / self.kBT)
                relative_residence_data = np.zeros((n_traj, self.num_shells[map_index] + 2))
                for traj_index in range(n_traj):
                    relative_residence_data[traj_index, :] = self.traj_shell_wise_residence(src_path, traj_index+1)[0]
                shell_wise_num_sites = self.traj_shell_wise_residence(src_path, traj_index+1)[1]
            
                exact_relative_residence = np.multiply(shell_wise_num_sites, shell_wise_pop_factors) / np.dot(shell_wise_num_sites, shell_wise_pop_factors)
                mean_relative_residence_data = np.mean(relative_residence_data, axis=0)
                sem_relative_residence_data = np.std(relative_residence_data, axis=0) / np.sqrt(n_traj)
            
                np.save(src_path / f'exact_relative_residence_{self.dopant_element_type_list[map_index]}.npy', exact_relative_residence)
                np.save(src_path / f'mean_relative_residence_data_{self.dopant_element_type_list[map_index]}.npy', mean_relative_residence_data)
                np.save(src_path / f'sem_relative_residence_data_{self.dopant_element_type_list[map_index]}.npy', sem_relative_residence_data)
        return None
    
    def plot_shell_wise_residence(self, src_path):
        for map_index, dopant_element_type in enumerate(self.dopant_element_type_list):
            if self.num_dopants[map_index]:
                exact_relative_residence = np.load(src_path / f'exact_relative_residence_{dopant_element_type}.npy')
                mean_relative_residence_data = np.load(src_path / f'mean_relative_residence_data_{dopant_element_type}.npy')
                sem_relative_residence_data = np.load(src_path / f'sem_relative_residence_data_{dopant_element_type}.npy')
            
                plt.switch_backend('Agg')
                fig = plt.figure()
                ax = fig.add_subplot(111)
            
                # TODO: Avoid all hard-coding
                shell_index_list = np.arange(len(self.relative_energies[map_index]))
                ax.plot(shell_index_list, mean_relative_residence_data, 'o-',
                         c='#0504aa', mfc='#0504aa', mec='black')
                ax.errorbar(shell_index_list, mean_relative_residence_data,
                             yerr=sem_relative_residence_data, fmt='o', capsize=3,
                             c='#0504aa', mfc='none', mec='none')
                for shell_index in shell_index_list:
                    ax.plot([shell_index - 0.1, shell_index + 0.1], [exact_relative_residence[shell_index]] * 2,
                             '-', c='#d62728')
                ax.set_xlabel('Shell Index')
                ax.set_ylabel('Relative Residence')
                ax.set_title('W10: [0.6596 eV, -0.0168 eV, -0.0154 eV, 0.0000 eV]')
                plt.tight_layout()
                plt.savefig(str(src_path / f'Relative Residence_Shell_wise_{dopant_element_type}.png'))
        return None
