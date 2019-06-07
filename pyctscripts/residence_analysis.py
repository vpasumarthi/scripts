#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import yaml

from PyCT import constants


class Residence(object):
    def __init__(self, src_path, temp, total_elements_per_unit_cell):
        self.src_path = src_path
        # Load simulation parameters
        sim_param_file_name = 'simulation_parameters.yml'
        sim_param_file_path = self.src_path / sim_param_file_name
        with open(sim_param_file_path, 'r') as stream:
            try:
                self.sim_params = yaml.load(stream)
            except yaml.YAMLError as exc:
                print(exc)

        self.kBT = constants.KB / constants.EV2J * temp  # kBT in eV
        self.total_elements_per_unit_cell = total_elements_per_unit_cell
        self.species_count = self.sim_params['species_count']
        self.num_total_species = np.sum(self.species_count)
        self.system_size = np.asarray(self.sim_params['system_size'])
        self.ndim = len(self.sim_params['pbc'])

        # doping parameters
        self.doping_params = self.sim_params['doping']
        self.num_dopants = self.doping_params['num_dopants']
        self.num_dopant_element_types = len(self.num_dopants)
        bulk_site_relative_energies = 0.0
        self.relative_energies = []
        self.num_shells = []
        self.substitution_element_type_list = []
        substitution_element_type_count = {}
        self.dopant_element_type_list = []
        for element_map in self.doping_params['doping_element_map']:
            substitution_element_type, dopant_element_type = element_map.split(':')
            self.substitution_element_type_list.append(substitution_element_type)
            self.dopant_element_type_list.append(dopant_element_type)
            if substitution_element_type in substitution_element_type_count:
                substitution_element_type_count[substitution_element_type] += 1
            else:
                substitution_element_type_count[substitution_element_type] = 1
            count_index = substitution_element_type_count[substitution_element_type] - 1
            mapping_relative_energies = self.sim_params['relative_energies']['doping'][substitution_element_type][count_index][:]
            mapping_relative_energies.append(bulk_site_relative_energies)
            self.relative_energies.append(mapping_relative_energies)
            self.num_shells.append(len(self.relative_energies[-1]) - 2)
        return None

    def traj_shell_wise_residence(self, traj_index, num_shells):
        site_indices_data = np.load(f'{self.src_path}/traj{traj_index}/site_indices.npy')[()]
        occupancy = np.load(f'{self.src_path}/traj{traj_index}/occupancy.npy')[()]
        time = np.load(f'{self.src_path}/traj{traj_index}/time_data.npy')[()]
        time_step_data = np.tile(np.diff(time)[:, None], self.num_total_species)

        shell_wise_site_count = np.zeros(num_shells+2)
        shell_wise_residence_time = np.zeros(num_shells+2)
        for shell_index in range(num_shells+2):
            if shell_index == num_shells + 1:
                shell_wise_site_indices_data = site_indices_data[site_indices_data[:, 3] > shell_index - 1][:, 0]
            else:
                shell_wise_site_indices_data = site_indices_data[site_indices_data[:, 3] == shell_index][:, 0]
            shell_wise_site_count[shell_index] = len(shell_wise_site_indices_data)
            shell_wise_occupancy_data = np.isin(occupancy[:-1], shell_wise_site_indices_data)
            # NOTE: conslidating shell residence time of all species
            shell_wise_residence_time[shell_index] = np.sum(time_step_data[shell_wise_occupancy_data])

        relative_residence_data = shell_wise_residence_time / np.sum(shell_wise_residence_time)
        return (relative_residence_data, shell_wise_site_count)

    def shell_wise_residence(self, n_traj):
        for map_index, relative_energies in enumerate(self.relative_energies):
            if self.num_dopants[map_index]:
                map_index_relative_energies = relative_energies[:]
                
                shell_wise_pop_factors = np.exp(- np.asarray(map_index_relative_energies) / self.kBT)
                relative_residence_data = np.zeros((n_traj, self.num_shells[map_index] + 2))
                for traj_index in range(n_traj):
                    (relative_residence_data[traj_index, :], shell_wise_num_sites) = self.traj_shell_wise_residence(traj_index+1, self.num_shells[map_index])
            
                exact_relative_residence_data = np.multiply(shell_wise_num_sites, shell_wise_pop_factors) / np.dot(shell_wise_num_sites, shell_wise_pop_factors)
                mean_relative_residence_data = np.mean(relative_residence_data, axis=0)
                sem_relative_residence_data = np.std(relative_residence_data, axis=0) / np.sqrt(n_traj)
            
                np.save(self.src_path / f'shell_exact_relative_residence_data_{self.dopant_element_type_list[map_index]}.npy', exact_relative_residence_data)
                np.save(self.src_path / f'shell_mean_relative_residence_data_{self.dopant_element_type_list[map_index]}.npy', mean_relative_residence_data)
                np.save(self.src_path / f'shell_sem_relative_residence_data_{self.dopant_element_type_list[map_index]}.npy', sem_relative_residence_data)
        return None

    def get_unit_cell_indices(self, site_indices):
        unit_cell_indices_shape = np.append(self.ndim, site_indices.shape)
        unit_cell_indices = np.zeros(unit_cell_indices_shape, int)
        unit_cell_element_indices = site_indices % self.total_elements_per_unit_cell
        total_filled_unit_cells = ((site_indices - unit_cell_element_indices)
                                   // self.total_elements_per_unit_cell)
        for index in range(self.ndim):
            unit_cell_indices[index] = total_filled_unit_cells / self.system_size[index+1:].prod()
            total_filled_unit_cells -= unit_cell_indices[index] * self.system_size[index+1:].prod()
        return unit_cell_indices

    def get_layer_wise_site_indices(self, traj_number, interface, layer_length_ratio, gradient_direction):
        site_indices_data = np.load(f'{self.src_path}/traj{traj_number}/site_indices.npy')[()]
        
        num_layers = len(layer_length_ratio)
        layer_wise_shell_site_indices = [np.empty(shape=(num_layers, map_index_num_shells+2), dtype=object) for map_index_num_shells in self.num_shells]
        bin_edges = np.cumsum(layer_length_ratio) * self.system_size[gradient_direction] / np.sum(layer_length_ratio)
        bin_edges = np.append(np.array([0]), bin_edges)
        for map_index, map_index_num_shells in enumerate(self.num_shells):
            for shell_index in range(map_index_num_shells+2):
                if shell_index == map_index_num_shells + 1:
                    shell_wise_site_indices_data = site_indices_data[(site_indices_data[:, 1] == map_index) & (site_indices_data[:, 3] > shell_index - 1)][:, 0]
                else:
                    shell_wise_site_indices_data = site_indices_data[(site_indices_data[:, 1] == map_index) & (site_indices_data[:, 3] == shell_index)][:, 0]
                if interface == 'flat':
                    unit_cell_indices = self.get_unit_cell_indices(shell_wise_site_indices_data)[gradient_direction]
                    for layer_index in range(num_layers):
                        layer_wise_shell_site_indices[map_index][layer_index, shell_index] = shell_wise_site_indices_data[(unit_cell_indices >= bin_edges[layer_index]) & (unit_cell_indices < bin_edges[layer_index+1])]
                elif interface =='bumpy':
                    if shell_index == 0 or shell_index == map_index_num_shells+1:
                        unit_cell_indices = self.get_unit_cell_indices(shell_wise_site_indices_data)[gradient_direction]
                        for layer_index in range(num_layers):
                            layer_wise_shell_site_indices[map_index][layer_index, shell_index] = shell_wise_site_indices_data[(unit_cell_indices >= bin_edges[layer_index]) & (unit_cell_indices < bin_edges[layer_index+1])]
                    else:
                        dopant_site_shell_index = 0
                        for layer_index in range(num_layers):
                            dopant_site_indices = np.isin(site_indices_data[:, 0], layer_wise_shell_site_indices[map_index][layer_index, dopant_site_shell_index])
                            dopant_element_indices = site_indices_data[dopant_site_indices, 2]
                            layer_wise_shell_site_indices[map_index][layer_index, shell_index] = site_indices_data[:, 0][np.isin(site_indices_data[:, 1], map_index) & np.isin(site_indices_data[:, 2], dopant_element_indices) & np.isin(site_indices_data[:, 3], shell_index)]

        map_index_layer_wise_site_indices = np.empty(shape=(self.num_dopant_element_types, num_layers), dtype=object)
        layer_wise_site_indices = np.empty(num_layers, dtype=object)
        layer_wise_num_sites = np.zeros((self.num_dopant_element_types, num_layers), int)
        for map_index in range(self.num_dopant_element_types):
            for layer_index in range(num_layers):
                map_index_layer_wise_site_indices[map_index, layer_index] = np.hstack(layer_wise_shell_site_indices[map_index][layer_index])
                layer_wise_num_sites[map_index, layer_index] = len(map_index_layer_wise_site_indices[map_index, layer_index])
        for layer_index in range(num_layers):
            layer_wise_site_indices[layer_index] = np.hstack(map_index_layer_wise_site_indices[map_index])
        return (layer_wise_shell_site_indices, layer_wise_site_indices, layer_wise_num_sites, site_indices_data)

    def traj_exact_layer_wise_residence(self, layer_wise_shell_site_indices):
        map_index_layer_based_pop_factors = []
        for map_index, map_index_relative_energies in enumerate(self.relative_energies):
            (num_layers, map_index_num_shells) = layer_wise_shell_site_indices[map_index].shape
            map_index_layer_based_pop_factors.append(np.empty(num_layers))
            shell_wise_pop_factors = np.exp(- np.asarray(map_index_relative_energies) / self.kBT)
            for layer_index in range(num_layers):
                layer_shell_wise_num_sites = np.zeros(map_index_num_shells)
                for shell_index in range(map_index_num_shells):
                    layer_shell_wise_num_sites[shell_index] = len(layer_wise_shell_site_indices[map_index][layer_index][shell_index])
                map_index_layer_based_pop_factors[map_index][layer_index] = np.dot(shell_wise_pop_factors, layer_shell_wise_num_sites)
        layer_based_pop_factors = np.asarray(map_index_layer_based_pop_factors).sum(axis=0)
        exact_relative_residence_data = layer_based_pop_factors / np.sum(layer_based_pop_factors)
        return exact_relative_residence_data

    def traj_layer_wise_residence(self, traj_number, site_indices_data, layer_wise_site_indices):
        occupancy = np.load(f'{self.src_path}/traj{traj_number}/occupancy.npy')[()]
        time = np.load(f'{self.src_path}/traj{traj_number}/time_data.npy')[()]
        time_step_data = np.tile(np.diff(time)[:, None], self.num_total_species)

        num_layers = len(layer_wise_site_indices)
        layer_wise_residence = np.zeros(num_layers)
        MINBINS = site_indices_data[-1, 0] + 1
        occupant_site_wise_residence = np.bincount(occupancy[:-1].reshape(-1), time_step_data.reshape(-1), MINBINS)
        for layer_index in range(num_layers):
            layer_wise_residence[layer_index] = occupant_site_wise_residence[np.unique(layer_wise_site_indices[layer_index])].sum()
        traj_relative_residence_data = layer_wise_residence / np.sum(layer_wise_residence)
        return traj_relative_residence_data

    def layer_wise_residence(self, n_traj, interface, return_num_accessible_sites):
        # NOTE: Assuming identical gradient direction and layer_length_ratio for all existing dopant element types
        sample_existing_map_index = (np.asarray(self.num_dopants) > 0).tolist().index(True)
        layer_length_ratio = self.doping_params['gradient'][sample_existing_map_index]['step_length_ratio']
        num_layers = len(layer_length_ratio)
        gradient_direction = self.doping_params['gradient'][sample_existing_map_index]['ld']

        relative_residence_data = np.zeros((n_traj, num_layers))
        exact_relative_residence_data = np.zeros((n_traj, num_layers))
        layer_wise_num_sites_data = np.zeros((n_traj, self.num_dopant_element_types, num_layers), int)
        for traj_index in range(n_traj):
            (layer_wise_shell_site_indices, layer_wise_site_indices, layer_wise_num_sites_data[traj_index], site_indices_data) = self.get_layer_wise_site_indices(traj_index+1, interface, layer_length_ratio, gradient_direction)
            exact_relative_residence_data[traj_index, :] = self.traj_exact_layer_wise_residence(layer_wise_shell_site_indices)
            relative_residence_data[traj_index, :] = self.traj_layer_wise_residence(traj_index+1, site_indices_data, layer_wise_site_indices)

        mean_relative_residence_data = np.mean(relative_residence_data, axis=0)
        sem_relative_residence_data = np.std(relative_residence_data, axis=0) / np.sqrt(n_traj)
        mean_exact_relative_residence_data = np.mean(exact_relative_residence_data, axis=0)
        sem_exact_relative_residence_data = np.std(exact_relative_residence_data, axis=0) / np.sqrt(n_traj)

        percent_deviation = np.divide((exact_relative_residence_data - relative_residence_data), exact_relative_residence_data) * 100
        mean_percent_deviation = np.mean(percent_deviation, axis=0)
        sem_percent_deviation = np.std(percent_deviation, axis=0) / np.sqrt(n_traj)

        mean_layer_wise_num_sites_data = np.mean(layer_wise_num_sites_data, axis=0)
        sem_layer_wise_num_sites_data = np.std(layer_wise_num_sites_data, axis=0) / np.sqrt(n_traj)

        np.save(self.src_path / f'layer_{interface}_mean_relative_residence_data.npy', mean_relative_residence_data)
        np.save(self.src_path / f'layer_{interface}_sem_relative_residence_data.npy', sem_relative_residence_data)
        np.save(self.src_path / f'layer_{interface}_mean_exact_relative_residence_data.npy', mean_exact_relative_residence_data)
        np.save(self.src_path / f'layer_{interface}_sem_exact_relative_residence_data.npy', sem_exact_relative_residence_data)
        np.save(self.src_path / f'layer_{interface}_mean_percent_deviation.npy', mean_percent_deviation)
        np.save(self.src_path / f'layer_{interface}_sem_percent_deviation.npy', sem_percent_deviation)
        if return_num_accessible_sites:
            np.save(self.src_path / f'layer_{interface}_mean_layer_wise_num_sites_data.npy', mean_layer_wise_num_sites_data)
            np.save(self.src_path / f'layer_{interface}_sem_layer_wise_num_sites_data.npy', sem_layer_wise_num_sites_data)
        return None

    def plot_shell_wise_residence(self, show_exact):
        # Plot specifications
        figure_dpi = 600

        # Font specifications
        font_family = 'sans-serif'
        font_name = 'Calibri'
        plt.rcParams['font.family'] = font_family
        plt.rcParams['font.sans-serif'] = [font_name]
        title_size = 18
        font_size = 16
        label_size = 12

        for map_index, dopant_element_type in enumerate(self.dopant_element_type_list):
            if self.num_dopants[map_index]:
                map_index_relative_energies = self.relative_energies[map_index][:]
                num_shells = len(map_index_relative_energies) - 2

                # show exact relative residence values for single species
                if show_exact:
                    exact_relative_residence = np.load(self.src_path / f'shell_exact_relative_residence_data_{self.dopant_element_type_list[map_index]}.npy')
                mean_relative_residence_data = np.load(self.src_path / f'shell_mean_relative_residence_data_{self.dopant_element_type_list[map_index]}.npy')
                sem_relative_residence_data = np.load(self.src_path / f'shell_sem_relative_residence_data_{self.dopant_element_type_list[map_index]}.npy')
            
                plt.switch_backend('Agg')
                fig = plt.figure()
                ax = fig.add_subplot(111)
            
                shell_index_list = np.arange(len(self.relative_energies[map_index]))
                ax.plot(shell_index_list, mean_relative_residence_data, 'o-',
                         c='#0504aa', mfc='#0504aa', mec='black', label='simulation')
                ax.errorbar(shell_index_list, mean_relative_residence_data,
                             yerr=sem_relative_residence_data, fmt='o', capsize=3,
                             c='#0504aa', mfc='none', mec='none')
                if show_exact:
                    for shell_index in shell_index_list:
                        if shell_index == shell_index_list[0]:
                            ax.plot([shell_index - 0.1, shell_index + 0.1], [exact_relative_residence[shell_index]] * 2,
                                     '-', c='#d62728', label='exact')
                        else:
                            ax.plot([shell_index - 0.1, shell_index + 0.1], [exact_relative_residence[shell_index]] * 2,
                                     '-', c='#d62728')

                x_ticks = np.arange(num_shells+2)
                x_tick_labels = [str(tick) for tick in x_ticks]
                plt.xticks(x_ticks, x_tick_labels, fontsize=label_size)
                plt.yticks(fontsize=label_size)

                ax.legend(fontsize=label_size)
                ax.set_xlabel('Shell Index', fontsize=font_size)
                ax.set_ylabel('Relative Residence', fontsize=font_size)
                ax.set_title(f'{dopant_element_type}{self.num_dopants[map_index]:02d}: {num_shells}shells; e{self.species_count[0]}h{self.species_count[1]}', fontsize=title_size)
                plt.tight_layout()
                plt.savefig(str(self.src_path / f'Relative Residence_Shell_wise_{dopant_element_type}.png'), dpi=figure_dpi)
        return None

    def plot_layer_wise_residence(self, interface, show_exact, plot_num_accessible_sites):
        # Plot specifications
        figure_dpi = 600

        # Font specifications
        font_family = 'sans-serif'
        font_name = 'Calibri'
        plt.rcParams['font.family'] = font_family
        plt.rcParams['font.sans-serif'] = [font_name]
        title_size = 18
        font_size = 16
        label_size = 12

        for map_index, dopant_element_type in enumerate(self.dopant_element_type_list):
            if self.num_dopants[map_index]:
                map_index_relative_energies = self.relative_energies[map_index][:]
                num_shells = len(map_index_relative_energies) - 2

                mean_relative_residence_data = np.load(self.src_path / f'layer_{interface}_mean_relative_residence_data.npy')
                sem_relative_residence_data = np.load(self.src_path / f'layer_{interface}_sem_relative_residence_data.npy')
        
                if plot_num_accessible_sites:
                    mean_layer_wise_num_sites_data = np.load(self.src_path / f'layer_{interface}_mean_layer_wise_num_sites_data.npy')
                    sem_layer_wise_num_sites_data = np.load(self.src_path / f'layer_{interface}_sem_layer_wise_num_sites_data.npy')
        
                # show exact relative residence values for single species
                if show_exact:
                    mean_exact_relative_residence_data = np.load(self.src_path / f'layer_{interface}_mean_exact_relative_residence_data.npy')
                    sem_exact_relative_residence_data = np.load(self.src_path / f'layer_{interface}_sem_exact_relative_residence_data.npy')
                    mean_percent_deviation = np.load(self.src_path / f'layer_{interface}_mean_percent_deviation.npy')
                    sem_percent_deviation = np.load(self.src_path / f'layer_{interface}_sem_percent_deviation.npy')

                plt.switch_backend('Agg')
                fig1 = plt.figure()
                ax1 = fig1.add_subplot(111)

                layer_length_ratio = self.doping_params['gradient'][map_index]['step_length_ratio']
                num_layers = len(layer_length_ratio)
                layer_index_list = np.arange(num_layers)
                ax1.plot(layer_index_list, mean_relative_residence_data, 'o-',
                         c='#0504aa', mfc='#0504aa', mec='black', label='simulation')
                ax1.errorbar(layer_index_list, mean_relative_residence_data,
                             yerr=sem_relative_residence_data, fmt='o', capsize=3,
                             c='#0504aa', mfc='none', mec='none')
                if show_exact:
                    ax1.plot(layer_index_list, mean_exact_relative_residence_data, 'o-',
                             c='#d62728', mfc='#d62728', mec='black', label='prediction (1 species)')
                    ax1.errorbar(layer_index_list, mean_exact_relative_residence_data,
                                 yerr=sem_exact_relative_residence_data, fmt='o', capsize=3,
                                 c='#d62728', mfc='none', mec='none')

                x_ticks = np.arange(num_layers)
                x_tick_labels = [str(tick) for tick in x_ticks]
                plt.xticks(x_ticks, x_tick_labels, fontsize=label_size)
                plt.yticks(fontsize=label_size)

                ax1.legend(fontsize=label_size)
                ax1.set_xlabel('Layer Index', fontsize=font_size)
                ax1.set_ylabel('Relative Residence', fontsize=font_size)
                ax1.set_title(f'{dopant_element_type}{self.num_dopants[map_index]:02d}: {num_shells}shells; e{self.species_count[0]}h{self.species_count[1]} in L{num_layers} ({interface})', fontsize=title_size)
                plt.tight_layout()
                plt.savefig(str(self.src_path / f'Relative Residence_Layer_wise_{interface}_{dopant_element_type}.png'), dpi=figure_dpi)

                if plot_num_accessible_sites:
                    # Layer-wise number of sites
                    fig2 = plt.figure()
                    ax2 = fig2.add_subplot(111)
    
                    ax2.plot(layer_index_list, mean_layer_wise_num_sites_data, 'o-',
                             c='#0504aa', mfc='#0504aa', mec='black', label='simulation')
                    ax2.errorbar(layer_index_list, mean_layer_wise_num_sites_data,
                                 yerr=sem_layer_wise_num_sites_data, fmt='o', capsize=3,
                                 c='#0504aa', mfc='none', mec='none')
    
                    x_ticks = np.arange(num_layers)
                    x_tick_labels = [str(tick) for tick in x_ticks]
                    plt.xticks(x_ticks, x_tick_labels, fontsize=label_size)
                    plt.yticks(fontsize=label_size)
    
                    ax2.legend(fontsize=label_size)
                    ax2.set_xlabel('Layer Index', fontsize=font_size)
                    ax2.set_ylabel('Number of accessible sites', fontsize=font_size)
                    ax2.set_title(f'{dopant_element_type}{self.num_dopants[map_index]:02d}: {num_shells}shells; e{self.species_count[0]}h{self.species_count[1]} in L{num_layers} ({interface})', fontsize=title_size)
                    plt.tight_layout()
                    plt.savefig(str(self.src_path / f'Layer_wise_Number_of_sites_{interface}_{dopant_element_type}.png'), dpi=figure_dpi)

                if show_exact:
                    fig3 = plt.figure()
                    ax3 = fig3.add_subplot(111)
                    ax3.plot(layer_index_list, mean_percent_deviation, 'o-',
                             c='#0504aa', mfc='#0504aa', mec='black')
                    ax3.errorbar(layer_index_list, mean_percent_deviation,
                                 yerr=sem_percent_deviation, fmt='o', capsize=3,
                                 c='#0504aa', mfc='none', mec='none')

                    x_ticks = np.arange(num_layers)
                    x_tick_labels = [str(tick) for tick in x_ticks]
                    plt.xticks(x_ticks, x_tick_labels, fontsize=label_size)
                    plt.yticks(fontsize=label_size)

                    ax3.set_xlabel('Layer Index', fontsize=font_size)
                    ax3.set_ylabel('Relative Residence Deviation (%)', fontsize=font_size)
                    ax3.set_title(f'{dopant_element_type}{self.num_dopants[map_index]:02d}: {num_shells}shells; e{self.species_count[0]}h{self.species_count[1]} in L{num_layers} ({interface})', fontsize=title_size)
                    plt.tight_layout()
                    plt.savefig(str(self.src_path / f'Relative Residence Deviation_Layer_wise_{interface}_{dopant_element_type}.png'), dpi=figure_dpi)
        return None
