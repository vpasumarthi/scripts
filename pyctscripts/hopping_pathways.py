#!/usr/bin/env python

from fractions import gcd

import numpy as np

from PyCT.io import read_poscar


class HoppingPathways(object):
    """Class definition to generate charge transfer pathways"""

    def __init__(self, input_coordinate_file_path, class_list):
        self.dst_path = input_coordinate_file_path.parent
        self.class_list = class_list
        poscar_info = read_poscar(input_coordinate_file_path)
        self.lattice_matrix = poscar_info['lattice_matrix']
        self.element_types = poscar_info['element_types']
        self.n_elements_per_unit_cell = poscar_info['num_elements']
        coordinate_type = poscar_info['coordinate_type']
        unit_cell_coords = poscar_info['coordinates']
        if coordinate_type == 'Direct':
            self.fractional_unit_cell_coords = unit_cell_coords
        elif coordinate_type == 'Cartesian':
            self.fractional_unit_cell_coords = np.dot(
                        unit_cell_coords, np.linalg.inv(self.lattice_matrix))

        # define derived parameters
        n_element_types = len(self.element_types)
        self.total_elements_per_unit_cell = self.n_elements_per_unit_cell.sum()
        element_type_index_list = np.repeat(np.arange(n_element_types),
                                            self.n_elements_per_unit_cell)

        # sort element wise coordinates in ascending order of z-coordinate
        start_index = 0
        for element_index in range(n_element_types):
            end_index = start_index + self.n_elements_per_unit_cell[
                                                                element_index]
            element_unit_cell_coords = self.fractional_unit_cell_coords[
                                    element_type_index_list == element_index]
            self.fractional_unit_cell_coords[start_index:end_index] = (
                element_unit_cell_coords[element_unit_cell_coords[:, 2].argsort()])
            start_index = end_index
        return None

    def generate_quantum_indices(self, system_size, system_element_index,
                                 n_elements_per_unit_cell):
        """Returns the quantum indices of the element"""
        quantum_indices = np.zeros(5, dtype=int)  # [0] * 5
        total_elements_per_unit_cell = n_elements_per_unit_cell.sum()
        unitcell_element_index = (system_element_index
                                  % total_elements_per_unit_cell)
        quantum_indices[3] = np.where(np.cumsum(n_elements_per_unit_cell)
                                      >= (unitcell_element_index + 1))[0][0]
        quantum_indices[4] = (
                            unitcell_element_index
                            - n_elements_per_unit_cell[:quantum_indices[3]].sum())
        n_filled_unit_cells = ((system_element_index - unitcell_element_index)
                               / total_elements_per_unit_cell)
        for index in range(3):
            quantum_indices[index] = (n_filled_unit_cells
                                      / system_size[index+1:].prod())
            n_filled_unit_cells -= (quantum_indices[index]
                                    * system_size[index+1:].prod())
        return quantum_indices

    def generate_avoid_element_indices(self, avoid_element_type, num_cells,
                                       center_site_element_type_index):
        if avoid_element_type:
            avoid_element_type_index = self.element_types.index(
                                                            avoid_element_type)
            system_element_index_offset_array = np.repeat(
                np.arange(0, self.total_elements_per_unit_cell * num_cells,
                          self.total_elements_per_unit_cell),
                self.n_elements_per_unit_cell[center_site_element_type_index])
            avoid_element_indices = (
                np.tile(self.n_elements_per_unit_cell[:avoid_element_type_index].sum()
                        + np.arange(0, self.n_elements_per_unit_cell[avoid_element_type_index]),
                        num_cells) + system_element_index_offset_array)
        else:
            avoid_element_indices = []
        return avoid_element_indices

    def generate_bridge_neighbor_list(self, num_center_elements,
                                      center_site_fract_coords,
                                      system_fract_coords, avoid_element_indices,
                                      bridge_cutoff_dist_limits):
        bridge_neighbor_list = np.empty(num_center_elements, dtype=object)
        for center_site_index, center_site_fract_coord in enumerate(
                                                    center_site_fract_coords):
            i_bridge_neighbor_list = []
            for neighbor_site_index, neighbor_site_fract_coord in enumerate(
                                                        system_fract_coords):
                if neighbor_site_index not in avoid_element_indices:
                    lattice_direction = (neighbor_site_fract_coord
                                         - center_site_fract_coord)
                    neighbor_displacement_vector = np.dot(
                            lattice_direction[None, :], self.lattice_matrix)
                    displacement = np.linalg.norm(neighbor_displacement_vector)
                    if (bridge_cutoff_dist_limits[0] < displacement
                            <= bridge_cutoff_dist_limits[1]):
                        i_bridge_neighbor_list.append(neighbor_site_index)
            bridge_neighbor_list[center_site_index] = np.asarray(
                                                        i_bridge_neighbor_list)
        return bridge_neighbor_list
    
    def generate_site_coordinates(self, center_site_element_type_index):
        # generate array of unit cell translational coordinates
        pbc = np.ones(3, int)
        num_cells = 3**sum(pbc)
        x_range = range(-1, 2) if pbc[0] == 1 else [0]
        y_range = range(-1, 2) if pbc[1] == 1 else [0]
        z_range = range(-1, 2) if pbc[2] == 1 else [0]
        system_size = np.array([len(x_range), len(y_range), len(z_range)])
        unitcell_translational_coords = np.zeros((num_cells, 3))  # Initialization
        index = 0
        for x_offset in x_range:
            for y_offset in y_range:
                for z_offset in z_range:
                    unitcell_translational_coords[index] = np.array(
                                                [x_offset, y_offset, z_offset])
                    index += 1

        # extract center site fractional coordinates
        num_center_elements = self.n_elements_per_unit_cell[
                                                center_site_element_type_index]
        center_site_indices = (self.n_elements_per_unit_cell[
                                        :center_site_element_type_index].sum()
                               + np.arange(num_center_elements))
        center_site_fract_coords = self.fractional_unit_cell_coords[
                                                        center_site_indices]
    
        # generate fractional coordinates for neighbor sites
        # and all system elements
        neighbor_site_fract_coords = np.zeros((num_cells * num_center_elements, 3))
        system_fract_coords = np.zeros((
                            num_cells * self.total_elements_per_unit_cell, 3))
        for i_cell in range(num_cells):
            neighbor_site_fract_coords[(i_cell * num_center_elements):(
                                            (i_cell + 1) * num_center_elements)] \
                                        = (center_site_fract_coords
                                           + unitcell_translational_coords[i_cell])
            system_fract_coords[(i_cell * self.total_elements_per_unit_cell):(
                                    (i_cell + 1) * self.total_elements_per_unit_cell)] \
                                        = (self.fractional_unit_cell_coords
                                           + unitcell_translational_coords[i_cell])
        site_coordinate_info = (system_size, num_cells, num_center_elements,
                                center_site_fract_coords,
                                neighbor_site_fract_coords,
                                system_fract_coords)
        return site_coordinate_info

    def generate_raw_pathway_data(
                    self, site_coordinate_info, bridge_neighbor_list,
                    bridge_cutoff_dist_limits, neighbor_cutoff_dist_limits,
                    round_lattice_parameters, desired_coordinate_parameters):
        (system_size, num_cells, num_center_elements, center_site_fract_coords,
         neighbor_site_fract_coords, system_fract_coords) = site_coordinate_info

        if round_lattice_parameters:
            base = round_lattice_parameters['base']
            prec = round_lattice_parameters['prec']

        # initialize class pair list
        if self.class_list:
            center_site_class_list = self.class_list[0]
            neighbor_site_class_list = np.tile(self.class_list[1], num_cells)
            class_pair_list = np.empty(num_center_elements, dtype=object)


        displacement_vector_list = np.empty(num_center_elements, dtype=object)
        lattice_direction_list = np.empty(num_center_elements, dtype=object)
        displacement_list = np.empty(num_center_elements, dtype=object)
        bridge_list = np.empty(num_center_elements, dtype=object)
        num_neighbors = np.zeros(num_center_elements, dtype=int)
        for center_site_index, center_site_fract_coord in enumerate(
                                                    center_site_fract_coords):
            i_displacement_vectors = []
            i_lattice_direction_list = []
            i_displacements = []
            i_bridge_list = []
            if self.class_list:
                i_class_pair_list = []
            for neighbor_site_index, neighbor_site_fract_coord in enumerate(
                                                neighbor_site_fract_coords):
                lattice_direction = (neighbor_site_fract_coord
                                     - center_site_fract_coord)
                neighbor_displacement_vector = np.dot(
                            lattice_direction[None, :], self.lattice_matrix)
                displacement = np.linalg.norm(neighbor_displacement_vector)
                if (neighbor_cutoff_dist_limits[0]
                        < displacement
                        <= neighbor_cutoff_dist_limits[1]):
                    i_displacement_vectors.append(neighbor_displacement_vector)
                    i_displacements.append(displacement)
                    num_neighbors[center_site_index] += 1
                    if round_lattice_parameters:
                        lattice_direction = np.round(
                                (base * np.round((lattice_direction) / base)),
                                prec)
                    i_lattice_direction_list.append(lattice_direction)
    
                    # print fractional coordinates in the desired super cell size
                    if desired_coordinate_parameters:
                        desired_system_size = desired_coordinate_parameters[
                                                        'desired_system_size']
                        dist_list = desired_coordinate_parameters['dist_list']
                        prec = desired_coordinate_parameters['prec']
                        dist = np.round(displacement, prec)
                        if dist in dist_list:
                            print(dist)
                            print('center class:',
                                  center_site_class_list[center_site_index])
                            print('neighbor class:',
                                  neighbor_site_class_list[neighbor_site_index])
                            print('num of bonds:',
                                  len(bridge_neighbor_list[center_site_index]))
                            print('center:',
                                  np.round(np.divide(center_site_fract_coord,
                                                     desired_system_size), 3))
                            print('neighbor:',
                                  np.round(np.divide(neighbor_site_fract_coord,
                                                     desired_system_size), 3))
    
                    # determine class pair list
                    if self.class_list:
                        i_class_pair_list.append(
                                str(center_site_class_list[center_site_index])
                                + ':' + str(neighbor_site_class_list[
                                                        neighbor_site_index]))
    
                    # determine bridging species
                    bridge_site_exists = 0
                    bridge_site_type = ''
                    for i_center_neighbor_se_index in bridge_neighbor_list[
                                                            center_site_index]:
                        i_center_neighbor_fract_coord = system_fract_coords[
                                                    i_center_neighbor_se_index]
                        bridgelattice_direction = (
                                            neighbor_site_fract_coord
                                            - i_center_neighbor_fract_coord)
                        bridgeneighbor_displacement_vector = np.dot(
                                            bridgelattice_direction[None, :],
                                            self.lattice_matrix)
                        bridgedisplacement = np.linalg.norm(
                                            bridgeneighbor_displacement_vector)
                        if (bridge_cutoff_dist_limits[0]
                                < bridgedisplacement
                                <= bridge_cutoff_dist_limits[1]):
                            bridge_site_exists = 1
                            bridge_site_index = i_center_neighbor_se_index
                            bridge_site_quantum_indices = (
                                    self.generate_quantum_indices(
                                                system_size, bridge_site_index,
                                                self.n_elements_per_unit_cell))
                            bridge_site_type += (
                                (', ' if bridge_site_type != '' else '')
                                + self.element_types[bridge_site_quantum_indices[3]])
                    if not bridge_site_exists:
                        bridge_site_type = 'space'
                    i_bridge_list.append(bridge_site_type)
    
            bridge_list[center_site_index] = np.asarray(i_bridge_list)
            displacement_vector_list[center_site_index] = np.asarray(
                                                        i_displacement_vectors)
            lattice_direction_list[center_site_index] = np.asarray(
                                                    i_lattice_direction_list)
            displacement_list[center_site_index] = np.asarray(i_displacements)
            if self.class_list:
                class_pair_list[center_site_index] = np.asarray(
                                                            i_class_pair_list)
        pathway_data = (displacement_list, bridge_list, lattice_direction_list,
                        num_neighbors, class_pair_list, center_site_class_list)
        return pathway_data

    def generate_pathway_list(self, cutoff_dist_key, cutoff,
                              avoid_element_type, precision_parameters,
                              print_parameters, desired_coordinate_parameters):
        """ generate pathway list for the given set of element types"""
        # define input parameters
        neighbor_cutoff = cutoff['neighbor']
        bridge_cutoff = cutoff['bridge']
        round_lattice_parameters = precision_parameters['round_lattice_parameters']
        equivalency_prec = precision_parameters['equivalency']
        pathway_prec = precision_parameters['pathway']
        print_equivalency = print_parameters['equivalency']
        print_pathway_list = print_parameters['pathway_list']
    
        # define derived parameters
        [center_element_type, _] = cutoff_dist_key.split(':')
        center_site_element_type_index = self.element_types.index(center_element_type)
        neighbor_cutoff_dist_limits = [0, neighbor_cutoff]
        bridge_cutoff_dist_limits = [0, bridge_cutoff]
        if round_lattice_parameters:
            base = round_lattice_parameters['base']
    
        site_coordinate_info = self.generate_site_coordinates(
                                                center_site_element_type_index)
        (_, num_cells, num_center_elements, center_site_fract_coords, _,
         system_fract_coords) = site_coordinate_info
    
        # generate list of element indices to avoid during bridge calculations
        avoid_element_indices = self.generate_avoid_element_indices(
                avoid_element_type, num_cells, center_site_element_type_index)

        # generate bridge neighbor list
        bridge_neighbor_list = self.generate_bridge_neighbor_list(
                                num_center_elements, center_site_fract_coords,
                                system_fract_coords, avoid_element_indices,
                                bridge_cutoff_dist_limits)

        # generate pathway data
        pathway_data = self.generate_raw_pathway_data(
                    site_coordinate_info, bridge_neighbor_list,
                    bridge_cutoff_dist_limits, neighbor_cutoff_dist_limits,
                    round_lattice_parameters, desired_coordinate_parameters)
        (displacement_list, bridge_list, lattice_direction_list,
         num_neighbors, class_pair_list, center_site_class_list) = pathway_data

        # determine irreducible form of lattice directions
        sorted_lattice_direction_list = np.empty(num_center_elements,
                                                 dtype=object)
        sorted_displacement_list = np.empty(num_center_elements, dtype=object)
        if self.class_list:
            sorted_class_pair_list = np.empty(num_center_elements,
                                              dtype=object)
        sorted_bridge_list = np.empty(num_center_elements, dtype=object)
        pathway_list = np.empty(num_center_elements, dtype=object)
        for i_center_element_index in range(num_center_elements):
            sorted_displacement_list[i_center_element_index] = (
                    displacement_list[i_center_element_index][
                        displacement_list[i_center_element_index].argsort()])
            sorted_bridge_list[i_center_element_index] = (
                    bridge_list[i_center_element_index][
                        displacement_list[i_center_element_index].argsort()])
            if self.class_list:
                sorted_class_pair_list[i_center_element_index] = (
                    class_pair_list[i_center_element_index][
                        displacement_list[i_center_element_index].argsort()])
            if round_lattice_parameters:
                lattice_direction_list[i_center_element_index] = np.round(
                                lattice_direction_list[i_center_element_index]
                                / base).astype(int)
                center_site_l_d_list = lattice_direction_list[
                                                        i_center_element_index]
                for index in range(num_neighbors[i_center_element_index]):
                    center_site_abs_l_d_list = abs(center_site_l_d_list[index])
                    nz = np.nonzero(center_site_abs_l_d_list)[0]
                    nz_center_site_abs_l_dlist = center_site_abs_l_d_list[nz]
                    if len(nz) == 1:
                        lattice_direction_list[i_center_element_index][
                            index] = (center_site_l_d_list[index]
                                      / center_site_abs_l_d_list[nz])
                    elif len(nz) == 2:
                        lattice_direction_list[i_center_element_index][
                            index] = (center_site_l_d_list[index]
                                      / gcd(nz_center_site_abs_l_dlist[0],
                                            nz_center_site_abs_l_dlist[1]))
                    else:
                        lattice_direction_list[i_center_element_index][
                            index] = (center_site_l_d_list[index]
                                      / gcd(gcd(nz_center_site_abs_l_dlist[0],
                                                nz_center_site_abs_l_dlist[1]),
                                                nz_center_site_abs_l_dlist[2]))
            sorted_lattice_direction_list[i_center_element_index] = (
                lattice_direction_list[i_center_element_index][
                        displacement_list[i_center_element_index].argsort()])
    
            # print equivalency of all center sites with their
            # respective class reference site
            if print_equivalency:
                if self.class_list:
                    ref_index = (np.argmax(center_site_class_list
                                           == center_site_class_list[
                                               i_center_element_index]))
                    ref_index = 0
                print(
                    np.array_equal(np.round(sorted_displacement_list[ref_index],
                                            equivalency_prec),
                                   np.round(sorted_displacement_list[
                                       i_center_element_index],
                                       equivalency_prec)))
    
            # generate center site pathway list
            if self.class_list:
                center_site_pathway_list = np.hstack((
                        np.round(sorted_lattice_direction_list[
                                                        i_center_element_index],
                                 pathway_prec),
                        np.round(sorted_displacement_list[i_center_element_index],
                                 pathway_prec)[:, None],
                        sorted_class_pair_list[i_center_element_index][:, None],
                        sorted_bridge_list[i_center_element_index][:, None]))
            else:
                center_site_pathway_list = np.hstack((
                        np.round(sorted_lattice_direction_list[i_center_element_index],
                                 pathway_prec),
                        np.round(sorted_displacement_list[i_center_element_index],
                                 pathway_prec)[:, None],
                        sorted_bridge_list[i_center_element_index][:, None]))
            pathway_list[i_center_element_index] = center_site_pathway_list
    
            if print_pathway_list:
                np.set_printoptions(suppress=True)
                print(center_site_pathway_list)
    
        pathway_file_name = ('pathway_list_' + cutoff_dist_key.replace(':','-')
                             + '_cutoff=' + str(neighbor_cutoff) + '.npy')
        pathway_file_path = self.dst_path / pathway_file_name
        np.save(pathway_file_path, pathway_list)
        return
