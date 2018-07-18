#!/usr/bin/env python

import numpy as np

from PyCT.io import read_poscar, write_poscar


def extract_cluster(src_file_path, dst_file_path, site_index_list, bond_limits,
                    terminating_element_type, terminating_bond_distance,
                    oxidation_list, bridge_search_depth, charge_neutral, prec):
    poscar_info = read_poscar(src_file_path)
    lattice_matrix = poscar_info['lattice_matrix']
    element_types = poscar_info['element_types']
    n_elements = poscar_info['num_elements']
    total_elements = poscar_info['total_elements']
    coordinate_type = poscar_info['coordinate_type']
    unit_cell_coords = poscar_info['coordinates']
    if coordinate_type == 'Direct':
        fractional_coords = unit_cell_coords
    elif coordinate_type == 'Cartesian':
        fractional_coords = np.dot(unit_cell_coords,
                                             np.linalg.inv(lattice_matrix))
    file_format = poscar_info['file_format']

    num_sites = len(site_index_list)
    element_types_consolidated = []
    unique_element_types = set(element_types)
    num_unique_element_types = len(unique_element_types)
    n_elements_consolidated = np.zeros(num_unique_element_types, int)
    n_elements_cumulative = n_elements.cumsum()
    unique_element_type_index = -1
    for element_type_index, element_type in enumerate(element_types):
        if element_type not in element_types_consolidated:
            unique_element_type_index += 1
            element_types_consolidated.append(element_type)
        n_elements_consolidated[unique_element_type_index] += n_elements[
                                                            element_type_index]

    # Generate Bonding Neighbor List
    # No PBC implemented here.
    cartesian_coords = np.dot(fractional_coords, lattice_matrix)
    bonding_neighbor_list_indices = np.empty(total_elements, dtype=object)
    element_type_list = []
    for center_site_index, center_coord in enumerate(cartesian_coords):
        element_type_index = np.where(
                            n_elements_cumulative > center_site_index)[0][0]
        center_element_type = element_types[element_type_index]
        element_type_list.append(center_element_type)
        neighbor_list_indices = []
        for neighbor_site_index, neighbor_coord in enumerate(cartesian_coords):
            element_type_index = np.where(
                            n_elements_cumulative > neighbor_site_index)[0][0]
            neighbor_element_type = element_types[element_type_index]
            ref_keys = [':'.join([center_element_type, neighbor_element_type]),
                        ':'.join([neighbor_element_type, center_element_type])]
            if ref_keys[0] in bond_limits:
                bond_limit = bond_limits[ref_keys[0]]
            elif ref_keys[1] in bond_limits:
                bond_limit = bond_limits[ref_keys[1]]
            else:
                bond_limit = 0
            if bond_limit != 0:
                displacement_vector = neighbor_coord - center_coord
                displacement = np.linalg.norm(displacement_vector)
                if (0 < displacement < bond_limit):
                    neighbor_list_indices.append(neighbor_site_index)
        bonding_neighbor_list_indices[
                        center_site_index] = np.asarray(neighbor_list_indices)

    # Generate Cluster
    cluster_element_indices = np.asarray(site_index_list)

    bridge_found = 0
    bridge_depth = 0
    search_index_lists = np.empty(num_sites, dtype=object)
    search_index_lists.fill(np.empty(0, int))
    for site_index in range(num_sites):
        search_index_lists[site_index] = np.append(
                search_index_lists[site_index], site_index_list[site_index])
    while (not bridge_found and bridge_depth < bridge_search_depth):
        for site_index, search_index_list in enumerate(search_index_lists):
            for search_index in search_index_list:
                search_index_lists[site_index] = np.append(
                                search_index_lists[site_index],
                                bonding_neighbor_list_indices[search_index])
            search_index_lists[site_index] = np.unique(search_index_lists[
                                                                site_index])
        bridge_indices = np.intersect1d(search_index_lists[0],
                                        search_index_lists[1])
        bridge_depth += 1
        if len(bridge_indices):
            bridge_found = 1

    # Add neighbor indices up to bridging species
    for site_index in range(num_sites):
        cluster_element_indices = np.append(cluster_element_indices,
                                            search_index_lists[site_index])
    cluster_element_indices = np.unique(cluster_element_indices)

    # Add terminating O sites
    for element_index in cluster_element_indices:
        element_type = element_type_list[element_index]
        if element_type != 'O':
            cluster_element_indices = np.append(
                                cluster_element_indices,
                                bonding_neighbor_list_indices[element_index])
    cluster_element_indices = np.unique(cluster_element_indices)

    # Generate input parameters to write_poscar
    n_elements_cluster = [0] * num_unique_element_types
    coordinates_cluster = []
    for element_index in cluster_element_indices:
        element_type = element_type_list[element_index]
        element_type_index = element_types_consolidated.index(element_type)
        n_elements_cluster[element_type_index] += 1
        coordinates_cluster.append(fractional_coords[element_index])
    non_zero_indices = [index for index in range(num_unique_element_types)
                        if n_elements_cluster[index] != 0]
    n_elements_cluster = [n_elements_cluster[index]
                          for index in non_zero_indices]
    element_types_cluster = [element_types_consolidated[index]
                             for index in non_zero_indices]

    # Generate coordinates of terminating H sites
    h_coordinates_list = []
    h_bond_parent_element_indices = []
    for element_index in cluster_element_indices:
        element_type = element_type_list[element_index]
        if element_type == 'O':
            element_cart_coordinates = cartesian_coords[element_index]
            bonded_indices = bonding_neighbor_list_indices[element_index]
            for bonded_index in bonded_indices:
                if bonded_index not in cluster_element_indices:
                    bonded_element_cart_coordinates = cartesian_coords[
                                                                bonded_index]
                    disp_vector = (bonded_element_cart_coordinates
                                   - element_cart_coordinates)
                    displacement = np.linalg.norm(disp_vector)
                    h_cart_coordinates = (element_cart_coordinates
                                          + (disp_vector / displacement
                                             * terminating_bond_distance))
                    h_fract_coordinates = np.dot(h_cart_coordinates,
                                                 np.linalg.inv(lattice_matrix))
                    h_coordinates_list.append(h_fract_coordinates)
                    h_bond_parent_element_indices.append(element_index)
    h_bond_parent_element_indices = np.asarray(h_bond_parent_element_indices)
    num_h_sites = len(h_coordinates_list)

    # Ensure cluster charge neutrality
    if charge_neutral:
        cluster_charge = 0
        for element_index, element_type in enumerate(element_types_cluster):
            element_charge = oxidation_list[element_type]
            n_elements = n_elements_cluster[element_index]
            cluster_charge += element_charge * n_elements
        cluster_charge += oxidation_list['H'] * num_h_sites
        if cluster_charge % 2 == 0:
            num_pairs = int(cluster_charge / 2)
            center_of_sites = (fractional_coords[site_index_list[0]]
                               + fractional_coords[site_index_list[1]]) / 2
            h_dir_list = []
            h_disp = []
            for h_coordinates in h_coordinates_list:
                dir_vector = h_coordinates - center_of_sites
                disp = np.linalg.norm(np.dot(dir_vector, lattice_matrix))
                h_dir_list.append(dir_vector)
                h_disp.append(disp)
            h_dir_list = np.asarray(h_dir_list)
            h_disp = np.asarray(h_disp)
            sort_indices = h_disp.argsort()
            sorted_h_dir_list = np.round(h_dir_list[sort_indices], prec)
            sorted_h_bond_parent_element_indices = (
                                h_bond_parent_element_indices[sort_indices])
            discard_indices = []
            num_pairs_discarded = 0
            max_index = num_h_sites - 1
            while (num_pairs_discarded != num_pairs and max_index >= 1):
                # Check if the targeted O site has more than 1 proton attached
                if sum(sorted_h_bond_parent_element_indices
                       == sorted_h_bond_parent_element_indices[max_index]) > 1:
                    if (np.array_equal(sorted_h_dir_list[max_index],
                                       -sorted_h_dir_list[max_index - 1])):
                        discard_indices.extend([sort_indices[max_index - 1],
                                                sort_indices[max_index]])
                        sorted_h_bond_parent_element_indices = np.delete(
                                        sorted_h_bond_parent_element_indices,
                                        np.array([max_index - 1, max_index]))
                        max_index -= 2
                        num_pairs_discarded += 1
                    else:
                        max_index -= 1
                else:
                    max_index -= 1
            
            cluster_charge -= 2 * num_pairs_discarded * oxidation_list['H']
            print('Final cluster charge: %d' % cluster_charge)
            h_coordinates_list = [h_coordinates_list[index]
                                  for index in range(num_h_sites)
                                  if index not in discard_indices]
            num_h_sites -= len(discard_indices)
    
    # Add terminating H coordinates
    if num_h_sites:
        element_types_cluster.append(terminating_element_type)
        n_elements_cluster.append(num_h_sites)
        coordinates_cluster.extend(h_coordinates_list)
    print(element_types_cluster)
    print(n_elements_cluster)
    write_poscar(src_file_path, dst_file_path, file_format,
                 element_types_cluster, n_elements_cluster, coordinate_type,
                 coordinates_cluster)
    return None
