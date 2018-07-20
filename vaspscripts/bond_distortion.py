#!/usr/bin/env python

import numpy as np

from PyCT.io import read_poscar, write_poscar


def bond_distortion(src_file_path, polyhedron_stretch, localized_element_type,
                    localized_site_number, neighbor_element_type_list,
                    neighbor_cutoff_list, stretch_percent_list,
                    polyhedron_neighbor_cutoff, polyhedron_neighbor_stretch_percent):
    poscar_info = read_poscar(src_file_path)
    lattice_matrix = poscar_info['lattice_matrix']
    element_types = poscar_info['element_types']
    n_elements = poscar_info['num_elements']
    coordinate_type = poscar_info['coordinate_type']
    unit_cell_coords = poscar_info['coordinates']
    if coordinate_type == 'Direct':
        fractional_coords = unit_cell_coords
    elif coordinate_type == 'Cartesian':
        fractional_coords = np.dot(unit_cell_coords,
                                  np.linalg.inv(lattice_matrix))
    file_format = poscar_info['file_format']

    element_types_consolidated = []
    unique_element_types = set(element_types)
    num_unique_element_types = len(unique_element_types)
    n_elements_consolidated = np.zeros(num_unique_element_types, int)
    unique_element_type_index = -1
    for element_type_index, element_type in enumerate(element_types):
        if element_type not in element_types_consolidated:
            unique_element_type_index += 1
            element_types_consolidated.append(element_type)
        n_elements_consolidated[unique_element_type_index] += n_elements[
                                                            element_type_index]
    localized_element_type_index = element_types_consolidated.index(
                                                        localized_element_type)
    localized_site_coords = (
        fractional_coords[n_elements_consolidated[
                                        :localized_element_type_index].sum()
                          + localized_site_number - 1])

    # generate array of unit cell translational coordinates
    pbc = np.ones(3, int)
    num_cells = 3**sum(pbc)
    x_range = range(-1, 2) if pbc[0] == 1 else [0]
    y_range = range(-1, 2) if pbc[1] == 1 else [0]
    z_range = range(-1, 2) if pbc[2] == 1 else [0]
    cell_translational_coords = np.zeros((num_cells, 3))  # Initialization
    index = 0
    for x_offset in x_range:
        for y_offset in y_range:
            for z_offset in z_range:
                cell_translational_coords[index] = np.array(
                                                [x_offset, y_offset, z_offset])
                index += 1
    localized_site_coords_imageconsolidated = (localized_site_coords
                                               + cell_translational_coords)
    for distort_element_type_index, distort_element_type in enumerate(
                                                neighbor_element_type_list):
        neighbor_element_type_index = element_types_consolidated.index(
                                                        distort_element_type)
        neighbor_site_coords = fractional_coords[
                n_elements_consolidated[:neighbor_element_type_index].sum()
                + range(n_elements_consolidated[neighbor_element_type_index])]
        neighbor_cutoff_dist_limits = [
                        0, neighbor_cutoff_list[distort_element_type_index]]

        if polyhedron_stretch:
            # generate polyhedron vertex site indices
            polyhedron_vertex_site_indices = []
            center_site_coord_list = []
            for neighbor_site_index, neighbor_site_coord in enumerate(
                                                            neighbor_site_coords):
                lattice_directions = (localized_site_coords_imageconsolidated
                                      - neighbor_site_coord)
                min_disp = np.linalg.norm(np.sum(lattice_matrix, axis=0))
                for i_cell in range(num_cells):
                    displacement = np.linalg.norm(
                                np.dot(lattice_directions[i_cell], lattice_matrix))
                    if displacement < min_disp:
                        min_disp = displacement
                        center_site_coords = localized_site_coords_imageconsolidated[i_cell]
                if (neighbor_cutoff_dist_limits[0] < min_disp
                        <= neighbor_cutoff_dist_limits[1]):
                    polyhedron_vertex_site_indices.append(neighbor_site_index)
                    center_site_coord_list.append(center_site_coords)

            # generate neighbor list
            cumulative_vertex_site_neighbor_list = []
            cumulative_vertex_site_neighbor_element_type_list = []
            cumulative_vertex_site_coord_list = []
            vertex_element_type = neighbor_element_type_list[0]
            vertex_element_type_index = element_types_consolidated.index(vertex_element_type)
            vertex_neighbor_element_types_consolidated = element_types_consolidated[:]
            # Avoid neighbors of same element type as vertex site
            vertex_neighbor_element_types_consolidated.remove(vertex_element_type)
            # Generate vertex neighbor site coordinate list
            vertex_neighbor_site_coords = {}
            for neighbor_element_type in vertex_neighbor_element_types_consolidated:
                neighbor_element_type_index = element_types_consolidated.index(neighbor_element_type)
                vertex_neighbor_site_coords[neighbor_element_type] = fractional_coords[
                        n_elements_consolidated[:neighbor_element_type_index].sum()
                        + range(n_elements_consolidated[neighbor_element_type_index])]
                neighbor_cutoff_dist_limits = [0, polyhedron_neighbor_cutoff]
            for vertex_site_index in polyhedron_vertex_site_indices:
                vertex_site_coords = (
                    fractional_coords[n_elements_consolidated[
                                                    :vertex_element_type_index].sum()
                                      + vertex_site_index])
                vertex_site_coords_imageconsolidated = (vertex_site_coords
                                                        + cell_translational_coords)
                vertex_site_neighbor_list = []
                vertex_site_neighbor_element_type_list = []
                vertex_site_coord_list = []
                for neighbor_element_type in vertex_neighbor_element_types_consolidated:
                    for neighbor_site_index, neighbor_site_coord in enumerate(
                                vertex_neighbor_site_coords[neighbor_element_type]):
                        # Avoid localized site from vertex neighbor list
                        if neighbor_element_type != localized_element_type or neighbor_site_index != localized_site_number-1:
                            lattice_directions = (vertex_site_coords_imageconsolidated
                                                  - neighbor_site_coord)
                            min_disp = np.linalg.norm(np.sum(lattice_matrix, axis=0))
                            for i_cell in range(num_cells):
                                displacement = np.linalg.norm(
                                            np.dot(lattice_directions[i_cell], lattice_matrix))
                                if displacement < min_disp:
                                    min_disp = displacement
                                    center_site_coords = vertex_site_coords_imageconsolidated[i_cell]
                            if (neighbor_cutoff_dist_limits[0] < min_disp
                                    <= neighbor_cutoff_dist_limits[1]):
                                vertex_site_neighbor_list.append(neighbor_site_index)
                                vertex_site_neighbor_element_type_list.append(neighbor_element_type)
                                vertex_site_coord_list.append(center_site_coords)
                cumulative_vertex_site_neighbor_list.append(vertex_site_neighbor_list)
                cumulative_vertex_site_neighbor_element_type_list.append(vertex_site_neighbor_element_type_list)
                cumulative_vertex_site_coord_list.append(vertex_site_coord_list)

            # generate distortion
            for vertex_index, vertex_site_index in enumerate(polyhedron_vertex_site_indices):
                neighbor_list = cumulative_vertex_site_neighbor_list[vertex_index]
                num_neighbors = len(neighbor_list)
                neighbor_element_type_list = cumulative_vertex_site_neighbor_element_type_list[vertex_index]
                vertex_site_coord_list = cumulative_vertex_site_coord_list[vertex_index]
                for i_neighbor in range(num_neighbors):
                    neighbor_element_type = neighbor_element_type_list[i_neighbor]
                    neighbor_element_type_index = element_types_consolidated.index(neighbor_element_type)
                    head_start = n_elements_consolidated[:neighbor_element_type_index].sum()
                    neighbor_site_index = neighbor_list[i_neighbor]
                    lattice_direction = (vertex_neighbor_site_coords[neighbor_element_type][neighbor_site_index]
                                        - vertex_site_coord_list[i_neighbor])
                    displacement = np.linalg.norm(np.dot(lattice_direction,
                                                         lattice_matrix))
                    unit_vector = lattice_direction / displacement
                    index = head_start + neighbor_site_index
                    new_coordinate = (
                        vertex_site_coord_list[i_neighbor] + unit_vector
                        * (displacement * (1 + polyhedron_neighbor_stretch_percent / 100)))
                    fractional_coords[index] = new_coordinate
        else:
            # generate neighbor list
            neighbor_list = []
            center_site_coord_list = []
            for neighbor_site_index, neighbor_site_coord in enumerate(
                                                            neighbor_site_coords):
                lattice_directions = (localized_site_coords_imageconsolidated
                                      - neighbor_site_coord)
                min_disp = np.linalg.norm(np.sum(lattice_matrix, axis=0))
                for i_cell in range(num_cells):
                    displacement = np.linalg.norm(
                                np.dot(lattice_directions[i_cell], lattice_matrix))
                    if displacement < min_disp:
                        min_disp = displacement
                        center_site_coords = localized_site_coords_imageconsolidated[i_cell]
                if (neighbor_cutoff_dist_limits[0] < min_disp
                        <= neighbor_cutoff_dist_limits[1]):
                    neighbor_list.append(neighbor_site_index)
                    center_site_coord_list.append(center_site_coords)
    
            # generate distortion
            num_neighbors = len(neighbor_list)
            head_start = n_elements_consolidated[:neighbor_element_type_index].sum()
            for i_neighbor in range(num_neighbors):
                lattice_direction = (neighbor_site_coords[neighbor_list[i_neighbor]]
                                    - center_site_coord_list[i_neighbor])
                displacement = np.linalg.norm(np.dot(lattice_direction,
                                                     lattice_matrix))
                unit_vector = lattice_direction / displacement
                index = head_start + neighbor_list[i_neighbor]
                new_coordinate = (
                    center_site_coord_list[i_neighbor] + unit_vector
                    * (displacement * (1 + stretch_percent_list[distort_element_type_index] / 100)))
                fractional_coords[index] = new_coordinate
    dst_path = src_file_path.parent
    dst_file_name = src_file_path.name + '.out'
    dst_file_path = dst_path / dst_file_name
    write_poscar(src_file_path, dst_file_path, file_format, element_types,
                 n_elements, coordinate_type, fractional_coords)
    return None
