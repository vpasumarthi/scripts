#!/usr/bin/env python

import numpy as np

from PyCT import io


def create_dst_file(dst_file_name):
    open(dst_file_name, 'w').close()
    return None


def basis_order(input_coordinate_file_name, dst_file_name, dst_element_types):
    poscar_info = io.read_poscar(input_coordinate_file_name)
    element_types = poscar_info['element_types']
    n_elements = poscar_info['num_elements']
    dst_file = open(dst_file_name, 'a')

    # Order of atoms for reading basis sets
    headspace = 7
    num_elements_per_line = 12
    dst_file.write(' $BAS\n')
    for element_type in dst_element_types:
        element_index = element_types.index(element_type)
        num_elements = n_elements[element_index]
        while num_elements != 0:
            if num_elements > num_elements_per_line:
                element_list = ', '.join(['%2s' % element_type]
                                         * num_elements_per_line)
                num_elements -= num_elements_per_line
            else:
                element_list = ', '.join(['%2s' % element_type] * num_elements)
                num_elements -= num_elements
            line = (''.join([' ' * headspace, element_list, ',\n']))
            dst_file.write(line)
    dst_file.write(' $END\n\n')
    return None


def geometry(input_coordinate_file_name, dst_file_name, dst_element_types):
    poscar_info = io.read_poscar(input_coordinate_file_name)
    coordinate_type = poscar_info['coordinate_type']
    #import pdb; pdb.set_trace()
    if coordinate_type == 'Cartesian':
        cartesian_coordinates = poscar_info['coordinates']
    else:
        lattice_matrix = poscar_info['lattice_matrix']
        fractional_coordinates = poscar_info['coordinates']
        cartesian_coordinates = np.dot(fractional_coordinates, lattice_matrix)
    element_types = poscar_info['element_types']
    n_elements = poscar_info['num_elements']
    dst_file = open(dst_file_name, 'a')

    # data group to specify geometry
    column_sep_width = 3
    dst_file.write(' $GEO\n')
    for element_type in dst_element_types:
        element_index = element_types.index(element_type)
        num_elements = n_elements[element_index]
        start_index = n_elements[:element_index].sum()
        end_index = start_index + num_elements
        element_type_coordinates = cartesian_coordinates[start_index:end_index]
        for i_element in range(num_elements):
            line = ''.join(['%2s' % element_type,
                            ' ' * column_sep_width,
                            '%9.6f' % element_type_coordinates[i_element][0],
                            ' ' * column_sep_width,
                            '%9.6f' % element_type_coordinates[i_element][1],
                            ' ' * column_sep_width,
                            '%9.6f' % element_type_coordinates[i_element][2],
                            '\n'])
            dst_file.write(line)
    dst_file.write('\n $END\n')
    dst_file.close()
    return None


def frz_orbitals(orbital_indices, dst_file_name):
    dst_file = open(dst_file_name, 'a')
    headspace = 15
    num_orbitals_per_line = 10
    start_line_parameter = 'IFRZ='
    len_parameter = len(start_line_parameter)
    start_line_headspace = headspace - len_parameter
    orbital_index_width = len(str(max(orbital_indices)))
    num_orbitals = num_unprinted_orbitals = len(orbital_indices)
    start_index = 0
    while num_unprinted_orbitals > 0:
        if num_unprinted_orbitals > num_orbitals_per_line:
            end_index = start_index + num_orbitals_per_line
            num_unprinted_orbitals -= num_orbitals_per_line
        else:
            end_index = num_orbitals
            num_unprinted_orbitals = 0
        orbital_index_list = ','.join(
                            [f'{orbital_indices[index]:>{orbital_index_width}}'
                             for index in range(start_index, end_index)])
        if start_index == 0:
            line = (''.join([' ' * start_line_headspace, start_line_parameter,
                             orbital_index_list, ',\n']))
        else:
            line = (''.join([' ' * headspace, orbital_index_list, ',\n']))
        start_index = end_index
        dst_file.write(line)
    dst_file.close()
    return None
