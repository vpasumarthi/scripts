#!/usr/bin/env python

import numpy as np

from PyCT import io


def hnd_input(input_coordinate_file_name, dst_file_name, dst_element_types):
    poscar_info = io.read_poscar(input_coordinate_file_name)
    coordinate_type = poscar_info[-2]
    if coordinate_type == 'Cartesian':
        cartesian_coordinates = poscar_info[-1]
    else:
        lattice_matrix = poscar_info[0]
        fractional_coordinates = poscar_info[-1]
        cartesian_coordinates = np.dot(fractional_coordinates, lattice_matrix)
    element_types = poscar_info[1]
    n_elements = poscar_info[2]
    open(dst_file_name, 'w').close()
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
