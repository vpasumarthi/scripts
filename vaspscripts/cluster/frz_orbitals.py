#!/usr/bin/env python


def frz_orbitals(orbital_indices, dst_file_name):
    open(dst_file_name, 'w').close()
    dst_file = open(dst_file_name, 'a')

    headspace = 15 # 10
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
            end_index = num_orbitals - 1
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