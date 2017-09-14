#!/usr/bin/env python

"""
Script to reorder input POSCAR by atom name
"""

import numpy as np
from shutil import copymode

system_name = 'BiVO4_2x2x2'
scaling_factor = 1.0
lineno_latticeParameters = range(3, 6)
atom_names = ['Bi', 'V', 'O']
num_atom_types = len(atom_names)
lineno_atom_numbers = 7
lineno_coordinate_type = 8
lineno_start_coordinates = 9
srcFileName = 'POSCAR'
dstFileName = 'POSCAR_Reorded'
srcFile = open(srcFileName, 'r')
open(dstFileName, 'w').close()
dstFile = open(dstFileName, 'a')
for lineIndex, line in enumerate(srcFile):
    lineNumber = lineIndex + 1
    if lineNumber < lineno_start_coordinates:
        if lineNumber == 1:
            line = system_name + '\n'
        elif lineNumber == 2:
            line = ' ' * 3 + '%1.14f\n' % scaling_factor
        elif lineNumber == 6:
            atom_names_input = line.split()
            line = (''.join(['%5s' % atom_names[index]
                    for index in range(num_atom_types)])) + '\n'
        elif lineNumber == 7:
            atom_numbers_input = map(int, line.split())
            atom_numbers = [0] * num_atom_types
            for atom_name_index in range(num_atom_types):
                atom_number_indices = ([index for index, value in
                                        enumerate(atom_names_input)
                                        if
                                        value == atom_names[atom_name_index]])
                for index in atom_number_indices:
                    atom_numbers[atom_name_index] += atom_numbers_input[index]
            line = ''.join(['%5s' % atom_numbers[index]
                            for index in range(num_atom_types)]) + '\n'
        dstFile.write(line)

coordinates_input = np.loadtxt('POSCAR', skiprows=lineno_start_coordinates-1)
num_atoms = sum(atom_numbers)
atom_numbers_input_cumsum = np.cumsum(atom_numbers_input)
coordinates_reordered = np.zeros((num_atoms, 3))
start_index = 0
end_index = 0
for atom_name_index in range(num_atom_types):
    end_index += atom_numbers[atom_name_index]
    atom_number_indices = [index for index, value in
                           enumerate(atom_names_input)
                           if value == atom_names[atom_name_index]]
    row_indices = []
    for index in atom_number_indices:
        if index == 0:
            row_indices.extend(range(atom_numbers_input_cumsum[index]))
        else:
            row_indices.extend(range(atom_numbers_input_cumsum[index-1],
                                     atom_numbers_input_cumsum[index]))
    coordinates_reordered[start_index:end_index, :] = (coordinates_input
                                                       [row_indices, :])
    start_index = end_index
np.savetxt(dstFile, coordinates_reordered, fmt='%20.16f%20.16f%20.16f')
srcFile.close()
dstFile.close()
copymode(srcFileName, dstFileName)
