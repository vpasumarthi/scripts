#!/usr/bin/env python

"""
Script to carve out unit cell from super cell
"""

import numpy as np
from shutil import copymode

system_size_input = [2, 2, 2]
num_unitcells = np.prod(system_size_input)
system_name = 'BiVO4_unitcell'
scaling_factor = 1.0
lineno_latticeParameters = range(3, 6)
lineno_start_coordinates = 9

srcFileName = 'POSCAR_Reorded'
dstFileName = 'POSCAR_unitcell'
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
        elif lineNumber in lineno_latticeParameters:
            lattice_vector = np.fromstring(line, sep=' ')
            lattice_vector_scaled = (lattice_vector /
                                     system_size_input[lineNumber - 3])
            line = ' ' + (''.join(['%22.16f' % lattice_vector_scaled[index]
                          for index in range(3)])) + '\n'
        elif lineNumber == 6:
            num_atom_types = len(line.split())
        elif lineNumber == 7:
            atom_numbers_input = map(int, line.split())
            line = ''.join(['%5s' % (atom_numbers_input[index] / num_unitcells)
                            for index in range(num_atom_types)]) + '\n'
        dstFile.write(line)

coordinates_input = np.loadtxt(srcFileName,
                               skiprows=lineno_start_coordinates-1)
num_atoms = sum(atom_numbers_input) / num_unitcells
coordinates_unitcell = (coordinates_input[
                        (coordinates_input[:, 0] > 1. / system_size_input[0]) &
                        (coordinates_input[:, 1] > 1. / system_size_input[1]) &
                        (coordinates_input[:, 2] > 1. / system_size_input[2]),
                        :])
print len(coordinates_unitcell)
coordinates_unitcell *= np.array(system_size_input)[None, :]

np.savetxt(dstFile, coordinates_unitcell, fmt='%20.16f%20.16f%20.16f')
srcFile.close()
dstFile.close()
copymode(srcFileName, dstFileName)
