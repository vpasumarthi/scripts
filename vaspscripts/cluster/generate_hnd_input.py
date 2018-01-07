#!/usr/bin/env python

from vaspscripts.cluster.hnd_input import hnd_input

input_coordinate_file_name = 'POSCAR01'
dst_file_name = input_coordinate_file_name + '.hndin'
dst_element_types = ['V', 'O', 'H']
hnd_input(input_coordinate_file_name, dst_file_name, dst_element_types)
