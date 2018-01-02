#!/usr/bin/env python

from vaspscripts.cluster.hnd_input import hnd_input

input_coordinate_file_name = 'POSCAR01.cart'
dst_file_name = 'POSCAR01.hndin'
dst_element_types = ['V', 'Bi', 'O', 'H']
hnd_input(input_coordinate_file_name, dst_file_name, dst_element_types)
