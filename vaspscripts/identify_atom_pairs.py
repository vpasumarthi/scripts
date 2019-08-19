#!/usr/bin/env python

import ase.io.vasp
from ase.atoms import symbols2numbers


def identify_atom_pairs(src_file_path, element_type):
    cell = ase.io.vasp.read_vasp(str(src_file_path))
    atomic_pairwise_distances = cell.get_all_distances(mic=True)
    atomic_number = symbols2numbers(element_type)[0]
    return None

