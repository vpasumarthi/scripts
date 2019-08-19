#!/usr/bin/env python

import ase.io.vasp


def identify_atom_pairs(src_file_path):
    cell = ase.io.vasp.read_vasp(str(src_file_path))
    atomic_pairwise_distances = cell.get_all_distances(mic=True)
    return None

