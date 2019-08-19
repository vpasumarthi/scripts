#!/usr/bin/env python

import ase.io.vasp
from ase.atoms import symbols2numbers
import numpy as np

def identify_atom_pairs(src_file_path, element_type):
    cell = ase.io.vasp.read_vasp(str(src_file_path))
    atomic_pairwise_distances = cell.get_all_distances(mic=True)
    atomic_number = symbols2numbers(element_type)[0]
    atomic_indices = np.where(cell.numbers == symbols2numbers('O')[0])[0]
    start_index = atomic_indices[0]
    end_index = atomic_indices[-1] + 1
    element_pairwise_distances = atomic_pairwise_distances[start_index:end_index, start_index:end_index]
    return None

