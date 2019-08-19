#!/usr/bin/env python

from itertools import combinations

import ase.io.vasp
from ase.atoms import symbols2numbers
import numpy as np


def identify_atom_pairs(src_file_path, element_type, desired_pairwise_distance, num_pairs_to_be_selected):
    cell = ase.io.vasp.read_vasp(str(src_file_path))
    atomic_pairwise_distances = cell.get_all_distances(mic=True)
    atomic_number = symbols2numbers(element_type)[0]
    atomic_indices = np.where(cell.numbers == symbols2numbers('O')[0])[0]
    start_index = atomic_indices[0]
    end_index = atomic_indices[-1] + 1
    element_pairwise_distances = atomic_pairwise_distances[start_index:end_index, start_index:end_index]

    rounding_digits = len(desired_pairwise_distance.split(".")[1])
    desired_pairwise_distance = float(desired_pairwise_distance)
    desired_pairs_temp = np.where(element_pairwise_distances.round(rounding_digits) == desired_pairwise_distance)
    desired_pairs = np.hstack((desired_pairs_temp[0][:, None], desired_pairs_temp[1][:, None]))
    num_desired_pairs = len(desired_pairs)

    # Identify mid-points of desired atomic pairs
    translation_vectors = np.zeros((num_desired_pairs, 3))
    for index in range(num_desired_pairs):
        atom1, atom2 = desired_pairs[index]
        translation_vectors[index, :] = cell.get_distance(atom1, atom2, mic=True, vector=True) / 2 

    cell_with_midpoints = cell.__getitem__(desired_pairs[:, 0].tolist())
    cell_with_midpoints.translate(translation_vectors)
    cell_with_midpoints.wrap()
    midpoint_pair_distances = cell_with_midpoints.get_all_distances(mic=True)

    # Analyze all possible comibinations among the midpoints
    all_combinations = list(combinations(range(num_desired_pairs), num_pairs_to_be_selected)) 
    num_combinations = len(all_combinations)
    atom_pair_combinations = np.zeros((num_combinations, 2 * num_pairs_to_be_selected), int)
    num_distances = len(list(combinations(range(num_pairs_to_be_selected), 2)))
    midpoint_pairwise_distances_array = np.zeros((num_combinations, num_distances)) 
    extract_indices = np.triu_indices(num_pairs_to_be_selected, 1)
    for combination_index, combination in enumerate(all_combinations):
        atom_pair_combinations[combination_index] = [desired_pairs[pair_index][atom_index] for pair_index in combination for atom_index in range(2)]
        midpoint_pairwise_distances = cell_with_midpoints.__getitem__(list(combination)).get_all_distances(mic=True)
        midpoint_pairwise_distance_array[combination_index, :] = midpoint_pairwise_distances[extract_indices] 
    return None

