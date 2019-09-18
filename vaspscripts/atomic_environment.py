#!/usr/bin/env python

from itertools import combinations

import ase.io.vasp
from ase.atoms import symbols2numbers
import numpy as np

def identify_desired_atom_pair_indices(src_file_path, element_type, desired_pairwise_distance):
    cell = ase.io.vasp.read_vasp(str(src_file_path))
    atomic_number = symbols2numbers(element_type)[0]
    atomic_indices = np.where(cell.numbers == symbols2numbers('O')[0])[0]
    start_index = atomic_indices[0]
    end_index = atomic_indices[-1] + 1
    num_element_type_atoms = len(atomic_indices)
    element_pairwise_distances = np.zeros((num_element_type_atoms, num_element_type_atoms))
    element_type_atomic_indices = np.arange(start_index, end_index)
    for atomic_index in element_type_atomic_indices:
        element_pairwise_distances[atomic_index, :] = cell.get_distances(atomic_index, element_type_atomic_indices, mic=True)

    rounding_digits = len(desired_pairwise_distance.split(".")[1])
    desired_pairwise_distance = float(desired_pairwise_distance)
    desired_pair_internal_indices_temp = np.where(element_pairwise_distances.round(rounding_digits) == desired_pairwise_distance)
    desired_pair_internal_indices = np.hstack((desired_pair_internal_indices_temp[0][:, None], desired_pair_internal_indices_temp[1][:, None]))

    # avoiding duplicate pairs
    desired_pair_internal_indices = desired_pair_internal_indices[desired_pair_internal_indices[:, 1] > desired_pair_internal_indices[:, 0]]
    desired_pair_indices = atomic_indices[desired_pair_internal_indices]
    return (cell, atomic_indices, desired_pair_indices)

def get_well_dispersed_pairs(src_file_path, element_type, desired_pairwise_distance, num_pairs_to_be_selected, input_config=None):
    (cell, _, desired_pair_indices) = identify_desired_atom_pair_indices(src_file_path, element_type, desired_pairwise_distance)
    num_desired_pairs = len(desired_pair_indices)

    # Identify mid-points of desired atomic pairs
    translation_vectors = np.zeros((num_desired_pairs, 3))
    for index in range(num_desired_pairs):
        atom1, atom2 = desired_pair_indices[index]
        translation_vectors[index, :] = cell.get_distance(atom1, atom2, mic=True, vector=True) / 2 

    cell_with_midpoints = cell.__getitem__(desired_pair_indices[:, 0].tolist())
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
        atom_pair_combinations[combination_index] = [desired_pair_indices[pair_index][atom_index] for pair_index in combination for atom_index in range(2)]
        midpoint_pairwise_distances = cell_with_midpoints.__getitem__(list(combination)).get_all_distances(mic=True)
        midpoint_pairwise_distances_array[combination_index, :] = midpoint_pairwise_distances[extract_indices]

    mean_midpoint_pairwise_distances = midpoint_pairwise_distances_array.mean(axis=1)
    std_midpoint_pairwise_distances = midpoint_pairwise_distances_array.std(axis=1)

    # Identify well-dispersed configuration
    min_std_pair_index = np.argmin(std_midpoint_pairwise_distances)
    best_mean_pairwise_distance = mean_midpoint_pairwise_distances[min_std_pair_index]
    min_std_pairwise_distance = std_midpoint_pairwise_distances[min_std_pair_index]
    best_well_dispersed_config = atom_pair_combinations[min_std_pair_index]
    print(f'Mean distance is {best_mean_pairwise_distance:.3f} at the minimum standard deviation of {min_std_pairwise_distance:.3f}')
    best_well_dispersed_fstring = [f'Pair-wise atomic indices (starting with 0) for the well-dispersed config are: ']
    atom_index = 0 
    for pair_index in range(num_pairs_to_be_selected):
        best_well_dispersed_fstring.append(f'({best_well_dispersed_config[atom_index]}, {best_well_dispersed_config[atom_index+1]})')
        atom_index += 2
        if pair_index != num_pairs_to_be_selected-1:
            best_well_dispersed_fstring.append(', ')
        else:
            best_well_dispersed_fstring.append('.')
    print("".join(best_well_dispersed_fstring))

    # Input configuration
    if input_config:
        input_config_index = np.where((atom_pair_combinations == input_config).all(axis=1))[0][0]
        print(f'Mean distance for the input configuration: {mean_midpoint_pairwise_distances[input_config_index]:.3f}')
        print(f'Standard deviation for the input configuraiton: {std_midpoint_pairwise_distances[input_config_index]:.3f}')

    # Save data arrays to dst path
    np.save(src_file_path.parent / 'midpoint_pairwise_distances_array.npy', midpoint_pairwise_distances_array)
    np.save(src_file_path.parent / 'atom_pair_combinations.npy', atom_pair_combinations)
    return None

def get_plane_analysis(src_file_path, element_type, desired_pairwise_distance):
    (cell, atomic_indices, desired_pair_indices) = identify_desired_atom_pair_indices(src_file_path, element_type, desired_pairwise_distance)

    # plane contributions of all points
    cell_lengths = np.linalg.norm(cell.cell, axis=1)
    element_positions = cell.positions[atomic_indices]
    # planes parallel to x/a + y/b = 1
    plane_contributions = element_positions[:, 0] / cell_lengths[0] + element_positions[:, 1] / cell_lengths[1]
    num_planes = 32  # identified from visual observation
    bins = np.linspace(0, 2, num_planes+1)  # 0 through one corner, 1 through diagonal, 2 through diagonally opposite corner
    atoms_sorted_by_plane = np.empty(num_planes, object)
    pair_atoms_in_plane = np.empty(num_planes, object)
    num_atoms_by_plane = np.zeros(num_planes)
    for plane_index in np.arange(num_planes):
        atoms_sorted_by_plane[plane_index] = atomic_indices[np.where((plane_contributions >= bins[plane_index]) & (plane_contributions < bins[plane_index+1]))[0]]
        num_atoms_by_plane[plane_index] = len(atoms_sorted_by_plane[plane_index])
    
        # identify pair atoms in plane
        pair_element_index_list = []
        for atom_index in atoms_sorted_by_plane[plane_index]:
            pair_element_index_tuple = np.where(desired_pair_indices == atom_index)
            if pair_element_index_tuple[0].shape[0]:
                pair_element_index_list.append([pair_element_index_tuple[0][0], pair_element_index_tuple[1][0]])
        pair_element_index_array = np.asarray(pair_element_index_list)
        pair_atoms_in_plane[plane_index] = desired_pair_indices[pair_element_index_array[:, 0], pair_element_index_array[:, 1]]
    return None
