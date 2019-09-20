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

def get_unique_pathways(compiled_pathways):
    rounding_digits = 3
    compiled_pathways[:, 2:] = np.round(compiled_pathways[:, 2:], rounding_digits)

    # eliminate mirror image symmetry around cartesian axes
    unique_pathways = []
    unique_distances = np.unique(compiled_pathways[:, -1])
    for distance in unique_distances:
        matching_pathway_row_indices = np.where(compiled_pathways[:, -1] == distance)[0]
        absolute_contribution = np.abs(compiled_pathways[matching_pathway_row_indices, 2:])
        unique_rows = np.unique(absolute_contribution, return_index=True, axis=0)[1]
        for row_index in unique_rows:
            unique_pathways.append(compiled_pathways[matching_pathway_row_indices[row_index]])
    unique_pathways = np.asarray(unique_pathways)
    return unique_pathways

def get_plane_analysis(src_file_path, cell_size, element_type,
                       desired_pairwise_distance, num_plane_separation,
                       ref_plane, selection_bounds):
    (cell, atomic_indices, desired_pair_indices) = identify_desired_atom_pair_indices(src_file_path, element_type, desired_pairwise_distance)

    # plane contributions of all points
    cell_lengths = np.linalg.norm(cell.cell, axis=1)
    element_positions = cell.positions[atomic_indices]
    # planes parallel to x/a + y/b = 1
    plane_contributions = element_positions[:, 0] / cell_lengths[0] + element_positions[:, 1] / cell_lengths[1]
    num_planes = 2 * (cell_size[0] + cell_size[1])  # plane intersects a and b axis at half-unit cell length
    bin_edge_shift = plane_contributions.min()
    bins = np.linspace(0, 2, num_planes+1)  + bin_edge_shift # 0 through one corner, 1 through diagonal, 2 through diagonally opposite corner
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
        if len(pair_element_index_list):
            pair_element_index_array = np.asarray(pair_element_index_list)
            pair_atoms_in_plane[plane_index] = desired_pair_indices[pair_element_index_array[:, 0], pair_element_index_array[:, 1]]

    if ref_plane == 'center':
        ref_plane_index = int(num_planes / 2 - 1)
    else:
        ref_plane_index = ref_plane
    upper_plane_index = ref_plane_index + num_plane_separation
    lower_plane_index = ref_plane_index - num_plane_separation

    # select one pair from ref plane
    ref_plane_atom_indices = pair_atoms_in_plane[ref_plane_index]
    ref_plane_atom_positions = cell.positions[ref_plane_atom_indices]
    xmin, xmax, ymin, ymax = selection_bounds
    pair_atoms_within_bounds = ref_plane_atom_indices[(ref_plane_atom_positions[:, 0] / cell_lengths[0] > xmin) & (ref_plane_atom_positions[:, 0] / cell_lengths[0] < xmax) & (ref_plane_atom_positions[:, 1] / cell_lengths[1] > ymin) & (ref_plane_atom_positions[:, 1] / cell_lengths[1] < ymax)]
    ref_pair_atom_of_choice = pair_atoms_within_bounds[0]
    pair_index_of_selected_atom = np.where(desired_pair_indices == ref_pair_atom_of_choice)[0][0]
    atom_indices_of_ref_pair = desired_pair_indices[pair_index_of_selected_atom, :]
    pair_atom1, pair_atom2 = atom_indices_of_ref_pair

    ## find in-plane neighbors
    # find atom indices of neighbors in ref plane
    neighbor_ref_plane_atom_indices = ref_plane_atom_indices.copy()
    for atom_index in atom_indices_of_ref_pair:
        ref_pair_atom_index = np.where(neighbor_ref_plane_atom_indices == atom_index)[0][0]
        neighbor_ref_plane_atom_indices = np.delete(neighbor_ref_plane_atom_indices, ref_pair_atom_index)

    # classify ref pair atom indices between up-the-plane, down-the-plane
    ref_pair_positions = cell.positions[atom_indices_of_ref_pair]
    ref_pair_plane_contributions = ref_pair_positions[:, :2] / np.tile(cell_lengths[:2], (2, 1))
    relative_plane_contributions = ref_pair_plane_contributions[1] / ref_pair_plane_contributions[0]
    if (relative_plane_contributions[0] < 1) & (relative_plane_contributions[1] > 1):
        pair_atom_index_up_the_plane = atom_indices_of_ref_pair[1]
        pair_atom_index_down_the_plane = atom_indices_of_ref_pair[0]
        plane_contribution_up_the_plane = ref_pair_plane_contributions[1]
        plane_contribution_down_the_plane = ref_pair_plane_contributions[0]
    elif (relative_plane_contributions[0] > 1) & (relative_plane_contributions[1] < 1):
        pair_atom_index_up_the_plane = atom_indices_of_ref_pair[0]
        pair_atom_index_down_the_plane = atom_indices_of_ref_pair[1]
        plane_contribution_up_the_plane = ref_pair_plane_contributions[0]
        plane_contribution_down_the_plane = ref_pair_plane_contributions[1]
    else:
        print(f'Pair atoms are not aligned along the direction of plane.')
        exit()

    # compute plane contributions for the remaining atoms
    neighbor_ref_plane_atom_positions = cell.positions[neighbor_ref_plane_atom_indices]
    num_neighbor_ref_plane_atoms = len(neighbor_ref_plane_atom_indices)
    neighbor_ref_plane_contributions = neighbor_ref_plane_atom_positions[:, :2] / np.tile(cell_lengths[:2], (num_neighbor_ref_plane_atoms, 1))

    # divide neighbor atoms into up-the-plane, down-the-plane sections
    neighbor_relative_contributions_up_the_plane = neighbor_ref_plane_contributions / np.tile(plane_contribution_up_the_plane, (num_neighbor_ref_plane_atoms, 1))
    neighbor_relative_contributions_down_the_plane = neighbor_ref_plane_contributions / np.tile(plane_contribution_down_the_plane, (num_neighbor_ref_plane_atoms, 1))
    neighbors_up_the_plane = neighbor_ref_plane_atom_indices[(neighbor_relative_contributions_up_the_plane[:, 0] < 1) &
                                                                 (neighbor_relative_contributions_up_the_plane[:, 1] > 1)]
    neighbors_down_the_plane = neighbor_ref_plane_atom_indices[(neighbor_relative_contributions_down_the_plane[:, 0] > 1) &
                                                                   (neighbor_relative_contributions_down_the_plane[:, 1] < 1)]

    # compute relative distances from the pair atoms
    neighbor_distances_up_the_plane = cell.get_distances(pair_atom_index_up_the_plane, neighbors_up_the_plane, mic=True)
    neighbor_distances_down_the_plane = cell.get_distances(pair_atom_index_down_the_plane, neighbors_down_the_plane, mic=True)
    neighbor_pair_atom_up_the_plane = neighbors_up_the_plane[np.argmin(neighbor_distances_up_the_plane)]
    neighbor_pair_atom_down_the_plane = neighbors_down_the_plane[np.argmin(neighbor_distances_down_the_plane)]
    pair_index_up_the_plane = np.where(desired_pair_indices == neighbor_pair_atom_up_the_plane)[0][0]
    pair_atoms_up_the_plane = desired_pair_indices[pair_index_up_the_plane, :]
    pair_index_down_the_plane = np.where(desired_pair_indices == neighbor_pair_atom_down_the_plane)[0][0]
    pair_atoms_down_the_plane = desired_pair_indices[pair_index_down_the_plane, :]

    # compute inter-pair distances with ref pair atoms: [(0, 0), (0, 1), (1, 0), (1, 1)]
    up_the_plane_distance_vectors_with_ref_pair = np.zeros((4, 3))
    up_the_plane_distance_vectors_with_ref_pair[:2, :] = cell.get_distances(atom_indices_of_ref_pair[0], pair_atoms_up_the_plane, mic=True, vector=True)
    up_the_plane_distance_vectors_with_ref_pair[2:, :] = cell.get_distances(atom_indices_of_ref_pair[1], pair_atoms_up_the_plane, mic=True, vector=True)
    up_the_plane_distances_with_ref_pair = np.linalg.norm(up_the_plane_distance_vectors_with_ref_pair, axis=1)

    down_the_plane_distance_vectors_with_ref_pair = np.zeros((4, 3))
    down_the_plane_distance_vectors_with_ref_pair[:2, :] = cell.get_distances(atom_indices_of_ref_pair[0], pair_atoms_down_the_plane, mic=True, vector=True)
    down_the_plane_distance_vectors_with_ref_pair[2:, :] = cell.get_distances(atom_indices_of_ref_pair[1], pair_atoms_down_the_plane, mic=True, vector=True)
    down_the_plane_distances_with_ref_pair = np.linalg.norm(down_the_plane_distance_vectors_with_ref_pair, axis=1)
    
    neighbor_pair_atoms = np.zeros((num_neighbor_ref_plane_atoms, 2), int)
    neighbor_pair_atoms[:, 0] = neighbor_ref_plane_atom_indices
    for index, neighbor_atom_index in enumerate(neighbor_ref_plane_atom_indices):
        row_index, column_index = np.where(desired_pair_indices == neighbor_atom_index)
        new_column_index = 1 if column_index[0] == 0 else 0
        neighbor_pair_atoms[index, 1] = desired_pair_indices[row_index[0], new_column_index]

    # find unique pathways within cutoff distance
    cutoff_distance = 7.0
    compiled_neighbor_pair_atoms = neighbor_pair_atoms.flatten()
    neighbor_distance_vectors_pair_atom1 = cell.get_distances(pair_atom1, compiled_neighbor_pair_atoms, mic=True, vector=True)
    neighbor_distances_pair_atom1 = np.linalg.norm(neighbor_distance_vectors_pair_atom1, axis=1)
    neighbor_distance_vectors_pair_atom2 = cell.get_distances(pair_atom2, compiled_neighbor_pair_atoms, mic=True, vector=True)
    neighbor_distances_pair_atom2 = np.linalg.norm(neighbor_distance_vectors_pair_atom2, axis=1)

    neighbor_array_indices_pair_atom1 = np.where(neighbor_distances_pair_atom1 < cutoff_distance)[0]
    neighbor_array_indices_pair_atom2 = np.where(neighbor_distances_pair_atom2 < cutoff_distance)[0]
    neighbor_atoms_from_pair_atom1 = compiled_neighbor_pair_atoms[neighbor_array_indices_pair_atom1]
    neighbor_atoms_from_pair_atom2 = compiled_neighbor_pair_atoms[neighbor_array_indices_pair_atom2]
    neighbor_distance_vectors_from_pair_atom1 = neighbor_distance_vectors_pair_atom1[neighbor_array_indices_pair_atom1]
    neighbor_distance_vectors_from_pair_atom2 = neighbor_distance_vectors_pair_atom2[neighbor_array_indices_pair_atom2]
    neighbor_distances_from_pair_atom1 = neighbor_distances_pair_atom1[neighbor_array_indices_pair_atom1]
    neighbor_distances_from_pair_atom2 = neighbor_distances_pair_atom2[neighbor_array_indices_pair_atom2]
    num_in_plane_pathways_pair_atom1 = len(neighbor_atoms_from_pair_atom1)
    num_in_plane_pathways_pair_atom2 = len(neighbor_atoms_from_pair_atom2)
    compiled_in_plane_pathways_pair_atom1 = np.hstack((pair_atom1 * np.ones(num_in_plane_pathways_pair_atom1)[:, None],
                                              neighbor_atoms_from_pair_atom1[:, None],
                                              neighbor_distance_vectors_from_pair_atom1,
                                              neighbor_distances_from_pair_atom1[:, None]))
    compiled_in_plane_pathways_pair_atom2 = np.hstack((pair_atom2 * np.ones(num_in_plane_pathways_pair_atom2)[:, None],
                                              neighbor_atoms_from_pair_atom2[:, None],
                                              neighbor_distance_vectors_from_pair_atom2,
                                              neighbor_distances_from_pair_atom2[:, None]))
    compiled_in_plane_pathways = np.vstack((compiled_in_plane_pathways_pair_atom1, compiled_in_plane_pathways_pair_atom2))
    unique_in_plane_pathways = get_unique_pathways(compiled_in_plane_pathways)

    ## find out-of-plane transport pathways

    # upper-plane
    upper_plane_atom_indices = pair_atoms_in_plane[upper_plane_index]
    upper_plane_atom_indices = upper_plane_atom_indices[(upper_plane_atom_indices != pair_atom1) &
                                                        (upper_plane_atom_indices != pair_atom2)]
    num_upper_plane_atoms = len(upper_plane_atom_indices)
    pair_atoms_in_upper_plane = np.zeros((num_upper_plane_atoms, 2), int)
    pair_atoms_in_upper_plane[:, 0] = upper_plane_atom_indices
    for index, atom_index in enumerate(upper_plane_atom_indices):
        row_index, column_index = np.where(desired_pair_indices == atom_index)
        new_column_index = 1 if column_index[0] == 0 else 0
        pair_atoms_in_upper_plane[index, 1] = desired_pair_indices[row_index[0], new_column_index]
    compiled_upper_plane_pair_atoms = pair_atoms_in_upper_plane.flatten()
    upper_plane_distance_vectors_pair_atom1 = cell.get_distances(pair_atom1, compiled_upper_plane_pair_atoms, mic=True, vector=True)
    upper_plane_distances_pair_atom1 = np.linalg.norm(upper_plane_distance_vectors_pair_atom1, axis=1)
    upper_plane_distance_vectors_pair_atom2 = cell.get_distances(pair_atom2, compiled_upper_plane_pair_atoms, mic=True, vector=True)
    upper_plane_distances_pair_atom2 = np.linalg.norm(upper_plane_distance_vectors_pair_atom2, axis=1)

    upper_plane_neighbor_array_indices_pair_atom1 = np.where(upper_plane_distances_pair_atom1 < cutoff_distance)[0]
    upper_plane_neighbor_array_indices_pair_atom2 = np.where(upper_plane_distances_pair_atom2 < cutoff_distance)[0]
    upper_plane_neighbor_atoms_from_pair_atom1 = compiled_upper_plane_pair_atoms[upper_plane_neighbor_array_indices_pair_atom1]
    upper_plane_neighbor_atoms_from_pair_atom2 = compiled_upper_plane_pair_atoms[upper_plane_neighbor_array_indices_pair_atom2]
    upper_plane_neighbor_distance_vectors_from_pair_atom1 = upper_plane_distance_vectors_pair_atom1[upper_plane_neighbor_array_indices_pair_atom1]
    upper_plane_neighbor_distance_vectors_from_pair_atom2 = upper_plane_distance_vectors_pair_atom2[upper_plane_neighbor_array_indices_pair_atom2]
    upper_plane_neighbor_distances_from_pair_atom1 = upper_plane_distances_pair_atom1[upper_plane_neighbor_array_indices_pair_atom1]
    upper_plane_neighbor_distances_from_pair_atom2 = upper_plane_distances_pair_atom2[upper_plane_neighbor_array_indices_pair_atom2]
    num_upper_plane_pathways_pair_atom1 = len(upper_plane_neighbor_atoms_from_pair_atom1)
    num_upper_plane_pathways_pair_atom2 = len(upper_plane_neighbor_atoms_from_pair_atom2)
    compiled_upper_plane_pathways_pair_atom1 = np.hstack((pair_atom1 * np.ones(num_upper_plane_pathways_pair_atom1)[:, None],
                                              upper_plane_neighbor_atoms_from_pair_atom1[:, None],
                                              upper_plane_neighbor_distance_vectors_from_pair_atom1,
                                              upper_plane_neighbor_distances_from_pair_atom1[:, None]))
    compiled_upper_plane_pathways_pair_atom2 = np.hstack((pair_atom2 * np.ones(num_upper_plane_pathways_pair_atom2)[:, None],
                                              upper_plane_neighbor_atoms_from_pair_atom2[:, None],
                                              upper_plane_neighbor_distance_vectors_from_pair_atom2,
                                              upper_plane_neighbor_distances_from_pair_atom2[:, None]))
    compiled_upper_plane_pathways = np.vstack((compiled_upper_plane_pathways_pair_atom1, compiled_upper_plane_pathways_pair_atom2))

    # lower-plane
    lower_plane_atom_indices = pair_atoms_in_plane[lower_plane_index]
    lower_plane_atom_indices = lower_plane_atom_indices[(lower_plane_atom_indices != pair_atom1) &
                                                        (lower_plane_atom_indices != pair_atom2)]
    num_lower_plane_atoms = len(lower_plane_atom_indices)
    pair_atoms_in_lower_plane = np.zeros((num_lower_plane_atoms, 2), int)
    pair_atoms_in_lower_plane[:, 0] = lower_plane_atom_indices
    for index, atom_index in enumerate(lower_plane_atom_indices):
        row_index, column_index = np.where(desired_pair_indices == atom_index)
        new_column_index = 1 if column_index[0] == 0 else 0
        pair_atoms_in_lower_plane[index, 1] = desired_pair_indices[row_index[0], new_column_index]
    compiled_lower_plane_pair_atoms = pair_atoms_in_lower_plane.flatten()
    lower_plane_distance_vectors_pair_atom1 = cell.get_distances(pair_atom1, compiled_lower_plane_pair_atoms, mic=True, vector=True)
    lower_plane_distances_pair_atom1 = np.linalg.norm(lower_plane_distance_vectors_pair_atom1, axis=1)
    lower_plane_distance_vectors_pair_atom2 = cell.get_distances(pair_atom2, compiled_lower_plane_pair_atoms, mic=True, vector=True)
    lower_plane_distances_pair_atom2 = np.linalg.norm(lower_plane_distance_vectors_pair_atom2, axis=1)

    lower_plane_neighbor_array_indices_pair_atom1 = np.where(lower_plane_distances_pair_atom1 < cutoff_distance)[0]
    lower_plane_neighbor_array_indices_pair_atom2 = np.where(lower_plane_distances_pair_atom2 < cutoff_distance)[0]
    lower_plane_neighbor_atoms_from_pair_atom1 = compiled_lower_plane_pair_atoms[lower_plane_neighbor_array_indices_pair_atom1]
    lower_plane_neighbor_atoms_from_pair_atom2 = compiled_lower_plane_pair_atoms[lower_plane_neighbor_array_indices_pair_atom2]
    lower_plane_neighbor_distance_vectors_from_pair_atom1 = lower_plane_distance_vectors_pair_atom1[lower_plane_neighbor_array_indices_pair_atom1]
    lower_plane_neighbor_distance_vectors_from_pair_atom2 = lower_plane_distance_vectors_pair_atom2[lower_plane_neighbor_array_indices_pair_atom2]
    lower_plane_neighbor_distances_from_pair_atom1 = lower_plane_distances_pair_atom1[lower_plane_neighbor_array_indices_pair_atom1]
    lower_plane_neighbor_distances_from_pair_atom2 = lower_plane_distances_pair_atom2[lower_plane_neighbor_array_indices_pair_atom2]
    num_lower_plane_pathways_pair_atom1 = len(lower_plane_neighbor_atoms_from_pair_atom1)
    num_lower_plane_pathways_pair_atom2 = len(lower_plane_neighbor_atoms_from_pair_atom2)
    compiled_lower_plane_pathways_pair_atom1 = np.hstack((pair_atom1 * np.ones(num_lower_plane_pathways_pair_atom1)[:, None],
                                              lower_plane_neighbor_atoms_from_pair_atom1[:, None],
                                              lower_plane_neighbor_distance_vectors_from_pair_atom1,
                                              lower_plane_neighbor_distances_from_pair_atom1[:, None]))
    compiled_lower_plane_pathways_pair_atom2 = np.hstack((pair_atom2 * np.ones(num_lower_plane_pathways_pair_atom2)[:, None],
                                              lower_plane_neighbor_atoms_from_pair_atom2[:, None],
                                              lower_plane_neighbor_distance_vectors_from_pair_atom2,
                                              lower_plane_neighbor_distances_from_pair_atom2[:, None]))
    compiled_lower_plane_pathways = np.vstack((compiled_lower_plane_pathways_pair_atom1, compiled_lower_plane_pathways_pair_atom2))

    compiled_out_of_plane_pathways = np.vstack((compiled_upper_plane_pathways, compiled_lower_plane_pathways))
    unique_out_of_plane_pathways = get_unique_pathways(compiled_out_of_plane_pathways)
    
    compiled_pathways = np.vstack((unique_in_plane_pathways, unique_out_of_plane_pathways))
    unique_pathways = get_unique_pathways(compiled_pathways)
    num_unique_pathways = len(unique_pathways)
    unique_neighbor_atom_indices = unique_pathways[:, 1].astype(int)
    unique_neighbor_pair_atoms = np.zeros((num_unique_pathways, 2), int)
    unique_neighbor_pair_atoms[:, 0] = unique_neighbor_atom_indices
    for index, atom_index in enumerate(unique_neighbor_atom_indices):
        row_index, column_index = np.where(desired_pair_indices == atom_index)
        new_column_index = 1 if column_index[0] == 0 else 0
        unique_neighbor_pair_atoms[index, 1] = desired_pair_indices[row_index[0], new_column_index]
    return (cell, unique_in_plane_pathways, unique_out_of_plane_pathways)
