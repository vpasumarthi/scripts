#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

from PyCT.io import read_poscar


def compute_distance(cartesian_coords, system_translational_vector_list,
                     index1, index2):
    center_coord = cartesian_coords[index1]
    neighbor_coord = cartesian_coords[index2]
    neighbor_image_coords = (system_translational_vector_list
                             + neighbor_coord)
    neighbor_image_displacement_vectors = (neighbor_image_coords
                                           - center_coord)
    neighbor_image_displacements = np.linalg.norm(
                            neighbor_image_displacement_vectors, axis=1)
    displacement = np.min(neighbor_image_displacements)
    min_index = np.argmin(neighbor_image_displacements)
    rel_pos_vector = neighbor_image_displacement_vectors[min_index]
    return np.append(rel_pos_vector, displacement)


def penalty_wise_spatial_distribution():
    pbc = [1, 1, 1]
    system_size = np.array([1, 1, 1])  #pseudo
    element_of_interest = 'V'
    dopant_site_number = 36
    system_data = np.loadtxt('W_552.dat')
    poscar_info = read_poscar('POSCAR')
    lattice_matrix = poscar_info['lattice_matrix']
    element_types = poscar_info['element_types']
    num_elements = poscar_info['num_elements']
    coordinate_type = poscar_info['coordinate_type']
    coords = poscar_info['coordinates']
    if coordinate_type == 'Direct':
        fractional_coords = np.copy(coords)
        cartesian_coords = np.dot(fractional_coords, lattice_matrix)
    
    x_range = range(-1, 2) if pbc[0] == 1 else [0]
    y_range = range(-1, 2) if pbc[1] == 1 else [0]
    z_range = range(-1, 2) if pbc[2] == 1 else [0]
    system_translational_vector_list = np.zeros((3**sum(pbc), 3))
    index = 0
    for x_offset in x_range:
        for y_offset in y_range:
            for z_offset in z_range:
                system_translational_vector_list[index] = (
                    np.dot(np.multiply(
                            np.array([x_offset, y_offset, z_offset]),
                            system_size), lattice_matrix))
                index += 1

    site_index_of_interest = element_types.index(element_of_interest)
    num_sites_of_interest = num_elements[site_index_of_interest]
    # (shell index, relpos_x, relpos_y, relpos_z, dist, rel_energy)
    distribution_data = np.zeros((num_sites_of_interest, 6))
    distribution_data[:, 0] = system_data[:, 1]
    dopant_site_index = sum(num_elements[:site_index_of_interest]) + dopant_site_number - 1
    for index, site_index in enumerate(range(num_sites_of_interest)):
        cell_site_index = site_index + sum(num_elements[:site_index_of_interest])
        distribution_data[index, 1:5] = compute_distance(
                            cartesian_coords, system_translational_vector_list,
                            dopant_site_index, cell_site_index)
    distribution_data[:, 5] = system_data[:, 3]

    plt.switch_backend('Agg')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(distribution_data[:, 1],
               distribution_data[:, 2],
               distribution_data[:, 3])
    pnt3d=ax.scatter(distribution_data[:, 1],
                     distribution_data[:, 2],
                     distribution_data[:, 3],
                     c=distribution_data[:, 4])
    cbar=plt.colorbar(pnt3d)
    cbar.set_label("Relative Energy (eV)")
    plt.savefig('Color_by_Distance_552.png', dpi=600)

    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    indices1 = np.where((distribution_data[:, 5] > -0.019) & (distribution_data[:, 5] < 0.005))[0]
    ax1.scatter(distribution_data[indices1, 1],
                distribution_data[indices1, 2],
                distribution_data[indices1, 3], edgecolor='black')
    pnt3d=ax1.scatter(distribution_data[indices1, 1],
                      distribution_data[indices1, 2],
                      distribution_data[indices1, 3],
                      c=distribution_data[indices1, 5],
                      edgecolor='black', linewidth=0.3, cmap='RdYlGn')
    cbar=plt.colorbar(pnt3d)
    cbar.set_label("Relative Energy (eV)")

    ax2 = fig.add_subplot(ax1, projection='3d')
    indices2 = np.where(distribution_data[:, 5] <= -0.019)[0]
    ax2.scatter(distribution_data[indices2, 1],
                distribution_data[indices2, 2],
                distribution_data[indices2, 3],
                c=cm.RdYlGn(0), edgecolor='black', linewidth=0.3)
    
    ax3 = fig.add_subplot(ax1, projection='3d')
    indices3 = np.where(distribution_data[:, 5] >= 0.005)[0]
    ax3.scatter(distribution_data[indices3, 1],
                distribution_data[indices3, 2],
                distribution_data[indices3, 3],
                c=cm.RdYlGn(255), edgecolor='black', linewidth=0.3)

    plt.savefig('Color_by_Relative_Energy_552.png', dpi=1200)
    return None
