#!/usr/bin/env python


def compute_layer_wise_residence(src_path, system_size, n_traj, gradient_ld,
                         layer_length_ratio):
    """Returns the layer wise residence of charge carriers
    :param src_path:
    :param system_size:
    :param n_traj:
    :param gradient_ld:
    :param layer_length_ratio:
    :return:
    """
    occupancy_data = read_occupancy(src_path, n_traj)
    
    return None