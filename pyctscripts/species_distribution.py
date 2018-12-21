#!/usr/bin/env python

import numpy as np


def get_species_distribution(src_path):
    occupancy = np.load(src_path / 'occupancy.npy')[()]
    time = np.load(src_path / 'time_data.npy')[()]
    unique_occupancies = np.unique(occupancy)
    bin_edges = np.hstack((unique_occupancies[0], unique_occupancies[:-1] + np.diff(unique_occupancies) / 2, unique_occupancies[-1]))
    hist = np.histogram(occupancy[:-1], bins=bin_edges, weights=np.diff(time)[:, None])[0]
    return None
