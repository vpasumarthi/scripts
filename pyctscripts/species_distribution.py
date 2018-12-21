#!/usr/bin/env python

import numpy as np


def get_species_distribution(src_path):
    occupancy = np.load(src_path / 'occupancy.npy')[()]
    time = np.load(src_path / 'time_data.npy')[()]
    return None
