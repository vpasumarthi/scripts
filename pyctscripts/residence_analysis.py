#!/usr/bin/env python

import numpy as np


def traj_shell_wise_residence(traj_index):
    
    return None

def shell_wise_residence(src_path, n_traj, kBT, shell_wise_penalties):
    num_shells = len(shell_wise_penalties)
    relative_residence_data = np.zeros((n_traj, num_shells))
    for traj_index in range(n_traj):
        relative_residence_data[traj_index, :] = traj_shell_wise_residence(src_path, traj_index+1)
    return None
