#!/usr/bin/env python

from pathlib import Path
import os

from PyCT.material_preprod import material_preprod

cwd = Path.cwd()
material_preprod(cwd)

with open(cwd / 'simulation_parameters.yml') as params_file:
    for line in params_file:
        if 'compute_mode' in line:
            compute_mode = line.strip().split()[-1]

sym_link_file_names = ['Run.py', 'simulation_parameters.yml']
n_traj = 100
if compute_mode == 'parallel':
    for traj_index in range(n_traj):
        for file_name in sym_link_file_names:
            src_file_path = cwd / file_name
            dst_file_path = cwd / f'traj{traj_index+1}' / file_name
            os.symlink(src_file_path, dst_file_path)
