#!/usr/bin/env python

import os

def generate_symlink(src_path):
    with open(src_path / 'simulation_parameters.yml') as params_file:
        for line in params_file:
            if 'compute_mode' in line:
                compute_mode = line.strip().split()[-1]
            elif 'n_traj' in line:
                n_traj = int(float(line.strip().split()[-1]))
    
    sym_link_file_names = ['Run.py', 'simulation_parameters.yml']
    if compute_mode == 'parallel':
        for traj_index in range(n_traj):
            for file_name in sym_link_file_names:
                src_file_path = src_path / file_name
                dst_file_path = src_path / f'traj{traj_index+1}' / file_name
                import pdb; pdb.set_trace()
                os.symlink(src_file_path, dst_file_path)
