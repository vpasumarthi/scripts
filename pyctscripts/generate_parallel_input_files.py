#!/usr/bin/env python

import os

def parallel_input_files(src_path):
    with open(src_path / 'simulation_parameters.yml') as params_file:
        for line in params_file:
            if 'compute_mode' in line:
                compute_mode = line.strip().split()[-1]
            elif 'n_traj' in line:
                n_traj = int(float(line.strip().split()[-1]))
    
    sym_link_file_names = ['Run.py', 'simulation_parameters.yml']
    slurm_search_term = "job-name"
    mail_type_search_term = "#SBATCH --mail-type="
    none_term = 'NONE'
    end_term = 'END'
    if compute_mode == 'parallel':
        for traj_index in range(n_traj):
            traj_dir_path = src_path / f'traj{traj_index+1}'
            for file_name in sym_link_file_names:
                src_file_path = src_path / file_name
                dst_file_path = traj_dir_path / file_name
                os.symlink(src_file_path, dst_file_path)

            # generate slurm files
            slurm_file_name = 'slurmscript'
            old_slurm_file_path = src_path / slurm_file_name
            new_slurm_file_path = traj_dir_path / slurm_file_name
            with open(old_slurm_file_path) as old_slurm_file, open(new_slurm_file_path, 'w') as new_slurm_file:
                for line in old_slurm_file:
                    if slurm_search_term in line:
                        new_slurm_file.write(f'{line[:-2]}-traj{traj_index+1}"\n')
                    elif mail_type_search_term in line:
                        if traj_index == n_traj - 1:
                            new_slurm_file.write(f'{mail_type_search_term}{end_term}\n')
                        else:
                            new_slurm_file.write(f'{mail_type_search_term}{none_term}\n')
                    else:
                        new_slurm_file.write(line)
