#!/usr/bin/env python

import os

def generate_parallel_input_files(src_path):
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
    return None

def generate_slurm_msd_script(src_path):
    old_slurm_file_path = src_path / 'slurmscript'
    new_slurm_file_path = src_path / 'slurmscript_msd'
    job_ids_file_path = src_path / 'job_ids.txt'
    
    insert_line_number = 13
    run_search_term = 'srun Run.py'
    msd_search_term = 'srun MSD.py'
    
    with open(job_ids_file_path, 'r') as job_ids_file:
       job_ids = job_ids_file.readline().strip()[:-1]
    
    with open(old_slurm_file_path, 'r') as old_slurm_file, open(new_slurm_file_path, 'w') as new_slurm_file:
        for line_index, line in enumerate(old_slurm_file):
            if line_index == insert_line_number - 1:
                new_slurm_file.write(f'#SBATCH --dependency=afterany:{job_ids}\n')
                new_slurm_file.write(line)
            elif run_search_term in line:
                new_slurm_file.write(f'srun MSD.py\n')
            elif msd_search_term not in line:
                new_slurm_file.write(line)
    return None
