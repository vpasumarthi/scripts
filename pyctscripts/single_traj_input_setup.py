#!/usr/bin/env python

# File to setup input files to perform single trajectory simulations

from pathlib import Path
from shutil import copy, move


start_traj = 1
end_traj = 100
copy_file_names = ["Run.py", "slurmscript", "simulation_parameters.yml"]
random_seed_search_term = "random_seed:"
work_dir_depth_search_term = "work_dir_depth:"
work_dir_depth_addon = 2
cwd = Path.cwd()

for traj_number in range(start_traj, end_traj+1):
    # Create directories for each trajectory
    traj_dir_path = cwd / f'{traj_number:03}'
    Path.mkdir(traj_dir_path, parents=True, exist_ok=True)

    # Copy input files to traj directories
    for file_name in copy_file_names:
        copy(cwd / file_name, traj_dir_path)

    # Modify simulation_parameters.yml file
    old_params_file_path = traj_dir_path / "simulation_parameters.yml"
    new_params_file_path = traj_dir_path / "simulation_parameters.yml.new"
    with open(old_params_file_path) as old_params_file, open(new_params_file_path, 'w') as new_params_file:
        for line in old_params_file:
            if random_seed_search_term in line:
                new_params_file.write(f'{random_seed_search_term} {traj_number}\n')
            elif work_dir_depth_search_term in line and "step" not in line:
                old_depth = int(line.strip().split()[1])
                new_depth = old_depth + work_dir_depth_addon
                new_params_file.write(f'{work_dir_depth_search_term} {new_depth}\n')
            else:
                new_params_file.write(line)
    move(new_params_file_path, old_params_file_path)

