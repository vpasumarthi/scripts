#!/usr/bin/env python

# File to setup input files to perform single trajectory simulations

from pathlib import Path
from shutil import copy


start_traj = 1
end_traj = 100
copy_file_names = ["Run.py", "MSD.py", "slurmscript", "simulation_parameters.yml"]
cwd = Path.cwd()

for traj_number in range(start_traj, end_traj+1):
    traj_dir_path = cwd / f'{traj_number:03}'
    Path.mkdir(traj_dir_path, parents=True, exist_ok=True)

    for file_name in copy_file_names:
        copy(cwd / file_name, traj_dir_path)

