#!/usr/bin/env python

# File to process output files from single trajectory simulations to a single output

from pathlib import Path


start_traj = 1
end_traj = 100
output_file_names = ["unwrapped_traj.dat", "random_state_info.txt"]
cwd = Path.cwd()

for file_name in output_file_names:
    output_file_path = cwd / file_name
    for traj_number in range(start_traj, end_traj+1):
        traj_dir_path = cwd / f'{traj_number:03}'
        traj_output_path = traj_dir_path / file_name
        with open(traj_output_path) as traj_output_file, open(output_file_path, 'a') as output_file:
            for line in traj_output_file:
                output_file.write(line)

