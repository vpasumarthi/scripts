#!/usr/bin/env python

# File to setup input files to perform single trajectory simulations

from pathlib import Path


start_traj = 1
end_traj = 100
cwd = Path.cwd()

for traj_number in range(start_traj, end_traj+1):
    traj_dir_path = cwd / f'{traj_number:03}'
    Path.mkdir(traj_dir_path, parents=True, exist_ok=True)

