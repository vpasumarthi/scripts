#!/usr/bin/env python

# File to generate input files for charge localization jobs

from pathlib import Path
from shutil import copy, copymode, move

import numpy as np


site_indices = np.loadtxt('site_indices.dat', dtype=int)
cwd = Path.cwd()
partial_run_index = 1
copy_file_names = ["POSCAR", "POTCAR", "INCAR", "KPOINTS", "slurmscript", "generate_bond_distortion.py"]
search_term = 'localized_site_number = '

for site_index in site_indices:
    # Generate input files
    site_dir_path = cwd / f'W36_ox6_V{site_index:02}_ox4'
    work_dir_path = site_dir_path / f'partial_run{partial_run_index:02}'
    Path.mkdir(work_dir_path, parents=True, exist_ok=True)

    for file_name in copy_file_names:
        copy(cwd / file_name, work_dir_path)
    old_file_path = work_dir_path / 'generate_bond_distortion.py'
    new_file_path = work_dir_path / 'generate_bond_distortion.py.new'
    with open(old_file_path) as old_file, open(new_file_path, 'w') as new_file:
        for line in old_file:
            if search_term in line:
                new_file.write(f'{search_term}{site_index}\n')
            else:
                new_file.write(line)
    copymode(old_file_path, new_file_path)
    move(new_file_path, old_file_path)
