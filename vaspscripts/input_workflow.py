#!/usr/bin/env python

# File to generate input files for charge localization jobs

from pathlib import Path
from shutil import copy

import numpy as np


site_indices = np.loadtxt('site_indices.dat', dtype=int)
cwd = Path.cwd()
partial_run_index = 1
copy_file_names = ["POSCAR", "POTCAR", "INCAR", "KPOINTS", "slurmscript", "generate_bond_distortion.py"]

for site_index in site_indices:
    # Generate input files
    site_dir_path = cwd / f'W36_ox6_V{site_index:02}_ox4'
    work_dir_path = site_dir_path / f'partial_run{partial_run_index:02}'
    Path.mkdir(work_dir_path, parents=True, exist_ok=True)

    for file_name in copy_file_names:
        copy(cwd / file_name, work_dir_path)
