#!/usr/bin/env python

# File to generate input files for charge localization jobs

from pathlib import Path

import numpy as np


site_indices = np.loadtxt('site_indices.dat', dtype=int)
cwd = Path.cwd()
partial_run_index = 1

for site_index in site_indices:
    site_dir_path = cwd / f'W36_ox6_V{site_index:02}_ox4'
    work_dir_path = site_dir_path / f'partial_run{partial_run_index:02}'
    Path.mkdir(work_dir_path, parents=True, exist_ok=True)
