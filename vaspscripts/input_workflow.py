# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# File to generate input files for charge localization jobs

from pathlib import Path
from shutil import copy, copymode, move
from os import chdir

import numpy as np


site_indices = np.loadtxt('site_indices.dat', dtype=int)
cwd = Path.cwd()
partial_run_index = 1
copy_file_names = ["POSCAR", "POTCAR", "INCAR", "KPOINTS", "slurmscript", "generate_bond_distortion.py"]
site_index_search_term = 'localized_site_number = '
exp_num_diff_lines = 4
ref_site_index = 20
ref_partial_run_index = 1
localized_site_element = 'V'
system_search_term = "SYSTEM          = "
magmom_search_term = "MAGMOM          = "
slurm_search_term = "#SBATCH --job-name="

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
            if site_index_search_term in line:
                new_file.write(f'{site_index_search_term}{site_index}\n')
            else:
                new_file.write(line)
    copymode(old_file_path, new_file_path)
    move(new_file_path, old_file_path)

    # Identify element distribution
    coord_file_path = work_dir_path / "POSCAR"
    with open(coord_file_path) as coord_file:
        for line_index, line in enumerate(coord_file):
            if line_index == 5:
                element_list = line.split()
            elif line_index == 6:
                num_elements_list = [int(num_elements) for num_elements in line.split()]

    # Generate bond distortion
    chdir(work_dir_path)
    exec(open(old_file_path).read())
    num_diff_lines = 0
    with open("POSCAR") as old_file, open("POSCAR.out") as new_file:
        for line1, line2 in zip(old_file, new_file):
            if line1 != line2 and '.' in line1 and '.' in line2:
                num_diff_lines += 1
    if num_diff_lines == exp_num_diff_lines:
        move("POSCAR.out", "POSCAR")
        ref_site_dir_path = cwd / f'W36_ox6_V{ref_site_index:02}_ox4'
        ref_work_dir_path = ref_site_dir_path / f'partial_run{ref_partial_run_index:02}'
        ref_coord_file_path = ref_work_dir_path / "POSCAR"
        old_coord_file_path = work_dir_path / "POSCAR"
        new_coord_file_path = work_dir_path / "POSCAR.new"
        with open(ref_coord_file_path) as ref_coords, open(old_coord_file_path) as old_coords, open(new_coord_file_path, 'w') as new_coords:
            for line_index, (line1, line2) in enumerate(zip(ref_coords, old_coords)):
                if line_index == 5 or line_index == 6:
                    new_coords.write(line1)
                else:
                    new_coords.write(line2)
        move(new_coord_file_path, old_coord_file_path)
    else:
        print(f'Unexpected bond distortion output in {work_dir_path}')
    chdir(cwd)

    # Modify INCAR file
    old_incar_file_path = work_dir_path / "INCAR"
    new_incar_file_path = work_dir_path / "INCAR.new"
    with open(old_incar_file_path) as old_incar, open(new_incar_file_path, 'w') as new_incar:
        for line in old_incar:
            if system_search_term in line:
                new_incar.write(f'{system_search_term}BVO-e-W36-ox6-V{site_index:02}-ox4  # System name\n')
            elif magmom_search_term in line:
                magmom_line = magmom_search_term
                for element_index, element_type in enumerate(element_list):
                    if element_type == localized_site_element:
                        magmom_line = f'{magmom_line} {site_index - 1}*0.0 1*1.0 {num_elements_list[element_index] - site_index}*0.0'
                    else:
                        magmom_line = f'{magmom_line} {num_elements_list[element_index]}*0.0'
                new_incar.write(magmom_line)
            else:
                new_incar.write(line)
    move(new_incar_file_path, old_incar_file_path)

    # Modify slurmscript
    old_slurm_file_path = work_dir_path / "slurmscript"
    new_slurm_file_path = work_dir_path / "slurmscript.new"
    with open(old_slurm_file_path) as old_slurm, open(new_slurm_file_path, 'w') as new_slurm:
        for line in old_slurm:
            if slurm_search_term in line:
                slurm_line = f'#SBATCH --job-name=BVO552-W36-ox6-V{site_index:02}-ox4-K111\n'
                new_slurm.write(slurm_line)
            else:
                new_slurm.write(line)
    move(new_slurm_file_path, old_slurm_file_path)

