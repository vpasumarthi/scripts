#!/usr/bin/env python

from pathlib import Path

from PyCT.material_preprod import material_preprod
from pyctscripts.parallel_input_files import generate_parallel_input_files

cwd = Path.cwd()
material_preprod(cwd)
parallel_input_files(cwd)
