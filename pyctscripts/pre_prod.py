#!/usr/bin/env python

from pathlib import Path

from PyCT.material_preprod import material_preprod
from pyctscripts.generate_parallel_input_files import parallel_input_files

cwd = Path.cwd()
material_preprod(cwd)
parallel_input_files(cwd)
