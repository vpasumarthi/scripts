#!/usr/bin/env python

from pathlib import Path

from PyCT.material_preprod import material_preprod
from pyctscripts.generate_symlink import generate_symlink

cwd = Path.cwd()
material_preprod(cwd)
generate_symlink(cwd)
