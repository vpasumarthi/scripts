#!/usr/bin/env python

from pathlib import Path

from PyCT.material_preprod import material_preprod

cwd = Path.cwd()
material_preprod(cwd)
