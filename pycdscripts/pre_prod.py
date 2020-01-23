# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from pathlib import Path

from PyCT.material_preprod import material_preprod
from pyctscripts.parallel_input_files import generate_parallel_input_files

cwd = Path.cwd()
material_preprod(cwd)
generate_parallel_input_files(cwd)
