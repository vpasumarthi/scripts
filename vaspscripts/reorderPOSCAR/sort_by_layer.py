# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from reorderPOSCAR import sort_by_layer

import os

# Frequently modified input parameters:
inputCoordinateFileName = 'POSCAR'
# fileFormatIndex: 0=VASP; 1=VESTA
fileFormatIndex = 0
sortElementType = 'Fe'

cwd = os.path.dirname(os.path.realpath(__file__))
inputFileLocation = os.path.join(cwd, inputCoordinateFileName)

sort_by_layer(inputFileLocation, fileFormatIndex, sortElementType)
