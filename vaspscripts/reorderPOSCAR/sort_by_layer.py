#!/usr/bin/env python

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
