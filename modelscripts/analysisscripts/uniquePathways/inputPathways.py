#!/usr/bin/env python

from generateUniquePathways import generateUniquePathways
import numpy as np
import os

# Frequently modified input parameters:
cutoffDistKey = 'O:O'
cutoff = 4
base = 0.05
prec = 2
nDim = 3
cwd = os.path.dirname(os.path.realpath(__file__))

inputCoordinateFileName = 'POSCAR'
inputFileLocation = os.path.join(cwd, inputCoordinateFileName)

generateUniquePathways(inputFileLocation, cutoffDistKey, cutoff, base, prec, cwd)