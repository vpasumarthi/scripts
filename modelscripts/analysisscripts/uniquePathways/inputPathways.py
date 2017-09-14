#!/usr/bin/env python

from generateUniquePathways import generateUniquePathways
import os

# Frequently modified input parameters:
cutoffDistKey = 'Zr:Zr'
cutoff = 6
base = 0.05
prec = 2
inputCoordinateFileName = 'CONTCAR'

cwd = os.path.dirname(os.path.realpath(__file__))
inputFileLocation = os.path.join(cwd, inputCoordinateFileName)

generateUniquePathways(inputFileLocation, cutoffDistKey, cutoff, base, prec, cwd)