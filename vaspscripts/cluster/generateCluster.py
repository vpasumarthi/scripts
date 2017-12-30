#!/usr/bin/env python

from cluster import cluster
import os

# Frequently modified input parameters:
inputCoordinateFileName = 'CONTCAR'
siteElementTypeList = ['V', 'V']
siteNumberList = [38, 89]
bondLimits = {'Bi:O': 2.60, 'V:O': 1.80}
terminatingElementType = 'H'
terminatingBondDistance = 0.96  # angstrom

cwd = os.path.dirname(os.path.realpath(__file__))
inputFileLocation = os.path.join(cwd, inputCoordinateFileName)

cluster(inputFileLocation, siteElementTypeList,
        siteNumberList, bondLimits, terminatingElementType,
        terminatingBondDistance)
