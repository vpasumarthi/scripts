#!/usr/bin/env python

from cluster import cluster
import os

# Frequently modified input parameters:
inputCoordinateFileName = 'CONTCAR'
siteElementTypeList = ['V', 'V']
siteNumberList = [14, 32]
bridgeSearchDepth = 2
terminatingElementType = 'H'
bondLimits = {'Bi:O': 2.6608, 'V:O': 2.84939}
terminatingBondDistance = 0.96  # angstrom

cwd = os.path.dirname(os.path.realpath(__file__))
inputFileLocation = os.path.join(cwd, inputCoordinateFileName)

cluster(inputFileLocation, siteElementTypeList,
        siteNumberList, bondLimits, bridgeSearchDepth, terminatingElementType,
        terminatingBondDistance)
