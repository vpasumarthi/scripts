#!/usr/bin/env python

from vaspscripts.cluster import cluster
import os

# Frequently modified input parameters:
inputCoordinateFileName = 'CONTCAR'
outputCoordinateFileName = 'POSCAR'
siteIndexList = [1134, 1037]  # site number - 1
bondLimits = {'Bi:O': 2.60, 'V:O': 1.80}
terminatingElementType = 'H'
terminatingBondDistance = 0.96  # angstrom
oxidationList = {'Bi': +3, 'V': +5, 'O': -2, 'H': +1}
bridgeSearchDepth = 1
chargeNeutral = 1
prec = 3

cwd = os.path.dirname(os.path.realpath(__file__))
inputFileLocation = os.path.join(cwd, inputCoordinateFileName)

cluster(inputFileLocation, outputCoordinateFileName, siteIndexList, bondLimits,
        terminatingElementType, terminatingBondDistance, oxidationList,
        bridgeSearchDepth, chargeNeutral, prec)
