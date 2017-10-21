#!/usr/bin/env python

from uniquePathways import uniquePathways
import numpy as np
import os

# Frequently modified input parameters:
inputCoordinateFileName = 'POSCAR'
cutoffDistKey = 'O:O'
avoidElementType = 'O'  # 'O' for cutoffDistKey = 'O:O'; '' for V
neighborCutoff = 4
bridgeCutoff = 2.51  # 2.62 for Zr:Zr; 3.50 for S:S; 2.51 for O:O
# {'base': 0.05, 'prec': 2} or {}
roundLatticeParameters = {'base': 0.05, 'prec': 2}
printPathwayList = 1
printEquivalency = 0
equivalencyPrec = 4
pathwayPrec = 5

# 'S:S': np.array([1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1])
# 'O:O': np.array([1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1])
unitCellCenterSiteClassList = unitCellNeighborSiteClassList = \
                    np.array([1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1])
# 'V:V':
# classList = []
# 'O:O' or 'S:S':
# classList = [unitCellCenterSiteClassList, unitCellNeighborSiteClassList]
classList = [unitCellCenterSiteClassList, unitCellNeighborSiteClassList]

# desiredCoordinateParameters = {'desiredSystemSize': np.array([2, 2, 2]),
#                                'distList': [2.96781, 2.99576],
#                                'prec': 5}
# desiredCoordinateParameters = {}
desiredCoordinateParameters = {}

cwd = os.path.dirname(os.path.realpath(__file__))
inputFileLocation = os.path.join(cwd, inputCoordinateFileName)

uniquePathways(inputFileLocation, cutoffDistKey, neighborCutoff, bridgeCutoff,
               cwd, pathwayPrec, equivalencyPrec, classList, avoidElementType,
               roundLatticeParameters, printPathwayList, printEquivalency,
               desiredCoordinateParameters)
