#!/usr/bin/env python

from generateUniquePathways import generateUniquePathways
import numpy as np
import os

# Frequently modified input parameters:
inputCoordinateFileName = 'CONTCAR'
cutoffDistKey = 'S:S'
avoidElementType = 'S' # 'O' for cutoffDistKey = 'O:O'; '' for V
neighborCutoff = 6
bridgeCutoff = 2.62 # 2.62 for Zr:Zr; 3.50 for S:S
roundLattice = 0
printStack = 0
printEquivalency = 1

# Input Parameters:
base = 0.05
prec = 2

unitCellCenterSiteClassList = unitCellNeighborSiteClassList = np.array([1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1])
# unitCellCenterSiteClassList = unitCellNeighborSiteClassList = np.array([1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1]) # 'O:O'
classList = [unitCellCenterSiteClassList, unitCellNeighborSiteClassList]

cwd = os.path.dirname(os.path.realpath(__file__))
inputFileLocation = os.path.join(cwd, inputCoordinateFileName)

generateUniquePathways(inputFileLocation, cutoffDistKey, neighborCutoff, bridgeCutoff, base, prec, cwd, 
                       classList, avoidElementType, roundLattice, printStack, printEquivalency)