#!/usr/bin/env python

from generateUniquePathways import generateUniquePathways
import numpy as np
import os

# Frequently modified input parameters:
inputCoordinateFileName = 'POSCAR'
cutoffDistKey = 'O:O'
avoidElementType = 'O' # 'O' for cutoffDistKey = 'O:O'; '' for V
neighborCutoff = 4
bridgeCutoff = 2.51 # 2.62 for Zr:Zr; 3.50 for S:S; 2.51 for O:O
roundLatticeParameters = {} #{'base': 0.05, 'prec': 2}
printPathwayList = 1
printEquivalency = 0
equivalencyPrec = 4
pathwayPrec = 5

# unitCellCenterSiteClassList = unitCellNeighborSiteClassList = np.array([1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1]) # 'S:S'
unitCellCenterSiteClassList = unitCellNeighborSiteClassList = np.array([1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1]) # 'O:O'
classList = [unitCellCenterSiteClassList, unitCellNeighborSiteClassList]
# classList = []

desiredCoordinateParameters = {} # {'desiredSystemSize': np.array([2, 2, 2]), 'distList': [2.96781, 2.99576], 'prec': 5}

cwd = os.path.dirname(os.path.realpath(__file__))
inputFileLocation = os.path.join(cwd, inputCoordinateFileName)

generateUniquePathways(inputFileLocation, cutoffDistKey, neighborCutoff, bridgeCutoff, cwd, pathwayPrec, equivalencyPrec, classList, 
                       avoidElementType, roundLatticeParameters, printPathwayList, printEquivalency, desiredCoordinateParameters)