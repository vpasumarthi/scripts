#!/usr/bin/env python

from generateUniquePathways import generateUniquePathways
import numpy as np
import os

# Frequently modified input parameters:
cutoffDistKey = 'S:S'
cutoff = 6
base = 0.05
prec = 2
inputCoordinateFileName = 'CONTCAR'
# class list definitions within 
unitCellCenterSiteClassList = unitCellNeighborSiteClassList = np.array([1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1])
# unitCellCenterSiteClassList = unitCellNeighborSiteClassList = np.array([1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1]) # 'O:O'
classList = [unitCellCenterSiteClassList, unitCellNeighborSiteClassList]

cwd = os.path.dirname(os.path.realpath(__file__))
inputFileLocation = os.path.join(cwd, inputCoordinateFileName)

generateUniquePathways(inputFileLocation, cutoffDistKey, cutoff, base, prec, cwd, classList)