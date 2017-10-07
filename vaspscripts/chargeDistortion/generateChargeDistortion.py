#!/usr/bin/env python

from chargeDistortion import chargeDistortion
import os

# Frequently modified input parameters:
inputCoordinateFileName = 'POSCAR'
localizedElementType = 'O'
localizedSiteNumber = 124
neighborElementTypeList = ['V', 'Bi']
neighborCutoffList = [2.00, 2.80]  # angstrom; V-O: 1.80; O-V: 1.80; O-Bi: 2.50
stretchPercentList = [8.00, 8.00]  # electron: 0.15; hole: 0.20

cwd = os.path.dirname(os.path.realpath(__file__))
inputFileLocation = os.path.join(cwd, inputCoordinateFileName)

chargeDistortion(inputFileLocation, localizedElementType, localizedSiteNumber,
                 neighborElementTypeList, neighborCutoffList,
                 stretchPercentList)
