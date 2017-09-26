#/usr/bin/env python

from chargeDistortion import chargeDistortion
import numpy as np
from shutil import copymode
import os

# Frequently modified input parameters:
inputCoordinateFileName = 'POSCAR'
localizedElementType = 'O'
localizedSiteNumber = 124
neighborElementTypeList = ['V', 'Bi']
neighborCutoffList = [1.80, 2.50] # angstrom; V-O: 1.80; O-V: 1.80; O-Bi: 2.50 
stretchLengthList = [0.15, 0.20] # electron: 0.15; hole: 0.20

cwd = os.path.dirname(os.path.realpath(__file__))
inputFileLocation = os.path.join(cwd, inputCoordinateFileName)

chargeDistortion(inputFileLocation, localizedElementType, localizedSiteNumber, neighborElementTypeList, 
                 neighborCutoffList, stretchLengthList)