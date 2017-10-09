#!/usr/bin/env python

from chargeDistortion import chargeDistortion
import os

# Frequently modified input parameters:
inputCoordinateFileName = 'POSCAR'
# fileFormatIndex: 0=VASP; 1=VESTA
fileFormatIndex = 0
localizedElementType = 'O'
localizedSiteNumber = 86
neighborElementTypeList = ['V', 'Bi']
neighborCutoffList = [1.80, 2.50]  # angstrom; V-O: 1.80; O-V: 1.80; O-Bi: 2.50
stretchPercentList = [8.00, 8.00]  # electron: 0.15; hole: 0.20

cwd = os.path.dirname(os.path.realpath(__file__))
inputFileLocation = os.path.join(cwd, inputCoordinateFileName)

chargeDistortion(inputFileLocation, fileFormatIndex, localizedElementType,
                 localizedSiteNumber, neighborElementTypeList,
                 neighborCutoffList, stretchPercentList)
