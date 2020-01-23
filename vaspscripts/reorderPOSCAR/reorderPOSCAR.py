# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as np


def readPOSCAR(srcFilePath, fileFormatIndex):
    # fileFormatIndex: 0=VASP; 1=VESTA
    latticeMatrix = np.zeros((3, 3))
    latticeParameterIndex = 0
    latticeParametersLineRange = range(3, 6)
    elementTypesLineNumber = 6 * fileFormatIndex
    numElementsLineNumber = 6 + fileFormatIndex
    coordStartLineNumber = 8 + fileFormatIndex
    inputFile = open(srcFilePath, 'r')
    for lineIndex, line in enumerate(inputFile):
        lineNumber = lineIndex + 1
        if lineNumber == 1 and not fileFormatIndex:
            elementTypes = line[:-1].split()
        elif lineNumber in latticeParametersLineRange:
            latticeMatrix[latticeParameterIndex, :] = np.fromstring(line,
                                                                    sep=' ')
            latticeParameterIndex += 1
        elif (lineNumber == elementTypesLineNumber
              and 'elementTypes' not in locals()):
            elementTypes = line.split()
        elif lineNumber == numElementsLineNumber:
            nElements = np.fromstring(line, dtype=int, sep=' ')
            totalElements = nElements.sum()
            fractionalCoords = np.zeros((totalElements, 3))
            elementIndex = 0
        elif ((lineNumber >= coordStartLineNumber)
              and (elementIndex < totalElements)):
            fractionalCoords[elementIndex, :] = np.fromstring(line, sep=' ')
            elementIndex += 1
    inputFile.close()
    POSCAR_INFO = np.array(
                    [latticeMatrix, elementTypes, nElements, fractionalCoords],
                    dtype=object)
    return POSCAR_INFO


def sort_by_layer(srcFilePath, fileFormatIndex, sortElementType, coordPrec=10,
                  diffPrec=8):
    coordStartLineNumber = 8 + fileFormatIndex
    [_, elementTypes, nElements, fractionalCoords] = readPOSCAR(
                                                srcFilePath, fileFormatIndex)
    elementTypes_consolidated = []
    uniqueElementTypes = set(elementTypes)
    numUniqueElementTypes = len(uniqueElementTypes)
    nElements_consolidated = np.zeros(numUniqueElementTypes, int)
    uniqueElementTypeIndex = -1
    for elementTypeIndex, elementType in enumerate(elementTypes):
        if elementType not in elementTypes_consolidated:
            uniqueElementTypeIndex += 1
            elementTypes_consolidated.append(elementType)
        nElements_consolidated[uniqueElementTypeIndex] += nElements[
                                                            elementTypeIndex]
    sortElementTypeIndex = elementTypes.index(sortElementType)
    numElements_sortElementType = nElements_consolidated[sortElementTypeIndex]
    startIndex = sum(nElements_consolidated[:sortElementTypeIndex])
    endIndex = startIndex + nElements[sortElementTypeIndex]
    sortElementTypeFractCoords = fractionalCoords[startIndex:endIndex]
    z_wise_reorderIndices = sortElementTypeFractCoords[:, 2].argsort()
    sortElementTypeFractCoords = sortElementTypeFractCoords[
                                                        z_wise_reorderIndices]
    zCoordinates = sortElementTypeFractCoords[:, 2]
    uniqueZCoordinates = np.unique(np.round(zCoordinates, coordPrec))
    diffZCoordinates = np.diff(uniqueZCoordinates)
    uniqueDiff = np.unique(np.round(diffZCoordinates, diffPrec))
    if len(uniqueDiff) == 2:
        layerCutoff = np.mean(uniqueDiff)
    layerIndices = np.zeros(numElements_sortElementType)
    for elementIndex in range(numElements_sortElementType):
        if elementIndex == 0:
            layerIndices[elementIndex] = 1
        else:
            zDiff = zCoordinates[elementIndex] - zCoordinates[elementIndex - 1]
            if zDiff < layerCutoff:
                layerIndices[elementIndex] = layerIndices[elementIndex - 1]
            else:
                layerIndices[elementIndex] = -layerIndices[elementIndex - 1]
    layer_wise_reorderIndices = layerIndices.argsort()
    sortElementTypeFractCoords = sortElementTypeFractCoords[
                                                    layer_wise_reorderIndices]
    print(np.unique(layerIndices, return_counts=True))
    headStart = (coordStartLineNumber - 1
                 + nElements_consolidated[:sortElementTypeIndex].sum())
    lineIndices = list(
        headStart + np.arange(numElements_sortElementType))
    newCoordinateList = np.copy(sortElementTypeFractCoords)
    newCoordinateList = np.asarray(newCoordinateList)
    writePOSCAR(srcFilePath, fileFormatIndex, lineIndices, newCoordinateList)
    return None


def writePOSCAR(srcFilePath, fileFormatIndex, lineIndices, newCoordinateList):
    dstFilePath = srcFilePath + '.out'
    srcFile = open(srcFilePath, 'r')
    open(dstFilePath, 'w').close()
    dstFile = open(dstFilePath, 'a')
    # fileFormatIndex: 0=VASP; 1=VESTA
    for lineIndex, line in enumerate(srcFile):
        if lineIndex in lineIndices:
            neighborIndex = lineIndices.index(lineIndex)
            if fileFormatIndex == 0:
                line = (
                    ''.join([
                        ' ' * 2,
                        '%18.16f' % newCoordinateList[neighborIndex][0],
                        ' ' * 2,
                        '%18.16f' % newCoordinateList[neighborIndex][1],
                        ' ' * 2,
                        '%18.16f' % newCoordinateList[neighborIndex][2]])
                    + '\n')
            elif fileFormatIndex == 1:
                line = (
                    ''.join([
                        ' ' * 5,
                        '%11.9f' % newCoordinateList[neighborIndex][0],
                        ' ' * 9,
                        '%11.9f' % newCoordinateList[neighborIndex][1],
                        ' ' * 9,
                        '%11.9f' % newCoordinateList[neighborIndex][2]])
                    + '\n')
        dstFile.write(line)
    srcFile.close()
    dstFile.close()
