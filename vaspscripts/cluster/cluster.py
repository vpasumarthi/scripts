#!/usr/bin/env python

import numpy as np


def readPOSCAR(srcFilePath):
    inputFile = open(srcFilePath, 'r')
    elementTypesLineNumber = 6
    for lineIndex, line in enumerate(inputFile):
        lineNumber = lineIndex + 1
        if lineNumber == elementTypesLineNumber:
            prospectiveElementTypes = line[:-1].split()
            elementTypesExist = prospectiveElementTypes[0].isalpha()
            break
    inputFile.close()

    inputFile = open(srcFilePath, 'r')
    latticeMatrix = np.zeros((3, 3))
    latticeParameterIndex = 0
    latticeParametersLineRange = range(3, 6)
    elementTypesLineNumber = 6 * elementTypesExist
    numElementsLineNumber = 6 + elementTypesExist
    coordinateTypeNumber = 7 + elementTypesExist
    coordStartLineNumber = 8 + elementTypesExist
    for lineIndex, line in enumerate(inputFile):
        lineNumber = lineIndex + 1
        if lineNumber == 1 and not elementTypesExist:
            elementTypes = line.split()
            print(elementTypes)
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
        elif lineNumber == coordinateTypeNumber:
            coordinateType = line.split()[0]
        elif lineNumber == coordStartLineNumber:
            elementIndex = 0
            fractionalCoords[elementIndex, :] = np.fromstring(line, sep=' ')
            coordinateStringLength = len(line.split()[0])
            if coordinateStringLength == 18:
                fileFormat = 'VASP'
            elif coordinateStringLength == 11:
                fileFormat = 'VESTA'
            else:
                fileFormat = 'unknown'
        elif ((lineNumber > coordStartLineNumber)
              and (elementIndex < totalElements - 1)):
            elementIndex += 1
            fractionalCoords[elementIndex, :] = np.fromstring(line, sep=' ')
    inputFile.close()
    POSCAR_INFO = np.array(
                    [latticeMatrix, elementTypes, nElements, coordinateType,
                     fractionalCoords, fileFormat], dtype=object)
    return POSCAR_INFO


def cluster(srcFilePath, siteElementTypeList, siteNumberList, bondLimits,
            bridgeSearchDepth, terminatingElementType,
            terminatingBondDistance):
    [latticeMatrix, elementTypes, nElements, coordinateType, fractionalCoords,
     fileFormat] = readPOSCAR(srcFilePath)
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
    elementWiseCoordinateList = [[] for _ in range(numUniqueElementTypes)]
    for siteIndex, siteElementType in enumerate(siteElementTypeList):
        siteElementTypeIndex = elementTypes_consolidated.index(siteElementType)
        siteCoordinates = fractionalCoords[
                            sum(nElements_consolidated[:siteElementTypeIndex])
                            + siteNumberList[siteIndex] - 1]
        elementWiseCoordinateList[siteElementTypeIndex].append(siteCoordinates)

    elementTypes_cluster = []
    nElements_cluster = []
    numCoordinates = 0
    for elementIndex in range(numUniqueElementTypes):
        elementIndexNumElements = len(elementWiseCoordinateList[elementIndex])
        if elementIndexNumElements != 0:
            elementTypes_cluster.append(
                                    elementTypes_consolidated[elementIndex])
            nElements_cluster.append(elementIndexNumElements)
            numCoordinates += elementIndexNumElements
    coordinates_cluster = [
            coordinates for elementWiseCoordinates in elementWiseCoordinateList
            for coordinates in elementWiseCoordinates]
    writePOSCAR(srcFilePath, fileFormat, elementTypes_cluster,
                nElements_cluster, coordinateType, coordinates_cluster)
    return None


def writePOSCAR(srcFilePath, fileFormat, elementTypes_cluster,
                nElements_cluster, coordinateType, coordinates_cluster):
    unmodifiedLineNumberLimit = 5
    dstFilePath = srcFilePath + '.out'
    srcFile = open(srcFilePath, 'r')
    open(dstFilePath, 'w').close()
    dstFile = open(dstFilePath, 'a')
    for lineIndex, line in enumerate(srcFile):
        lineNumber = lineIndex + 1
        if lineNumber <= unmodifiedLineNumberLimit:
            dstFile.write(line)
        else:
            break
    srcFile.close()
    elementTypesLine = (' ' * 3 + (' ' * 4).join(elementTypes_cluster) + '\n')
    dstFile.write(elementTypesLine)
    # import pdb; pdb.set_trace()
    nElementsLine = (' ' * 3 + (' ' * 4).join(map(str, nElements_cluster))
                     + '\n')
    dstFile.write(nElementsLine)
    dstFile.write(coordinateType + '\n')
    for elementCoordinates in coordinates_cluster:
        if fileFormat == 'VASP' or fileFormat == 'unknown':
            line = (
                ''.join([
                    ' ' * 2,
                    '%18.16f' % elementCoordinates[0],
                    ' ' * 2,
                    '%18.16f' % elementCoordinates[1],
                    ' ' * 2,
                    '%18.16f' % elementCoordinates[2]])
                + '\n')
        elif fileFormat == 'VESTA':
            line = (
                ''.join([
                    ' ' * 5,
                    '%11.9f' % elementCoordinates[0],
                    ' ' * 9,
                    '%11.9f' % elementCoordinates[1],
                    ' ' * 9,
                    '%11.9f' % elementCoordinates[2]])
                + '\n')
        dstFile.write(line)
    dstFile.close()
