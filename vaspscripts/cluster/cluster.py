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
                    [latticeMatrix, elementTypes, nElements, fractionalCoords,
                     fileFormat], dtype=object)
    return POSCAR_INFO


def cluster(srcFilePath, siteElementTypeList,
            siteNumberList, bridgeSearchDepth, terminatingElementType,
            terminatingBondDistance):
    [latticeMatrix, elementTypes, nElements, fractionalCoords, fileFormat] = (
                                                    readPOSCAR(srcFilePath))
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
    import pdb; pdb.set_trace()
    localizedElementTypeIndex = elementTypes_consolidated.index(
                                                        localizedElementType)
    localizedSiteCoords = (
        fractionalCoords[nElements_consolidated[
                                            :localizedElementTypeIndex].sum()
                         + localizedSiteNumber - 1])

    # generate array of unit cell translational coordinates
    pbc = np.ones(3, int)
    numCells = 3**sum(pbc)
    xRange = range(-1, 2) if pbc[0] == 1 else [0]
    yRange = range(-1, 2) if pbc[1] == 1 else [0]
    zRange = range(-1, 2) if pbc[2] == 1 else [0]
    cellTranslationalCoords = np.zeros((numCells, 3))  # Initialization
    index = 0
    for xOffset in xRange:
        for yOffset in yRange:
            for zOffset in zRange:
                cellTranslationalCoords[index] = np.array([xOffset,
                                                           yOffset,
                                                           zOffset])
                index += 1
    localizedSiteCoords_imageconsolidated = (localizedSiteCoords
                                             + cellTranslationalCoords)
    lineIndices = []
    newCoordinateList = []
    for distortElementTypeIndex, distortElementType in enumerate(
                                                    neighborElementTypeList):
        neighborElementTypeIndex = elementTypes_consolidated.index(
                                                            distortElementType)
        neighborSiteCoords = fractionalCoords[
                    nElements_consolidated[:neighborElementTypeIndex].sum()
                    + range(nElements_consolidated[neighborElementTypeIndex])]
        neighborCutoffDistLimits = [
                                0,
                                neighborCutoffList[distortElementTypeIndex]]

        # generate neighbor list
        neighborList = []
        centerSiteCoordList = []
        for neighborSiteIndex, neighborSiteCoord in enumerate(
                                                        neighborSiteCoords):
            latticeDirections = (localizedSiteCoords_imageconsolidated
                                 - neighborSiteCoord)
            minDisp = np.linalg.norm(np.sum(latticeMatrix, axis=0))
            for iCell in range(numCells):
                displacement = np.linalg.norm(np.dot(latticeDirections[iCell],
                                                     latticeMatrix))
                if displacement < minDisp:
                    minDisp = displacement
                    centerSiteCoords = localizedSiteCoords_imageconsolidated[
                                                                        iCell]
            if (neighborCutoffDistLimits[0] < minDisp
                    <= neighborCutoffDistLimits[1]):
                neighborList.append(neighborSiteIndex)
                centerSiteCoordList.append(centerSiteCoords)

        # generate distortion
        numNeighbors = len(neighborList)
        headStart = (coordStartLineNumber - 1
                     + nElements_consolidated[:neighborElementTypeIndex].sum())
        for iNeighbor in range(numNeighbors):
            latticeDirection = (neighborSiteCoords[neighborList[iNeighbor]]
                                - centerSiteCoordList[iNeighbor])
            displacement = np.linalg.norm(np.dot(latticeDirection,
                                                 latticeMatrix))
            unitVector = latticeDirection / displacement
            lineIndices.append(headStart + neighborList[iNeighbor])
            newCoordinateList.append(
                centerSiteCoordList[iNeighbor] + unitVector
                * (displacement
                   * (1 + stretchPercentList[distortElementTypeIndex] / 100)))
    newCoordinateList = np.asarray(newCoordinateList)
    writePOSCAR(srcFilePath, fileFormatIndex, elementTypes_cluster,
                nElements_cluster, coordinates_cluster)
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
