import numpy as np
from shutil import copymode
import os

def readPOSCAR(inputFilePath):
    latticeMatrix = np.zeros((3, 3))
    latticeParameterIndex = 0
    latticeParametersLineRange = range(3, 6)
    inputFile = open(inputFilePath, 'r')
    for lineIndex, line in enumerate(inputFile):
        lineNumber = lineIndex + 1
        if lineNumber in latticeParametersLineRange:
            latticeMatrix[latticeParameterIndex, :] = np.fromstring(line, sep=' ')
            latticeParameterIndex += 1
        elif lineNumber == 6:
            elementTypes = line.split()
        elif lineNumber == 7:
            nElements = np.fromstring(line, dtype=int, sep=' ')
            totalElements = nElements.sum()
            fractionalCoords = np.zeros((totalElements, 3))
            elementIndex = 0
        elif lineNumber > 8 and elementIndex < totalElements:
            fractionalCoords[elementIndex, :] = np.fromstring(line, sep=' ')
            elementIndex += 1
    inputFile.close()
    POSCAR_INFO = np.array([latticeMatrix, elementTypes, nElements, fractionalCoords], dtype=object)
    return POSCAR_INFO

def writePOSCAR(outputFilePath):
    pass

def chargeDistortion(inputFileLocation, localizedElementType, localizedSiteNumber, neighborElementType, 
                     neighborCutoff, stretchLength):
    [latticeMatrix, elementTypes, nElements, fractionalCoords] = readPOSCAR(inputFileLocation)
    elementTypes_consolidated = []
    uniqueElementTypes = set(elementTypes)
    numUniqueElementTypes = len(uniqueElementTypes)
    nElements_consolidated = np.zeros(numUniqueElementTypes, int)
    uniqueElementTypeIndex = -1
    for elementTypeIndex, elementType in enumerate(elementTypes):
        if elementType not in elementTypes_consolidated:
            uniqueElementTypeIndex += 1
            elementTypes_consolidated.append(elementType)
        nElements_consolidated[uniqueElementTypeIndex] += nElements[elementTypeIndex]
    localizedElementTypeIndex = elementTypes_consolidated.index(localizedElementType)
    localizedSiteCoords = fractionalCoords[nElements_consolidated[:localizedElementTypeIndex].sum() + localizedSiteNumber - 1]
    neighborElementTypeIndex = elementTypes_consolidated.index(neighborElementType)
    neighborSiteCoords = fractionalCoords[nElements_consolidated[:neighborElementTypeIndex].sum() + range(nElements_consolidated[neighborElementTypeIndex])]
    neighborCutoffDistLimits = [0, neighborCutoff]
    
    # generate array of unit cell translational coordinates
    pbc = np.ones(3, int)
    numCells = 3**sum(pbc)
    xRange = range(-1, 2) if pbc[0] == 1 else [0]
    yRange = range(-1, 2) if pbc[1] == 1 else [0]
    zRange = range(-1, 2) if pbc[2] == 1 else [0]
    cellTranslationalCoords = np.zeros((numCells, 3)) # Initialization
    index = 0
    for xOffset in xRange:
        for yOffset in yRange:
            for zOffset in zRange:
                cellTranslationalCoords[index] = np.array([xOffset, yOffset, zOffset])
                index += 1
    localizedSiteCoords_imageconsolidated = localizedSiteCoords + cellTranslationalCoords
    
    # generate neighbor list
    neighborList = []
    centerSiteCoordList = []
    for neighborSiteIndex, neighborSiteCoord in enumerate(neighborSiteCoords):
        latticeDirections = localizedSiteCoords_imageconsolidated - neighborSiteCoord
        minDisp = np.linalg.norm(np.sum(latticeMatrix, axis=0))
        for iCell in range(numCells):
            displacement = np.linalg.norm(np.dot(latticeDirections[iCell], latticeMatrix))
            if displacement < minDisp:
                minDisp = displacement
                centerSiteCoords = localizedSiteCoords_imageconsolidated[iCell]
        if neighborCutoffDistLimits[0] < minDisp <= neighborCutoffDistLimits[1]:
            neighborList.append(neighborSiteIndex)
            centerSiteCoordList.append(centerSiteCoords)

    # generate distortion
    newNeighborSiteCoords = np.copy(neighborSiteCoords)
    numNeighbors = len(neighborList)
    for iNeighbor in range(numNeighbors):
        latticeDirection = neighborSiteCoords[neighborList[iNeighbor]] - centerSiteCoordList[iNeighbor]
        displacement = np.linalg.norm(np.dot(latticeDirection, latticeMatrix))
        unitVector = latticeDirection / displacement
        newNeighborSiteCoords[neighborList[iNeighbor]] = centerSiteCoordList[iNeighbor] + unitVector * (displacement + stretchLength)
    
    import pdb; pdb.set_trace()
    writePOSCAR(outputFilePath)
    return