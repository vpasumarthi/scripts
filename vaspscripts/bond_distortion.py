#!/usr/bin/env python

import numpy as np

from PyCT.io import read_poscar, write_poscar


def bond_distortion(srcFilePath, localizedElementType, localizedSiteNumber,
                    neighborElementTypeList, neighborCutoffList,
                    stretchPercentList):
    poscar_info = read_poscar(srcFilePath)
    latticeMatrix = poscar_info['lattice_matrix']
    elementTypes = poscar_info['element_types']
    nElements = poscar_info['num_elements']
    coordinate_type = poscar_info['coordinate_type']
    unit_cell_coords = poscar_info['coordinates']
    if coordinate_type == 'Direct':
        fractionalCoords = unit_cell_coords
    elif coordinate_type == 'Cartesian':
        fractionalCoords = np.dot(unit_cell_coords,
                                  np.linalg.inv(latticeMatrix))
    file_format = poscar_info['file_format']

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
        headStart = nElements_consolidated[:neighborElementTypeIndex].sum()
        for iNeighbor in range(numNeighbors):
            latticeDirection = (neighborSiteCoords[neighborList[iNeighbor]]
                                - centerSiteCoordList[iNeighbor])
            displacement = np.linalg.norm(np.dot(latticeDirection,
                                                 latticeMatrix))
            unitVector = latticeDirection / displacement
            index = headStart + neighborList[iNeighbor]
            newCoordinate = (
                centerSiteCoordList[iNeighbor] + unitVector
                * (displacement
                   * (1 + stretchPercentList[distortElementTypeIndex] / 100)))
            fractionalCoords[index] = newCoordinate
    dstFilePath = srcFilePath + '.out'
    write_poscar(srcFilePath, dstFilePath, file_format,
                 elementTypes, nElements, coordinate_type,
                 fractionalCoords)
    return None
