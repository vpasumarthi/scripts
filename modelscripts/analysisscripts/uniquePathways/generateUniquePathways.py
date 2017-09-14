import numpy as np
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
	    	nElementsPerUnitCell = np.fromstring(line, dtype=int, sep=' ')
	    	totalElementsPerUnitCell = nElementsPerUnitCell.sum()
	    	fractionalUnitCellCoords = np.zeros((totalElementsPerUnitCell, 3))
	    	elementIndex = 0
	    elif lineNumber > 8 and elementIndex < totalElementsPerUnitCell:
	    	fractionalUnitCellCoords[elementIndex, :] = np.fromstring(line, sep=' ')
	    	elementIndex += 1
	inputFile.close()
	output = np.array([latticeMatrix, elementTypes, nElementsPerUnitCell, fractionalUnitCellCoords], dtype=object)
	return output
		
def generateQuantumIndices(systemSize, systemElementIndex, nElementsPerUnitCell):
	"""Returns the quantum indices of the element"""
	#assert systemElementIndex >= 0, 'System Element Index cannot be negative'
	#assert systemElementIndex < systemSize.prod() * self.material.totalElementsPerUnitCell, 'System Element Index out of range for the given system size'
	quantumIndices = np.zeros(5, dtype=int)#[0] * 5
	totalElementsPerUnitCell = nElementsPerUnitCell.sum()
	unitcellElementIndex = systemElementIndex % totalElementsPerUnitCell
	quantumIndices[3] = np.where(np.cumsum(nElementsPerUnitCell) >= (unitcellElementIndex + 1))[0][0]
	quantumIndices[4] = unitcellElementIndex - nElementsPerUnitCell[:quantumIndices[3]].sum()
	nFilledUnitCells = (systemElementIndex - unitcellElementIndex) / totalElementsPerUnitCell
	for index in range(3):
		quantumIndices[index] = nFilledUnitCells / systemSize[index+1:].prod()
		nFilledUnitCells -= quantumIndices[index] * systemSize[index+1:].prod()
	return quantumIndices

def generateUniquePathways(inputFileLocation, cutoffDistKey, neighborCutoff, bridgeCutoff, outdir, 
						classList=[], avoidElementType='', roundLatticeParameters={}, printPathwayList=0, 
						printEquivalency=0):
	""" generate unique pathways for the given set of element types"""
	neighborCutoffDistLimits = [0, neighborCutoff]
	bridgeCutoffDistLimits = [0, bridgeCutoff]
	[latticeMatrix, elementTypes, nElementsPerUnitCell, fractionalUnitCellCoords] = readPOSCAR(inputFileLocation)
	nElementTypes = len(elementTypes)
	totalElementsPerUnitCell = nElementsPerUnitCell.sum()
	elementTypeIndexList = np.repeat(np.arange(nElementTypes), nElementsPerUnitCell)
	if roundLatticeParameters:
		base = roundLatticeParameters['base']
		prec = roundLatticeParameters['prec']
	
	startIndex = 0
	for elementIndex in range(nElementTypes):
	    endIndex = startIndex + nElementsPerUnitCell[elementIndex] 
	    elementUnitCellCoords = fractionalUnitCellCoords[elementTypeIndexList==elementIndex]
	    fractionalUnitCellCoords[startIndex:endIndex] = elementUnitCellCoords[elementUnitCellCoords[:,2].argsort()]
	    startIndex = endIndex
	
	[centerElementType, neighborElementType] = cutoffDistKey.split(':')
	centerSiteElementTypeIndex = elementTypes.index(centerElementType) 
	neighborSiteElementTypeIndex = elementTypes.index(neighborElementType)
	
	pbc = np.ones(3, int)
	numCells = 3**sum(pbc)
	xRange = range(-1, 2) if pbc[0] == 1 else [0]
	yRange = range(-1, 2) if pbc[1] == 1 else [0]
	zRange = range(-1, 2) if pbc[2] == 1 else [0]
	localSystemSize = np.array([len(xRange), len(yRange), len(zRange)])
	unitcellTranslationalCoords = np.zeros((numCells, 3)) # Initialization
	index = 0
	for xOffset in xRange:
		for yOffset in yRange:
			for zOffset in zRange:
				unitcellTranslationalCoords[index] = np.array([xOffset, yOffset, zOffset])
				index += 1
	
	if avoidElementType:
		avoidElementTypeIndex = elementTypes.index(avoidElementType)
		systemElementIndexOffsetArray = (np.repeat(np.arange(0, totalElementsPerUnitCell * numCells, totalElementsPerUnitCell), 
												   nElementsPerUnitCell[centerSiteElementTypeIndex]))                
		avoidElementIndices = (np.tile(nElementsPerUnitCell[:avoidElementTypeIndex].sum() + 
									   np.arange(0, nElementsPerUnitCell[avoidElementTypeIndex]), numCells) + systemElementIndexOffsetArray)
	
	numCenterElements = nElementsPerUnitCell[centerSiteElementTypeIndex]
	pathwayList = np.empty(numCenterElements, dtype=object)
	centerSiteIndices = nElementsPerUnitCell[:centerSiteElementTypeIndex].sum() + np.arange(numCenterElements)
	centerSiteFractCoords = fractionalUnitCellCoords[centerSiteIndices]
	neighborSiteFractCoords = np.zeros((numCenterElements * numCells, 3))
	supercellFractCoords = np.zeros((numCells * totalElementsPerUnitCell, 3))
	for iCell in range(numCells):
		supercellFractCoords[(iCell * totalElementsPerUnitCell):((iCell + 1) * totalElementsPerUnitCell)] = fractionalUnitCellCoords + unitcellTranslationalCoords[iCell]
	neighborList = np.empty(numCenterElements, dtype=object)
	for centerSiteIndex, centerSiteFractCoord in enumerate(centerSiteFractCoords):
		iNeighborList = []
		for neighborSiteIndex, neighborSiteFractCoord in enumerate(supercellFractCoords):
			if avoidElementType:
				if neighborSiteIndex not in avoidElementIndices:
					latticeDirection = neighborSiteFractCoord - centerSiteFractCoord
					neighborDisplacementVector = np.dot(latticeDirection[None, :], latticeMatrix)
					displacement = np.linalg.norm(neighborDisplacementVector)
					if bridgeCutoffDistLimits[0] < displacement <= bridgeCutoffDistLimits[1]:
						debug = 0
						if debug:
#                                     print np.round(centerSiteFractCoord, 3)
							if np.array_equal(np.round(centerSiteFractCoord / 2, 3), np.array([0.376, 0.299, 0.272])):
								print neighborSiteIndex
								print displacement / ANG2BOHR
						iNeighborList.append(neighborSiteIndex)
		neighborList[centerSiteIndex] = np.asarray(iNeighborList)
	for iCell in range(numCells):
		neighborSiteFractCoords[(iCell * numCenterElements):((iCell + 1) * numCenterElements)] = centerSiteFractCoords + unitcellTranslationalCoords[iCell]
		
	if classList:
		centerSiteClassList = classList[0]
		neighborSiteClassList = np.tile(classList[1], numCells)
		classPairList = np.empty(numCenterElements, dtype=object)

	displacementVectorList = np.empty(numCenterElements, dtype=object)
	latticeDirectionList = np.empty(numCenterElements, dtype=object)
	displacementList = np.empty(numCenterElements, dtype=object)
	numNeighbors = np.zeros(numCenterElements, dtype=int)
	bridgeList = np.empty(numCenterElements, dtype=object)
	for centerSiteIndex, centerSiteFractCoord in enumerate(centerSiteFractCoords):
		iDisplacementVectors = []
		iLatticeDirectionList = []
		iDisplacements = []
		iBridgeList = []
		if classList:
			iClassPairList = []
		for neighborSiteIndex, neighborSiteFractCoord in enumerate(neighborSiteFractCoords):
			latticeDirection = neighborSiteFractCoord - centerSiteFractCoord
			if roundLatticeParameters:
				roundedLatticeDirection = np.round(base * np.round((latticeDirection) / base), prec)
			else:
				roundedLatticeDirection = latticeDirection
			neighborDisplacementVector = np.dot(latticeDirection[None, :], latticeMatrix)
			displacement = np.linalg.norm(neighborDisplacementVector)
			if neighborCutoffDistLimits[0] < displacement <= neighborCutoffDistLimits[1]:
				iDisplacementVectors.append(neighborDisplacementVector)
				iLatticeDirectionList.append(roundedLatticeDirection)
				iDisplacements.append(displacement)
				debug = 0
				if debug:
					debugDistList = [2.96781, 2.99576] #[2.85413, 2.86002, 3.00761, 3.02054]
					debugDist = np.round(displacement / ANG2BOHR, 5)
					if debugDist in debugDistList:
						print debugDist
						print 'center class:', centerSiteClassList[centerSiteIndex]
						print 'neighbor class:', neighborSiteClassList[neighborSiteIndex]
						print 'num of bonds:', len(neighborList[centerSiteIndex])
						print 'center:', np.round(centerSiteFractCoord / 2, 3)
						print 'neighbor:', np.round(neighborSiteFractCoord / 2, 3)
				numNeighbors[centerSiteIndex] += 1
				if classList:
					iClassPairList.append(str(centerSiteClassList[centerSiteIndex]) + ':' + str(neighborSiteClassList[neighborSiteIndex]))
				bridgeSiteExists = 0
				bridgeSiteType = ''
				for iCenterNeighborSEIndex in neighborList[centerSiteIndex]:
					iCenterNeighborFractCoord = supercellFractCoords[iCenterNeighborSEIndex]
					bridgelatticeDirection = neighborSiteFractCoord - iCenterNeighborFractCoord
					bridgeneighborDisplacementVector = np.dot(bridgelatticeDirection[None, :], latticeMatrix)
					bridgedisplacement = np.linalg.norm(bridgeneighborDisplacementVector)
					if bridgeCutoffDistLimits[0] < bridgedisplacement <= bridgeCutoffDistLimits[1]:
						bridgeSiteExists = 1
						bridgeSiteIndex = iCenterNeighborSEIndex
						bridgeSiteQuantumIndices = generateQuantumIndices(localSystemSize, bridgeSiteIndex, nElementsPerUnitCell)
						bridgeSiteType += (', ' if bridgeSiteType != '' else '') + elementTypes[bridgeSiteQuantumIndices[3]]
				if not bridgeSiteExists:
					bridgeSiteType = 'space'
				iBridgeList.append(bridgeSiteType)
		bridgeList[centerSiteIndex] = np.asarray(iBridgeList)
		displacementVectorList[centerSiteIndex] = np.asarray(iDisplacementVectors)
		latticeDirectionList[centerSiteIndex] = np.asarray(iLatticeDirectionList)
		displacementList[centerSiteIndex] = np.asarray(iDisplacements)
		if classList:
			classPairList[centerSiteIndex] = np.asarray(iClassPairList)
	
	from fractions import gcd
	sortedLatticeDirectionList = np.empty(numCenterElements, dtype=object)
	sortedDisplacementList = np.empty(numCenterElements, dtype=object)
	if classList:
		sortedClassPairList = np.empty(numCenterElements, dtype=object)
	sortedBridgeList = np.empty(numCenterElements, dtype=object)
	for iCenterElementIndex in range(numCenterElements):
		sortedDisplacementList[iCenterElementIndex] = displacementList[iCenterElementIndex][displacementList[iCenterElementIndex].argsort()]
		if classList:
			sortedClassPairList[iCenterElementIndex] = classPairList[iCenterElementIndex][displacementList[iCenterElementIndex].argsort()]
		sortedBridgeList[iCenterElementIndex] = bridgeList[iCenterElementIndex][displacementList[iCenterElementIndex].argsort()]
		if roundLatticeParameters:
			latticeDirectionList[iCenterElementIndex] = (latticeDirectionList[iCenterElementIndex] / base).astype(int)
			iCenterLDList = latticeDirectionList[iCenterElementIndex]
			for index in range(numNeighbors[iCenterElementIndex]):
				iAbsCenterLDList = abs(iCenterLDList[index])
				nz = np.nonzero(iAbsCenterLDList)[0]
				nzAbsCenterLDList = iAbsCenterLDList[nz]
				if len(nz) == 1:
					latticeDirectionList[iCenterElementIndex][index] = iCenterLDList[index] / iAbsCenterLDList[nz]
				elif len(nz) == 2:
					latticeDirectionList[iCenterElementIndex][index] = iCenterLDList[index] / gcd(nzAbsCenterLDList[0], nzAbsCenterLDList[1])
				else:
					latticeDirectionList[iCenterElementIndex][index] = iCenterLDList[index] / gcd(gcd(nzAbsCenterLDList[0], nzAbsCenterLDList[1]), nzAbsCenterLDList[2])
		sortedLatticeDirectionList[iCenterElementIndex] = latticeDirectionList[iCenterElementIndex][displacementList[iCenterElementIndex].argsort()]
		np.set_printoptions(suppress=True)

		# print equivalency of all center sites with their respective class reference site
		if printEquivalency:
			refIndex = np.argmax(centerSiteClassList == centerSiteClassList[iCenterElementIndex])
			print np.array_equal(np.round(sortedDisplacementList[refIndex], 4), np.round(sortedDisplacementList[iCenterElementIndex], 4))
		if classList:
			iPathwayList = np.hstack((np.round(sortedLatticeDirectionList[iCenterElementIndex], 4), np.round(sortedDisplacementList[iCenterElementIndex], 4)[:, None], sortedClassPairList[iCenterElementIndex][:, None], sortedBridgeList[iCenterElementIndex][:, None]))
		else:
			iPathwayList = np.hstack((np.round(sortedLatticeDirectionList[iCenterElementIndex], 4), np.round(sortedDisplacementList[iCenterElementIndex], 4)[:, None], sortedBridgeList[iCenterElementIndex][:, None]))
		pathwayList[iCenterElementIndex] = iPathwayList
		if printPathwayList:
			print iPathwayList
# 	import pdb; pdb.set_trace()
	latticeDirectionListFileName = 'latticeDirectionList_' + centerElementType + '-' + neighborElementType + '_cutoff=' + str(neighborCutoff)
	displacementListFileName = 'displacementList_' + centerElementType + '-' + neighborElementType + '_cutoff=' + str(neighborCutoff)
	pathwayFileName = 'pathwayList_' + centerElementType + '-' + neighborElementType + '_cutoff=' + str(neighborCutoff)
	latticeDirectionListFilePath = os.path.join(outdir, latticeDirectionListFileName) + '.npy'
	displacementListFilePath = os.path.join(outdir, displacementListFileName) + '.npy'
	pathwayFilePath = os.path.join(outdir, pathwayFileName) + '.npy'
	np.save(latticeDirectionListFilePath, sortedLatticeDirectionList)
	np.save(displacementListFilePath, sortedDisplacementList)
	np.save(pathwayFilePath, pathwayList)
	return