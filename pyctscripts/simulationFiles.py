#!/usr/bin/env python

from pathlib import Path
from shutil import copymode

import numpy as np
import yaml


class simulationFiles(object):
    """Class definition to generate simulation files"""

    # constants
    MIN2SEC = 60
    HR2SEC = 60 * MIN2SEC

    def __init__(self, simParamFileName):
        # Load simulation parameters
        with open(simParamFileName, 'r') as stream:
            try:
                params = yaml.load(stream)
            except yaml.YAMLError as exc:
                print(exc)
        for key, value in params.items():
            setattr(self, key, value)

    def slurmFiles(self, srcDirPath, varSpeciesTypeIndex, varSpeciesCountList,
                   kmcPrec=1.00E+04):
        # keywords
        commentKey = '##'
        jobNameKey = ('--job-name="' + self.material + '-'
                      + 'x'.join(str(element) for element in self.systemSize))
        outputKey = ('--output=' + self.material + '-'
                     + 'x'.join(str(element) for element in self.systemSize))
        partitionKey = '--partition='
        timeKey = '--time='

        chargeComb = self.ionChargeType[0] + self.speciesChargeType[0]
        nRuns = len(varSpeciesCountList)
        speciesCountList = [0] * len(self.speciesCount)
        for iRun in range(nRuns):
            # estimate simulation run time in sec
            nonVarSpeciesTypeIndex = int(not varSpeciesTypeIndex)
            speciesCountList[nonVarSpeciesTypeIndex] = self.speciesCount[
                                                        nonVarSpeciesTypeIndex]
            speciesCountList[varSpeciesTypeIndex] = varSpeciesCountList[iRun]

            # determine destination path
            parentDir1 = 'SimulationFiles'
            parentDir2 = ('ionChargeType=' + self.ionChargeType
                          + ';speciesChargeType=' + self.speciesChargeType)
            parentDir3 = (
                    str(speciesCountList[0])
                    + ('electron' if speciesCountList[0] == 1 else 'electrons')
                    + ',' + str(speciesCountList[1])
                    + ('hole' if speciesCountList[1] == 1 else 'holes'))
            parentDir4 = str(self.Temp) + 'K'
            workDir = (('%1.2E' % self.tFinal) + 'SEC,'
                       + ('%1.2E' % self.timeInterval) + 'TimeInterval,'
                       + ('%1.2E' % self.nTraj) + 'Traj')
            systemDirectoryPath = Path.cwd()
            workDirPath = systemDirectoryPath.joinpath(parentDir1, parentDir2,
                                                       parentDir3, parentDir4,
                                                       workDir)
            Path.mkdir(workDirPath, parents=True, exist_ok=True)
            dstFilePath = workDirPath.joinpath(self.dstFileName)

            # generate run file
            srcFilePath = srcDirPath / self.srcFileName
            with srcFilePath.open('r') as srcFile, \
                    dstFilePath.open('w') as dstFile:
                for line in srcFile:
                    if jobNameKey in line and line[:2] != commentKey:
                        line = ('#SBATCH ' + jobNameKey + '_' + chargeComb
                                + '_' + str(varSpeciesCountList[iRun]) + '"\n')
                    elif outputKey in line and line[:2] != commentKey:
                        line = ('#SBATCH ' + outputKey + '_' + chargeComb + '_'
                                + str(varSpeciesCountList[iRun]) + '.out\n')
                    elif partitionKey in line and line[:2] != commentKey:
                        line = ('#SBATCH ' + partitionKey + self.partitionValue
                                + '\n')
                    elif timeKey in line and line[:2] != commentKey:
                        numDays = numHours = numMins = numSec = 0
                        if self.partitionValue is 'debug':
                            numHours = 1
                        elif self.partitionValue is 'mdupuis2':
                            numDays = self.mdSlurmJobMaxTimeLimit
                        else:
                            estRunTime = self.runTime(speciesCountList,
                                                      varSpeciesTypeIndex,
                                                      kmcPrec)
                            if estRunTime > self.gcSlurmJobMaxTimeLimit:
                                numHours = self.gcSlurmJobMaxTimeLimit
                            else:
                                numHours = estRunTime // self.HR2SEC
                                numMins = ((estRunTime // self.MIN2SEC)
                                           % self.MIN2SEC)
                        timeLimit = f'{numDays:02d}-{numHours:02d}:{numMins:02d}:{numSec:02d}'
                        line = '#SBATCH ' + timeKey + timeLimit + '\n'
                    dstFile.write(line)
            copymode(self.srcFileName, str(dstFilePath))

    def runTime(self, speciesCountList, varSpeciesTypeIndex, kmcPrec):
        kTotal = np.dot(self.kTotalPerSpecies, speciesCountList)
        timeStep = 1 / kTotal
        kmcSteps = int(np.ceil(self.tFinal / timeStep / kmcPrec) * kmcPrec)
        numStatesPerStep = np.dot(self.numStatesPerSpecies,
                                  speciesCountList)
        totalStatesPerTraj = numStatesPerStep * kmcSteps
        estRunTime = int(self.timePerState[varSpeciesTypeIndex] * self.nTraj
                         * totalStatesPerTraj
                         + (self.addOnTimeLimit * self.HR2SEC))
        return estRunTime
