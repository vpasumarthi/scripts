#!/usr/bin/env python

from pathlib import Path

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

    def slurmFiles(self, varSpeciesTypeIndex, varSpeciesCountList,
                   kmcPrec=1.00E+04):
        # keywords
        jobNameKey = ('--job-name="' + self.material + '-'
                      + 'x'.join(str(element) for element in self.systemSize))
        outputKey = ('--output=' + self.material + '-'
                     + 'x'.join(str(element) for element in self.systemSize))

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

            # generate slurm file
            with dstFilePath.open('w') as dstFile:
                dstFile.write('#!/bin/sh\n')
                dstFile.write('#SBATCH ' + jobNameKey + '_' + chargeComb + '_'
                              + str(varSpeciesCountList[iRun]) + '"\n')
                dstFile.write('#SBATCH ' + outputKey + '_' + chargeComb + '_'
                              + str(varSpeciesCountList[iRun]) + '.out\n')
                dstFile.write(f'#SBATCH --partition={self.partitionValue}\n')
                if self.partitionValue == 'mdupuis2':
                    dstFile.write('#SBATCH --clusters=chemistry\n')
                numDays = numHours = numMins = numSec = 0
                if self.partitionValue == 'debug':
                    numHours = 1
                elif self.partitionValue == 'mdupuis2':
                    numDays = self.mdSlurmJobMaxTimeLimit
                else:
                    estRunTime = self.runTime(speciesCountList,
                                              varSpeciesTypeIndex, kmcPrec)
                    if estRunTime > self.gcSlurmJobMaxTimeLimit:
                        numHours = self.gcSlurmJobMaxTimeLimit
                    else:
                        numHours = estRunTime // self.HR2SEC
                        numMins = (estRunTime // self.MIN2SEC) % self.MIN2SEC
                timeLimit = f'{numDays:02d}-{numHours:02d}:{numMins:02d}:{numSec:02d}'
                dstFile.write(f'#SBATCH --time={timeLimit}\n')
                dstFile.write(f'#SBATCH --nodes={self.numNodes}\n')
                dstFile.write(f'#SBATCH --tasks-per-node={self.numTasksPerNode}\n')
                if self.exclusive:
                    dstFile.write('#SBATCH --exclusive\n')
                if self.mem:
                    dstFile.write(f'#SBATCH --mem={self.mem}\n')
                dstFile.write(
                    "\n# Job description:\n"
                    "# run KMC simulation followed by performing MSD analysis"
                    " of the output trajectories\n\n"
                    "echo \"SLURM_JOBID=\"$SLURM_JOBID\n"
                    "echo \"SLURM_JOB_NODELIST\"=$SLURM_JOB_NODELIST\n"
                    "echo \"SLURM_NNODES\"=$SLURM_NNODES\n"
                    "echo \"SLURMTMPDIR=\"$SLURMTMPDIR\n\n"
                    "echo \"working directory = \"$SLURM_SUBMIT_DIR\n\n"
                    "module load python\n"
                    "source activate py36\n"
                    "module list\n"
                    "ulimit -s unlimited\n\n"
                    "# The initial srun will trigger the SLURM prologue on"
                    " the compute nodes.\n"
                    "NPROCS=`srun --nodes=${SLURM_NNODES}"
                    " bash -c 'hostname' |wc -l`\n"
                    "echo NPROCS=$NPROCS\n"
                    "echo \"Launch mymodel with srun\"\n\n"
                    "#The PMI library is necessary for srun\n"
                    "export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so\n")
                if self.submitRun:
                    dstFile.write("srun Run.py\n")
                if self.submitMSD:
                    dstFile.write("srun MSD.py\n")
                dstFile.write("\necho \"All Done!\"\n")
        return None

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
