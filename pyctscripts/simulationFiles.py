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

    def dstPath(self, speciesCountList):
        # determine destination path
        childDir1 = 'SimulationFiles'
        childDir2 = ('ionChargeType=' + self.system['ionChargeType']
                     + ';speciesChargeType='
                     + self.system['speciesChargeType'])
        childDir3 = (
                    str(speciesCountList[0])
                    + ('electron' if speciesCountList[0] == 1 else 'electrons')
                    + ',' + str(speciesCountList[1])
                    + ('hole' if speciesCountList[1] == 1 else 'holes'))
        childDir4 = str(self.system['Temp']) + 'K'
        workDir = (('%1.2E' % self.run['tFinal']) + 'SEC,'
                   + ('%1.2E' % self.run['timeInterval']) + 'TimeInterval,'
                   + ('%1.2E' % self.run['nTraj']) + 'Traj')
        systemDirectoryPath = Path.cwd()
        workDirPath = (systemDirectoryPath / childDir1 / childDir2 / childDir3
                       / childDir4 / workDir)
        workDirDepth = len(workDirPath.parts) - len(systemDirectoryPath.parts)
        return (workDirPath, workDirDepth)

    def simParmFiles(self, varSpeciesTypeIndex, varSpeciesCountList):
        nRuns = len(varSpeciesCountList)
        speciesCountList = [0] * len(self.system['speciesCount'])
        for iRun in range(nRuns):
            nonVarSpeciesTypeIndex = int(not varSpeciesTypeIndex)
            speciesCountList[nonVarSpeciesTypeIndex] = (
                        self.system['speciesCount'][nonVarSpeciesTypeIndex])
            speciesCountList[varSpeciesTypeIndex] = varSpeciesCountList[iRun]

            (workDirPath, workDirDepth) = self.dstPath(speciesCountList)
            self.system['workDirDepth'] = workDirDepth
            Path.mkdir(workDirPath, parents=True, exist_ok=True)
            dstFilePath = workDirPath.joinpath(self.system['dstFileName'])

            # generate simulation parameter file
            with dstFilePath.open('w') as dstFile:
                dstFile.write('# System parameters:\n')
                yaml.dump(self.system, dstFile)
                dstFile.write('\n')
                dstFile.write('# Run parameters:\n')
                yaml.dump(self.run, dstFile, default_flow_style=False)
                dstFile.write('\n')
                dstFile.write('# MSD parameters:\n')
                yaml.dump(self.msd, dstFile, default_flow_style=False)
        return None

    def runTime(self, speciesCountList, varSpeciesTypeIndex, kmcPrec):
        kTotal = np.dot(self.kTotalPerSpecies, speciesCountList)
        timeStep = 1 / kTotal
        kmcSteps = int(np.ceil(self.run['tFinal'] / timeStep / kmcPrec)
                       * kmcPrec)
        numStatesPerStep = np.dot(self.numStatesPerSpecies,
                                  speciesCountList)
        totalStatesPerTraj = numStatesPerStep * kmcSteps
        estRunTime = int(self.timePerState[varSpeciesTypeIndex]
                         * self.run['nTraj'] * totalStatesPerTraj
                         + (self.addOnTimeLimit * self.HR2SEC))
        return estRunTime

    def slurmFiles(self, varSpeciesTypeIndex, varSpeciesCountList,
                   kmcPrec=1.00E+04):
        # keywords
        jobNameKey = ('--job-name="' + self.system['material'] + '-'
                      + 'x'.join(str(element)
                                 for element in self.system['systemSize']))
        outputKey = ('--output=' + self.system['material'] + '-'
                     + 'x'.join(str(element)
                                for element in self.system['systemSize']))

        chargeComb = (self.system['ionChargeType'][0]
                      + self.system['speciesChargeType'][0])
        nRuns = len(varSpeciesCountList)
        speciesCountList = [0] * len(self.system['speciesCount'])
        for iRun in range(nRuns):
            # estimate simulation run time in sec
            nonVarSpeciesTypeIndex = int(not varSpeciesTypeIndex)
            speciesCountList[nonVarSpeciesTypeIndex] = (
                        self.system['speciesCount'][nonVarSpeciesTypeIndex])
            speciesCountList[varSpeciesTypeIndex] = varSpeciesCountList[iRun]
            (workDirPath, _) = self.dstPath(speciesCountList)
            Path.mkdir(workDirPath, parents=True, exist_ok=True)
            dstFilePath = workDirPath.joinpath(self.slurm['dstFileName'])

            # generate slurm file
            with dstFilePath.open('w') as dstFile:
                dstFile.write('#!/bin/sh\n')
                dstFile.write('#SBATCH ' + jobNameKey + '_' + chargeComb + '_'
                              + str(varSpeciesCountList[iRun]) + '"\n')
                dstFile.write('#SBATCH ' + outputKey + '_' + chargeComb + '_'
                              + str(varSpeciesCountList[iRun]) + '.out\n')
                dstFile.write(f"#SBATCH --partition={self.slurm['partitionValue']}\n")
                if self.slurm['partitionValue'] == 'mdupuis2':
                    dstFile.write('#SBATCH --clusters=chemistry\n')
                numDays = numHours = numMins = numSec = 0
                if self.slurm['partitionValue'] == 'debug':
                    numHours = 1
                elif self.slurm['partitionValue'] == 'mdupuis2':
                    numDays = self.slurm['mdSlurmJobMaxTimeLimit']
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
                dstFile.write(f"#SBATCH --nodes={self.slurm['numNodes']}\n")
                dstFile.write(f"#SBATCH --tasks-per-node={self.slurm['numTasksPerNode']}\n")
                if self.slurm['exclusive']:
                    dstFile.write('#SBATCH --exclusive\n')
                if self.slurm['mem']:
                    dstFile.write(f"#SBATCH --mem={self.slurm['mem']}\n\n")
                if self.slurm['email']:
                    dstFile.write(f"#SBATCH --mail-user={self.slurm['email']}\n")
                    dstFile.write("#SBATCH --mail-type=END\n")
                dstFile.write(
                    "#SBATCH --constraint=IB\n"
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
                if self.slurm['submitRun']:
                    dstFile.write("srun Run.py\n")
                if self.slurm['submitMSD']:
                    dstFile.write("srun MSD.py\n")
                dstFile.write("\necho \"All Done!\"\n")
        return None

    def runFiles(self, varSpeciesTypeIndex, varSpeciesCountList):
        nRuns = len(varSpeciesCountList)
        speciesCountList = [0] * len(self.system['speciesCount'])
        for iRun in range(nRuns):
            nonVarSpeciesTypeIndex = int(not varSpeciesTypeIndex)
            speciesCountList[nonVarSpeciesTypeIndex] = (
                        self.system['speciesCount'][nonVarSpeciesTypeIndex])
            speciesCountList[varSpeciesTypeIndex] = varSpeciesCountList[iRun]

            (workDirPath, _) = self.dstPath(speciesCountList)
            dstFilePath = workDirPath.joinpath(self.run['dstFileName'])

            # generate simulation parameter file
            with dstFilePath.open('w') as dstFile:
                dstFile.write(
                    "#!/usr/bin/env python\n\n"
                    "from pathlib import Path\n\n"
                    "from PyCT.materialRun import materialRun\n\n"
                    "dstPath = Path.cwd()\n"
                    "materialRun(dstPath)\n")
        return None
