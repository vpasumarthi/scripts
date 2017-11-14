#!/usr/bin/env python

import matplotlib.pyplot as plt
import importlib

import numpy as np

fileName = 'unwrappedTraj.dat'
numTrajRecorded = int(1.00E+02)
tFinal = 1.00E-04
timeInterval = 1.00E-08
positionArray = np.loadtxt('unwrappedTraj.dat')
numPathStepsPerTraj = int(tFinal / timeInterval) + 1
nSpecies = int(positionArray.shape[1] / 3)

plt.switch_backend('Agg')
importlib.import_module('mpl_toolkits.mplot3d').Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for trajIndex in range(1):
    headStart = trajIndex * numPathStepsPerTraj
    for step in range(numPathStepsPerTraj):
        inertiaTensor = np.zeros((3, 3))
        stepPosition = positionArray[headStart + step]
        for speciesIndex in range(nSpecies):
            inertiaTensor[0][0] += (stepPosition[speciesIndex * 3 + 1]**2
                                    + stepPosition[speciesIndex * 3 + 2]**2)
            inertiaTensor[1][1] += (stepPosition[speciesIndex * 3 + 0]**2
                                    + stepPosition[speciesIndex * 3 + 2]**2)
            inertiaTensor[2][2] += (stepPosition[speciesIndex * 3 + 0]**2
                                    + stepPosition[speciesIndex * 3 + 1]**2)
            inertiaTensor[0][1] -= (stepPosition[speciesIndex * 3 + 0]
                                    * stepPosition[speciesIndex * 3 + 1])
            inertiaTensor[1][0] -= (stepPosition[speciesIndex * 3 + 0]
                                    * stepPosition[speciesIndex * 3 + 1])
            inertiaTensor[1][2] -= (stepPosition[speciesIndex * 3 + 1]
                                    * stepPosition[speciesIndex * 3 + 2])
            inertiaTensor[2][1] -= (stepPosition[speciesIndex * 3 + 1]
                                    * stepPosition[speciesIndex * 3 + 2])
            inertiaTensor[0][2] -= (stepPosition[speciesIndex * 3 + 0]
                                    * stepPosition[speciesIndex * 3 + 2])
            inertiaTensor[2][0] -= (stepPosition[speciesIndex * 3 + 0]
                                    * stepPosition[speciesIndex * 3 + 2])
        eigVal, eigVec = np.linalg.eig(inertiaTensor)
        print(eigVal)
        absVec = np.multiply(eigVec, eigVal[:, np.newaxis])
#         print(absVec)
#         print(absVec[0])
#         import pdb; pdb.set_trace()
#         ax.plot_surface(absVec[0], absVec[1], absVec[2])
#         plt.show()
# print(absVec)
# print(absVec[0])
import pdb; pdb.set_trace()
startPos = np.zeros(3)
posStack0 = np.vstack((startPos, absVec[0]))
ax.plot(posStack0[:, 0], posStack0[:, 1], posStack0[:, 2])
posStack1 = np.vstack((startPos, absVec[1]))
ax.plot(posStack1[:, 0], posStack1[:, 1], posStack1[:, 2])
posStack2 = np.vstack((startPos, absVec[2]))
ax.plot(posStack2[:, 0], posStack2[:, 1], posStack2[:, 2])
ax.plot_surface(absVec[0], absVec[1], absVec[2],
                rstride=4, cstride=4, color='b')
# plt.show()
plt.savefig('ellipsoid.png')
