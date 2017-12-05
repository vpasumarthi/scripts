#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from PyCT import constants


def trajAnalysis(dstPath, dispPrec):
    positionArray = (np.loadtxt(dstPath.joinpath('unwrappedTraj.dat'))
                     / constants.ANG2BOHR)
    numPositions = len(positionArray)
    desired_indices = [0]
    for stepIndex in range(1, numPositions):
        if not np.array_equal(
                            np.round(positionArray[stepIndex, :], dispPrec),
                            np.round(positionArray[stepIndex-1, :], dispPrec)):
            desired_indices.append(stepIndex)
    newPositionArray = np.copy(positionArray[desired_indices])

    dispVecArray = np.diff(newPositionArray, axis=0)
    dispArray = np.linalg.norm(dispVecArray, axis=1)

    numSteps = len(newPositionArray) - 1
    rattleList = []
    numRattles = 1
    for stepIndex in range(2, numSteps+1):
        if np.array_equal(np.round(newPositionArray[stepIndex, :], dispPrec),
                          np.round(newPositionArray[stepIndex-2, :], dispPrec)):
            numRattles += 1
        else:
            if numRattles > 2:
                rattleDist = np.linalg.norm(newPositionArray[stepIndex-1, :]
                                            - newPositionArray[stepIndex-2, :])
                escapeDist = np.linalg.norm(newPositionArray[stepIndex, :]
                                            - newPositionArray[stepIndex-1, :])
                rattleList.append([np.round(rattleDist, dispPrec),
                                   numRattles,
                                   escapeDist])
                numRattles = 1
    rattleArray = np.asarray(rattleList)

    # round displacements to given precision
    dispArray = np.round(dispArray, dispPrec)

    # collect to bins
    [unique, counts] = np.unique(dispArray, return_counts=True)

    plt.switch_backend('Agg')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    proc_indices = np.arange(len(unique))
    xtick_items = ['%1.4f' % item for item in unique]
    plt.bar(proc_indices, counts, align='center', alpha=0.5, edgecolor='black')
    plt.xticks(proc_indices, xtick_items, rotation='vertical')

    for i, v in enumerate(counts):
        ax.text(i - 0.2, v + 100, str(v), color='green', rotation='vertical',
                fontweight='bold')
    addRectangle = 1
    if addRectangle:
        ax.add_patch(patches.Rectangle((13.5, -10), 4, 500, fill=False,
                                       color='red'))
    ax.set_xlabel('Hopping Distance')
    ax.set_ylabel('Counts')
    ax.set_title('Histogram of processes')
    filename = 'process_histogram'
    figureName = filename + '.png'
    figurePath = dstPath / figureName
    plt.tight_layout()
    plt.savefig(str(figurePath))
    return None
