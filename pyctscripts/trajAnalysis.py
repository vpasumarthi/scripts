#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def trajAnalysis(dstPath, dispPrec):
    # CONSTANTS
    EPSILON0 = 8.854187817E-12  # Electric constant in F.m-1
    ANG = 1E-10  # Angstrom in m

    # FUNDAMENTAL ATOMIC UNITS
    # Source: http://physics.nist.gov/cuu/Constants/Table/allascii.txt
    EMASS = 9.10938356E-31  # Electron mass in Kg
    ECHARGE = 1.6021766208E-19  # Elementary charge in C
    HBAR = 1.054571800E-34  # Reduced Planck's constant in J.sec
    KE = 1 / (4 * np.pi * EPSILON0)

    # DERIVED ATOMIC UNITS
    # Bohr radius in m
    BOHR = HBAR**2 / (EMASS * ECHARGE**2 * KE)

    # CONVERSIONS
    ANG2BOHR = ANG / BOHR

    positionArray = (np.loadtxt(dstPath.joinpath('unwrappedTraj.dat'))
                     / ANG2BOHR)
    numPositions = len(positionArray)
    desired_indices = [0]
    for stepIndex in range(1, numPositions):
        if not np.array_equal(
                            np.round(positionArray[stepIndex, :], dispPrec),
                            np.round(positionArray[stepIndex-1, :], dispPrec)):
            desired_indices.append(stepIndex)
    positionArray = positionArray[desired_indices]

    dispVecArray = np.diff(positionArray, axis=0)
    dispArray = np.linalg.norm(dispVecArray, axis=1)

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
