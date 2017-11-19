#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


def plotEnergyDiagram(inputDataFilePath, colorInfo, markerInfo, linestyleInfo,
                      markerSize, fontSize):
    plt.switch_backend('Agg')
    energyData = np.loadtxt(inputDataFilePath)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    numData = len(energyData)
    rcData = np.linspace(0, 1, numData)
    ax.plot(rcData[0:int(numData / 2) + 1],
            energyData[0:int(numData / 2) + 1], c=colorInfo,
            marker=markerInfo, ls=linestyleInfo, markersize=markerSize)
    ax.plot(rcData[int(numData / 2):numData],
            energyData[int(numData / 2):numData], c=colorInfo,
            marker=markerInfo, ls=linestyleInfo, markersize=markerSize)
    ax.set_xlabel('Reaction Coordinate', fontsize=fontSize)
    ax.set_ylabel('Energy (eV)', fontsize=fontSize)
    plt.xticks(fontsize=fontSize)
    plt.yticks(fontsize=fontSize)
    plt.show()
    figurePath = inputDataFilePath.with_suffix('.png')
    plt.savefig(str(figurePath))
