#!/usr/bin/python

# Source: https://github.com/minillinim/ellipsoid/blob/master/ellipsoid.py

from __future__ import division
import matplotlib.pyplot as plt
import importlib
import numpy as np
from numpy import linalg


class EllipsoidTool:
    """Some stuff for playing with ellipsoids"""
    def __init__(self): pass

    def getMinVolEllipse(self, P=None, tolerance=0.01):
        """ Find the minimum volume ellipsoid which holds all the points

        Based on work by Nima Moshtagh
        http://www.mathworks.com/matlabcentral/fileexchange/9542
        and also by looking at:
        http://cctbx.sourceforge.net/current/python/scitbx.math.minimum_covering_ellipsoid.html
        Which is based on the first reference anyway!

        Here, P is a numpy array of N dimensional points like this:
        P = [[x,y,z,...], <-- one point per line
             [x,y,z,...],
             [x,y,z,...]]

        Returns:
        (center, radii, rotation)

        """
        (N, d) = np.shape(P)
        d = float(d)

        # Q will be our working array
        Q = np.vstack([np.copy(P.T), np.ones(N)])
        QT = Q.T

        # initializations
        err = 1.0 + tolerance
        u = (1.0 / N) * np.ones(N)

        # Khachiyan Algorithm
        while err > tolerance:
            V = np.dot(Q, np.dot(np.diag(u), QT))
            # M the diagonal vector of an NxN matrix
            M = np.diag(np.dot(QT, np.dot(linalg.inv(V), Q)))
            j = np.argmax(M)
            maximum = M[j]
            step_size = (maximum - d - 1.0) / ((d + 1.0) * (maximum - 1.0))
            new_u = (1.0 - step_size) * u
            new_u[j] += step_size
            err = np.linalg.norm(new_u - u)
            u = new_u

        # center of the ellipse
        center = np.dot(P.T, u)

        # the A matrix for the ellipse
        A = linalg.inv(
                       np.dot(P.T, np.dot(np.diag(u), P)) -
                       np.array([[a * b for b in center] for a in center])
                       ) / d

        # Get the values we'd like to return
        _, s, rotation = linalg.svd(A)
        radii = 1.0/np.sqrt(s)

        return (center, radii, rotation)

    def getEllipsoidVolume(self, radii):
        """Calculate the volume of the blob"""
        return 4./3.*np.pi*radii[0]*radii[1]*radii[2]

    def plotEllipsoid(self, center, radii, rotation, ax=None, plotAxes=False,
                      cageColor='b', cageAlpha=0.2):
        """Plot an ellipsoid"""
        make_ax = ax is None
        if make_ax:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        u = np.linspace(0.0, 2.0 * np.pi, 100)
        v = np.linspace(0.0, np.pi, 100)

        # cartesian coordinates that correspond to the spherical angles:
        x = radii[0] * np.outer(np.cos(u), np.sin(v))
        y = radii[1] * np.outer(np.sin(u), np.sin(v))
        z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
        # rotate accordingly
        for i in range(len(x)):
            for j in range(len(x)):
                [x[i, j], y[i, j], z[i, j]] = (
                                np.dot([x[i, j], y[i, j], z[i, j]], rotation)
                                + center)

        if plotAxes:
            # make some purdy axes
            axes = np.array([[radii[0], 0.0, 0.0],
                             [0.0, radii[1], 0.0],
                             [0.0, 0.0, radii[2]]])
            # rotate accordingly
            for i in range(len(axes)):
                axes[i] = np.dot(axes[i], rotation)

            # plot axes
            for p in axes:
                X3 = np.linspace(-p[0], p[0], 100) + center[0]
                Y3 = np.linspace(-p[1], p[1], 100) + center[1]
                Z3 = np.linspace(-p[2], p[2], 100) + center[2]
                ax.plot(X3, Y3, Z3, color=cageColor)

        # plot ellipsoid
        ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color=cageColor,
                          alpha=cageAlpha)
        ax.set_xlim(xmin=-1000, xmax=1000)
        ax.set_ylim(ymin=-1000, ymax=1000)
        ax.set_zlim(zmin=-1000, zmax=1000)
        if make_ax:
            plt.show()
            plt.close(fig)
            del fig


if __name__ == "__main__":
    # make 100 random points
    # P = np.reshape([random()*100 for i in range(300)],(100,3))
    fileName = 'unwrappedTraj.dat'
    numTrajRecorded = int(1.00E+02)
    tFinal = 1.00E-04
    timeInterval = 1.00E-08
    positionArray = np.loadtxt('unwrappedTraj.dat')
    numPathStepsPerTraj = int(tFinal / timeInterval) + 1
    nSpecies = int(positionArray.shape[1] / 3)
    dataArray = np.zeros((numPathStepsPerTraj, numTrajRecorded * nSpecies, 3))

    for trajIndex in range(numTrajRecorded):
        headStart = trajIndex * numPathStepsPerTraj
        for step in range(numPathStepsPerTraj):
            stepPosition = positionArray[headStart + step]
            for speciesIndex in range(nSpecies):
                dataArray[step, trajIndex * nSpecies + speciesIndex, :] = (
                        stepPosition[speciesIndex * 3: (speciesIndex + 1) * 3])
    plt.switch_backend('Agg')
    importlib.import_module('mpl_toolkits.mplot3d').Axes3D
    # find the ellipsoid
    # import pdb; pdb.set_trace()
    index = 0
    for step in range(numPathStepsPerTraj):
        if step % 100 == 0 and step != 0:
            P = dataArray[step]
            ET = EllipsoidTool()
            (center, radii, rotation) = ET.getMinVolEllipse(P, .01)

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            # plot points
            ax.scatter(P[:, 0], P[:, 1], P[:, 2], color='g', marker='*', s=100)

            # plot ellipsoid
            ET.plotEllipsoid(center, radii, rotation, ax=ax, plotAxes=True)

            plt.show()
            plt.savefig('ellipsoid_new' + str(index) + '.png')
            plt.close(fig)
            del fig
            index += 1
