import numpy as np
import matplotlib.pyplot as plt
from Homework7_Template import M33AnalyticOrbit
import os

def makePlot():
    orbit = M33AnalyticOrbit("orbit.txt")
    orbit.OrbitIntegration(0, 0.5, 10)

    diff2vecs = lambda v1, v2: np.sqrt(np.sum([i**2 for i in (v1-v2)]))

    M33orbit = np.genfromtxt("orbit.txt", names=True)
    newM33Pos = np.array([M33orbit['x'], M33orbit['y'], M33orbit['z']])
    newM33velo = np.array([M33orbit['vx'], M33orbit['vy'], M33orbit['vz']])
    # Old Data
    M33Data = np.genfromtxt(os.path.abspath("output/Orbit_M33.txt"), names=True)
    M31Data = np.genfromtxt(os.path.abspath("output/Orbit_M31.txt"), names=True)
    
    posVecM33 = np.array([M33Data['x'], M33Data['y'], M33Data['z']])
    posVecM31 = np.array([M31Data['x'], M31Data['y'], M31Data['z']])
    veloVecM31 = np.array([M31Data['vx'], M31Data['vy'], M31Data['vz']])
    veloVecM33 = np.array([M33Data['vx'], M33Data['vy'], M33Data['vz']])

    differencePM31M33 = []
    for i in range(len(M31Data['t'])):
        diffAtTimet = diff2vecs(np.array(posVecM33[:, i]), np.array(posVecM31[:, i]))
        differencePM31M33.append(diffAtTimet)
    differencePM31M33 = np.array(differencePM31M33)
    differenceVM31M33 = []
    for i in range(len(M31Data['t'])):
        diffAtTimet =  diff2vecs(np.array(veloVecM33[:, i]), np.array(veloVecM31[:, i]))
        differenceVM31M33.append(diffAtTimet)
    differenceVM31M33 = np.array(differenceVM31M33)

    fig, ax = plt.subplots(1, 2, sharex=True)
    trange = np.arange(0, 10, step=0.5)
    ax[0].plot(trange, differencePM31M33, color="k")
    ax[0].plot(trange, newM33Pos, color="cyan")
    ax[1].plot(trange, differenceVM31M33, color="k")
    ax[1].plot(trange, newM33velo, color="cyan")
    
    ax[1].set_xlabel("time (Gyr)")
    ax[0].set_ylabel("Position Difference (kpc)")
    ax[1].set_ylabel("Velocity Difference (kpc/Gyr)")
    plt.show()
makePlot()