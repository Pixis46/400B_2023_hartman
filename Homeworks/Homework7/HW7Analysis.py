import numpy as np
import matplotlib.pyplot as plt
from Homework7_Template import M33AnalyticOrbit
import os

def makePlot():
    orbit = M33AnalyticOrbit("orbit.txt")
    step = 0.1
    orbit.OrbitIntegration(0, 0.1, 10)
    
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

    trange = np.arange(0, 10, step=step)
    newM33PosMag = []
    for i in range(len(trange)):
        newM33PosMag.append(np.sqrt(np.sum(np.array([i**2 for i in newM33Pos[:, i]]))))
    newM33PosMag = np.array(newM33PosMag)
    newM33Vmag = []
    for i in range(len(trange)):
        newM33Vmag.append(np.sqrt(np.sum(np.array([i**2 for i in newM33velo[:, i]]))))
    newM33Vmag = np.array(newM33Vmag)

    fig, ax = plt.subplots(2, 1, sharex=True)
    ax[0].plot(M31Data['t'], differencePM31M33, color="k", label="Simulation")
    ax[0].plot(trange, newM33PosMag, color="cyan", label="Integrator")
    ax[1].plot(M31Data['t'], differenceVM31M33, color="k")
    ax[1].plot(trange, newM33Vmag, color="cyan")
    
    fig.suptitle("Simulation Data vs. Orbit Integrator")
    ax[1].set_xlabel("time (Gyr)")
    ax[0].set_ylabel("Position Difference (kpc)")
    ax[1].set_ylabel("Velocity Difference (kpc/Gyr)")
    ax[0].legend(loc="upper right")
    plt.savefig("HW7Figures.png")
makePlot()