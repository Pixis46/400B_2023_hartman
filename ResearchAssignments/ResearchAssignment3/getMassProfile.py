import numpy as np
import matplotlib.pyplot as plt
from CenterOfMass import CenterOfMass
from MassProfile import MassProfile
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os

# This program will be used to spherically average mass of the MW-M31 Merger Remnant
# and fit hernquist profiles to resulting dark matter halos from the final snapshot of data
# Combined MW and M31 data by copy-pasting the two files together

def sphericalAvg(data, COM, radius, step=0.1):
    '''
    This function spherically averages the mass over a distribution
    by taking the mass in a speherical shell in steps of radius from the center
    of mass and giving the total mass.
    Parameters:
        data: (array) An array of values representing particle positions in the format
              [x, y, z, m]
        COM: (array) The center of mass position of the data
        radius: (array) The radii to measure at in kpc
        step: (float) The step size to use to measure the density
    Returns:
        Density: (array) The density at the given radius in Msun
    '''
    posVec = np.array([data[0], data[1], data[2]])
    rVals = []
    for i in range(len(data[1])):
        sepfromCOM = np.sqrt(np.sum([i**2 for i in posVec[:, i] - COM]))
        rVals.append(sepfromCOM)
    rVals = np.array(rVals)

    density = []
    for r in radius:
        sliceMin = r
        sliceMax = r + step
        shellIndex = np.where((rVals > sliceMin) & (rVals < sliceMax))[0]
        shellMass = np.sum(data[3][shellIndex])
        print(f"SliceMin: {r}\nSliceMax: {r + step}\nSlice")
        shellVol = (4/3)*np.pi*(sliceMax**3 - sliceMin**3)
        density.append(shellMass / (shellVol))
    return density

def fitHernquist(data): # TODO
    '''
    Determines how well-fit to a Hernquist Profile the given data is.

    Parameters:
        data: (array) The density data to compare to a Hernquist Profile
    Returns:
        HernMass: (array) The best-fit Hernquist Profile
        a: (float) The best-fit Hernquit scale radius in kpc
        res: (array) The fit residuals
    '''
    hernfunc = MassProfile.HernquistMass

combCM = CenterOfMass(os.path.abspath("./MWM31_800.txt"), 1) # The DM halo of the combined system
COM = combCM.COM_P()
positions = np.array([combCM.x, combCM.y, combCM.z, combCM.m])
densityProf = []
radii = np.arange(0, 30.5, 0.1)
densityProf = sphericalAvg(positions, COM, radii)
densityProf = np.array(densityProf)
print(densityProf)
plt.plot(radii, densityProf)