import numpy as np
import matplotlib.pyplot as plt
from CenterOfMass_Template import CenterOfMass
from MassProfile import MassProfile
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# This program will be used to spherically average mass of the MW-M31 Merger Remnant
# and fit hernquist profiles to resulting dark matter halos from the final snapshot of data
# Combined MW and M31 data by copy-pasting the two files together

def sphericalAvg(dataFile, radius, step=0.1):
    '''
    This function spherically averages the mass over a distribution
    by taking the mass in a speherical shell in steps of radius from the center
    of mass and giving the total mass.
    Parameters:
        dataFile: (string) the file holding the particle data to measure
        radius: (float) The radius to measure at in kpc
        step: (float) The step size to use to measure the density
    Returns:
        Mass: (float) The Mass at the given radius in Msun
    '''
    combCM = CenterOfMass(dataFile, 1) # The DM halo of the combined system
    COM = combCM.COM_P()
    posVec = np.array([combCM['x'], combCM['y'], combCM['z']])
    rVals = []
    for i in range(len(combCM['x'])):
        sepfromCOM = np.sqrt(np.sum([i**2 for i in posVec[:, i] - COM]))
        rVals.append(sepfromCOM)
    rVals = np.array(rVals)

    sliceMin = radius
    sliceMax = radius + step
    shellIndex = np.where((rVals > sliceMin) & (rVals < sliceMax))[0]
    shellMass = np.sum(combCM.m[shellIndex])
    shellVol = (4/3)*np.pi*(sliceMax**3 - sliceMin**3)
    density = shellMass / (shellVol)

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

densityProf = []
radii = np.arange(0, 30.5, 0.1)
for r in radii:
    densityProf.append(sphericalAvg("MWM31_800.txt", r))
densityProf = np.array(densityProf)
plt.plot(radii, densityProf)