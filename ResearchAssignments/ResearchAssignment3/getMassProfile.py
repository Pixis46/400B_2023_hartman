import numpy as np
import matplotlib.pyplot as plt
from ReadFile import Read
from CenterOfMass_Template import CenterOfMass
from MassProfile import MassProfile

# This program will be used to spherically average mass of the MW-M31 Merger Remnant
# and fit hernquist profiles to resulting dark matter halos from the final snapshot of data

def sphericalAvg(data, radius, step=0.1): # TODO
    '''
    This function spherically averages the mass over a distribution
    by taking the mass in a speherical shell in steps of radius from the center
    of mass and giving the total mass.
    Parameters:
        data: (array) the distribution of particles to measure
        radius: (float) The radius to measure at in kpc
        step: (float) The step size to use to measure the density
    Returns:
        Mass: (float) The Mass at the given radius in Msun
    '''

def fitHernquist(): # TODO
    '''
    Determines how well-fit to a Hernquist Profile the given data is.

    Parameters:
        data: (array) The density data to compare to a Hernquist Profile
    Returns:
        HernMass: (array) The best-fit Hernquist Profile
        a: (float) The Hernquit scale radius in kpc
        res: (array) The fit residuals
    '''


# Collect Data
MW = CenterOfMass("MW_800.txt")
M31 = CenterOfMass("M31_800.txt")

# TODO: Combine the data together