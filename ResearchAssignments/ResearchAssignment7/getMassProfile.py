import numpy as np
import matplotlib.pyplot as plt
from CenterOfMass_Template import CenterOfMass
from MassProfile import MassProfile
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from ReadFile import Read
import os

# This program will be used to spherically average mass of the MW-M31 Merger Remnant
# and fit hernquist profiles to resulting dark matter halos from three snapshots of data: 300, 450, 800


# This code is designed to answer the question of how dark matter halo profiles (both density and mass)
# change with time throughout a merger.


def concatFiles(filename1, filename2, snapNum): # NEW
    '''
    combines the data from the two files into a new one with the header format as follows:
    Time (time)
    Total (total combined particle value)
    mass in 1e10, x,y,z in kpc, and vx,vy,vz in km/s
    The names of the data columns are listed on the next line
    Data

    These two files are assumed to be from the same SnapNumber and in the
    format above
    
    Parameters:
        filename1: (str) The path to the first file
        filename2: (str) The path to the second file
        snapNum: (int) The SnapNumber of the files. Both are assumed to be the same
    Returns:
        None
    '''
    time1, num1, data1 = Read(filename1)
    time2, num2, data2 = Read(filename2)
    concatData = np.concatenate([data1, data2])
    headers = ", ".join(data1.dtype.names)
    path = os.path.abspath(f"./ResearchAssignments/ResearchAssignment3/MWM31_{snapNum}.txt")
    np.savetxt(path, concatData, fmt="%1.9f", comments="", header=f"Time\t{time1.value}\nTotal\t{num1+num2}\nmass in 1e10, x,y,z in kpc and vx,vy,vz in km/s\n#{headers}")

def HernquistMass(r, a, Mhalo=1e10):
        '''
        Computes the mass enclosed within a given radius using a theoretical
        Hernquist (1990) Profile:
        density = (Ma)/(2*pi*r*(r+a)^3)
        M(r) = (M_halo*r^2)/((a+r)^2) 
        with r as the radius, a as the scale height, and Mhalo as the mass
        of the dark matter halo

        Parameters:
            r: (float or array of floats) The radius to compute the enclosed
               Mass at.
            a: (float) The scale height for the Hernquist profile
            Mhalo: (float) The mass of the dark matter halo of the galaxy
        
        Returns:
            EnclosedM: (float or array of floats)
                       the total enclosed mass at the given radii
                       in Solar Masses.
        '''
        M = (Mhalo*(r**2))/((a+r)**2)# Hernquist M
        return M

def HernquistDensity(r, a, Mhalo=1e10):
    '''
    Computes the mass density at a given radius using the Hernquist (1990)
    Density Profile:
    (Mhalo*a) / (2*np.pi*r*((r+a)**3))
    With r as the radius, a as the Hernquist scale height, and Mhalo as
    the total Dark Matter Halo Mass.

    Parameters:
        r: (array of floats) The radii to compute density at in kpc
        a: (float) The Hernquist scale height in kpc
        Mhalo: (float) The mass of the dark matter halo in Msun
    Returns:
        rho: (array of floats) The density at the given radii in Msun/kpc^3
    '''
    rho = (Mhalo*a) / (2*np.pi*r*((r+a)**3))
    return rho

def sphericalAvg(data, COM, radius): # NOT NEW, turns out this is just the Hernquist Mass code oops
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

    mass = []
    for r in radius:
        shellIndex = np.where((rVals < r))[0]
        shellMass = np.sum(data[3][shellIndex])
        #print(f"SliceMin: {r}\nSliceMax: {r + step}\nIndex: {shellIndex}\nSliceMass: {shellMass}")
        mass.append(shellMass)
    return np.array(mass)

def densProfile(data, COM, radius, step=0.1): # NEW
    '''
    This function determines the density profile of given data
    via spherical averaging. It takes an array of radii and, for each
    radius, sums the mass of all particles within a spherical shell
    defined by radius r+step, in kpc, on the outside, and r on the inside.
    The total mass is divided by the volume of the spherical shell to recover
    the density at the radius given.

    Parameters:
        data: (array of arrays of floats) the particles and associated position and mass data
              in an array. Positions in kpc and masses in 1e10 Msun
        COM: (array of floats) The position of the Center of Mass of the combined
             Milky Way-M31 system, in units of kpc
        radius: (array of floats) The radii to compute the density at, in kpc
        step: (float) The step to use to find the spherical shell to compute
              density in kpc. Default is 0.1 kpc
    Returns:
        density: (array of floats) An array of density values for the given radii
    '''
    posVec = np.array([data[0], data[1], data[2]])
    rVals = []
    for i in range(len(data[1])):
        # Compute vector magnitude of particle positions
        sepfromCOM = np.sqrt(np.sum([i**2 for i in posVec[:, i] - COM]))
        rVals.append(sepfromCOM)
    rVals = np.array(rVals)
    density = []
    for r in radius:
        sliceMax = r + step # define max shell radius
        shellIndex = np.where((rVals > r) & (rVals < sliceMax))[0]
        shellMass = np.sum(data[3][shellIndex])
        shellVol = (4/3) * np.pi*(sliceMax**3 - r**3)
        density.append(shellMass/shellVol)
    return np.array(density)


def fitHernquist(toFit, radii, haloMass, dens=False): # NEW
    '''
    Determines how well-fit to a Hernquist Profile the given data is.

    Parameters:
        toFit: (array) The Mass data to compare to a Hernquist Profile
        radii: (array) The radii that the toFit data were collected at
        Mhalo: (float) total halo mass in the galaxy
    Returns:
        HernMass: (array) The best-fit Hernquist Profile
        a: (float) The best-fit Hernquit scale radius in kpc
        res: (array) The fit residuals
    '''
    if not dens:
        hernfunc = lambda r, a:HernquistMass(r, a, Mhalo=haloMass)
        popt, pcov = curve_fit(hernfunc, radii, toFit)
        print(popt, pcov)
        return HernquistMass(radii, popt[0], haloMass), np.around(popt[0], 2), pcov[0][0]
    else:
        hernfunc = lambda r, a: HernquistDensity(r, a, Mhalo=haloMass)
        popt, pcov = curve_fit(hernfunc, radii, toFit)
        print(popt, pcov)
        return HernquistDensity(radii, popt[0], haloMass), np.around(popt[0], 2), pcov[0][0]
    
def pltSnapNum(combCM, COM, snapNum): # NEW
    '''
    This function creates the plot of Mass enclosed and density versus radius.
    The plots are two panel plots with the density distribution on the right and
    mass on the left, with the theoretical Hernquist profiles plotted and the fit
    parameters listed.
    Parameters:
        combCM: (CenterOfMass object) The CenterOfMass object that contains data on
                particle positions and masses. Positions in kpc and masses in 1e10 Msun
        COM: (array of floats) The center of mass of the system, computed with disk stars.
             positions given in kpc.
        snapNum: (int) The SnapNumber of the data. This is representative of time, with
                 time=SnapNum*10/0.7
    Returns:
        None
    '''
    positions = np.array([combCM.x, combCM.y, combCM.z, 1e10*combCM.m])
    radii = np.arange(0, 30, 0.1)
    massProf = sphericalAvg(positions, COM, radii)
    densProf = densProfile(positions, COM, radii, step = 1.220) # Radius

    Mhalo = np.sum(combCM.m)*1e10 # Msun
    hernFit, fita, pcov = fitHernquist(massProf, radii, Mhalo, dens=False)
    #densFit, fitadens, pcovdens = fitHernquist(densProf[8:], radii[8:], Mhalo, dens=True)
    densFit = HernquistDensity(radii[8:], fita, Mhalo)

    fig, ax = plt.subplots(1, 2, sharex=True, figsize=(11, 6))
    ax[0].semilogy(radii, massProf, color="k", label="Mass Data")
    ax[0].semilogy(radii, hernFit, color="g", ls="--", label="Hernquist Mass")
    ax[0].set_xlabel("r (kpc)")
    ax[0].set_ylabel(r"$M_{enclosed}$ ($M_{\odot}$)")
    ax[0].text(0.5, 0.45, "Merger Remnant", transform=ax[0].transAxes, size=13)
    ax[0].text(0.5, 0.41, r"a = " + f"{fita} $\pm$ {np.around(3*np.sqrt(pcov), 2)} kpc", transform=ax[0].transAxes, size=11)
    ax[0].text(0.5, 0.37, r"t = " + f"{np.around(snapNum*10/0.7, 0)} Myr", transform=ax[0].transAxes, size=11)
    ax[0].legend(loc="lower right")

    ax[1].semilogy(radii, densProf, color="k", label="Density Data")
    ax[1].semilogy(radii[8:], densFit, color="r", ls="--", label="Hernquist Density")
    # ax[1].text(0.1, 0.9, f"a={fitadens}", transform=ax[1].transAxes, size=11)
    # ax[1].text(0.1, 0.85, f"aerr={np.around(np.sqrt(pcovdens), 2)}", transform=ax[1].transAxes, size=11)
    ax[1].set_ylabel(r"$\rho(r)$($M_{\odot}$/$kpc^3$)")
    ax[1].set_xlabel("r (kpc)")

    plt.show()

def main():
    # Generate combined data files
    concatFiles(os.path.abspath("./ResearchAssignments/ResearchAssignment3/MW_300.txt"), os.path.abspath("./ResearchAssignments/ResearchAssignment3/M31_300.txt"), 300)
    concatFiles(os.path.abspath("./ResearchAssignments/ResearchAssignment3/MW_450.txt"), os.path.abspath("./ResearchAssignments/ResearchAssignment3/M31_450.txt"), 450)
    concatFiles(os.path.abspath("./ResearchAssignments/ResearchAssignment3/MW_800.txt"), os.path.abspath("./ResearchAssignments/ResearchAssignment3/M31_800.txt"), 800)
    
    #collect data
    # SnapNum 300 data
    disk300CoM = CenterOfMass(os.path.abspath("./ResearchAssignments/ResearchAssignment3/MWM31_300.txt"), 2)
    CoM300 = disk300CoM.COM_P(delta=0.11) # set above 0.1 to keep every point
    halo300CoM = CenterOfMass(os.path.abspath("./ResearchAssignments/ResearchAssignment3/MWM31_300.txt"), 1)
    # SnapNum 450 data
    diskCM450 = CenterOfMass(os.path.abspath("./ResearchAssignments/ResearchAssignment3/MWM31_450.txt"), 2) # Disk for CoM calcs
    CoM450 = diskCM450.COM_P(delta=0.1)
    halo450CoM = CenterOfMass(os.path.abspath("./ResearchAssignments/ResearchAssignment3/MWM31_450.txt"), 1)
    # SnapNum 800 data
    diskCM800 = CenterOfMass(os.path.abspath("./ResearchAssignments/ResearchAssignment3/MWM31_800.txt"), 2)
    CoM800 = diskCM800.COM_P(delta=0.1)
    halo800CoM = CenterOfMass(os.path.abspath("./ResearchAssignments/ResearchAssignment3/MWM31_800.txt"), 1)
    
    # 300 Plot
    pltSnapNum(halo300CoM, CoM300, 300)
    # 450 plot
    pltSnapNum(halo450CoM, CoM450, 450)
    # 800 plot
    pltSnapNum(halo800CoM, CoM800, 800)

if __name__ == "__main__":
    main()