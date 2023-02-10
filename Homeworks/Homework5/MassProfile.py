import astropy.units as u
import numpy as np
from ReadFile import Read
from CenterOfMass_Template import CenterOfMass
from astropy.constants import G
import matplotlib.pyplot as plt

class MassProfile:
    def __init__(self, galaxy, snap):
        strSnapNum = '000' + str(snap)
        # remove all bu the last 3 digits
        strSnapNum = strSnapNum[-3:]
        self.filename = f"{galaxy}_{strSnapNum}.txt"
        self.gname = galaxy
        self.time, self.total, self.data = Read(self.filename)
        self.x = self.data["x"]*u.kpc
        self.y = self.data["y"]*u.kpc
        self.z = self.data["z"]*u.kpc
        self.m = self.data["m"]
        self.vx = self.data["vx"]*u.km/u.s
        self.vy = self.data["vy"]*u.km/u.s
        self.vz = self.data["vz"]*u.km/u.s


    def MassEnclosed(self, prtclType, r):
        '''
        This function returns the mass enclosed within a disk
        of given radius the centered at the COM of the galaxy.
        This function will return 0 if there are no particles
        of the given type in the galaxy.
        Paramaters:
            prtclType: (integer) The type of particle to check
                       the mass of. Must be either 1, 2, or 3, with
                       1 corresponding to the Dark Matter Halo,
                       2 corresponding to disk stars, and
                       3 corresponding to bulge stars
            r: (float or array of floats) the radius or radii
               to check.
        Returns:
            masses: (Astropy Quantity or Array of Astropy Quantities)
                    The total mass of particles of the given type within
                    the corresponding radii in solar masses. If there
                    is only one radius given, a single quantity is
                    returned. Returns 0 or an array of 0s if there are
                    no particles of the given type in the snapshot.
        '''
        # Check the type of r and make it an array
        # Check if the prtclType has prtcls in the galaxy
        r_vals = []
        prtclIndex = np.where(self.data["type"] == prtclType)
        try:
            r_len = len(r) # Check for array
            r_vals = r
            if len(prtclIndex[0]) == 0:
                # Return an array of 0s
                return np.zeros(len(r))
        except:
            r_vals.append(r) # Make the float an array
            r_vals = np.array(r_vals)
            if len(prtclIndex)[0] == 0:
                # Return a float like what was entered
                return 0.
        COM = CenterOfMass(self.filename, prtclType)
        COMPos = COM.COM_P(0.1)
        
        # get prtcls of the correct type
        typeIndex = np.where(self.data["type"] == prtclType)[0]
        xType = self.x[typeIndex]
        yType = self.y[typeIndex]
        zType = self.z[typeIndex]
        mType = self.m[typeIndex]
        rType = np.sqrt((xType - COMPos[0])**2 + (yType - COMPos[1])**2 + (zType - COMPos[2])**2)
        # Get masses
        masses = []
        for dist in r_vals:
            insideIndex = np.where(rType.value < dist)
            totMass = np.sum(mType[insideIndex]*1e10)
            masses.append(totMass)
        masses = np.array(masses)*u.Msun
        if len(masses) == 1:
            return masses[0]
        else:
            return masses
    
    def MassEnclosedTotal(self, r):
        '''
        Determines the total mass of all particles enclosed within the
        given radii.

        Parameters:
            r: (Array of Floats) The radii to use for mass determination.
               The units are assumed to be in kpc
        
        Returns:
            totMass: (Array of Astropy Quantities) The total mass enclosed
                     within each given radius in solar masses.
        '''
        Enclosed1 = self.MassEnclosed(1, r)
        Enclosed2 = self.MassEnclosed(2, r)
        Enclosed3 = self.MassEnclosed(3, r)
        totMass = Enclosed1 + Enclosed2 + Enclosed3
        return totMass

    def HernquistMass(self, r, a, Mhalo): # TODO
        '''
        Computes the mass enclosed within a given radius using a theoretical
        Hernquist Profile:
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
            EnclosedM: (Astropy Quantity or Array of Astropy Quantities)
                       the total enclosed mass at the given radii
                       in Solar Masses.
        '''
        return
    
    def HernquistVCirc(self, r, a, Mhalo):# TODO
        '''
        Computes the Circular Speed at a given radius using a Hernquist
        profile for the galaxy.
        Parameters:
            r: (float or Array of floats) The radii to compute the
               circular speed at. Assumed to be in kpc
            a: (float) The Hernquist scale height
            Mhalo: (float) The dark matter halo mass of the galaxy
        Returns:
            speeds: (Astropy Quantity or Array of Astropy Quantities)
                    The circular speeds at the given radii in km/s.
                    Returned as a single Astropy Quantity if the input
                    r is a float, or an array of Astropy Quantities if the
                    input was an array of floats.
        '''
        return

    def CircularVelocity(self, prtclType, r):
        '''
        Determines the circular orbit velocities for particles of
        the given particle type at the given radii.
        Parameters:
            prtclType: (integer) The type of particle to check
                the mass of. Must be either 1, 2, or 3, with
                1 corresponding to the Dark Matter Halo,
                2 corresponding to disk stars, and
                3 corresponding to bulge stars
            r: (Array of floats) radii to compute the cicular velocity at.
               Assumed to be in kpc
        Returns:
            speeds: (Array of Astropy Quantities) The circular speeds at
                    each radius in km/s
        '''
        prtclTypeMass = self.MassEnclosed(prtclType, r)
        Grav = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        r = r*u.kpc
        # compute and return circular speeds
        # newtonian gravity is assumed so mv^2/r = GMm/r^2 -> v^2 = GM/r
        return np.sqrt(Grav*prtclTypeMass/r)
    
    def CircularVelocityTotal(self, r):
        '''
        Returns the circular velocity at a given radii
        including the mass constribution from each component of the
        galaxy (the halo, disk, and bulge mass)
        Parameters:
            r: (Array of floats) an array of radii to get velocities at.
               Assumed to be given in kpc
        Returns:
            speeds: (Array of Astropy Quantities) Circular orbital speeds
                    for each given radius
        '''
        M = self.MassEnclosedTotal(r)
        Grav = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        r = r*u.kpc
        # compute and return circular speeds
        # newtonian gravity is assumed so mv^2/r = GMm/r^2 -> v^2 = GM/r
        return np.sqrt(Grav*M/r)


if __name__ == "__main__":
    test_r = np.arange(0.1, 25, step=0.05)
    # MW
    MW = MassProfile("MW", 0)
    MWhaloProfile = MW.MassEnclosed(1, test_r)
    MWdiskProfile = MW.MassEnclosed(2, test_r)
    MWBulgeProfile = MW.MassEnclosed(3, test_r)
    MWtotProfile = MW.MassEnclosedTotal(test_r)
    # M31
    M31 = MassProfile("M31", 0)
    M31haloProfile = M31.MassEnclosed(1, test_r)
    M31diskProfile = M31.MassEnclosed(2, test_r)
    M31BulgeProfile = M31.MassEnclosed(3, test_r)
    M31totProfile = M31.MassEnclosedTotal(test_r)
    # M33
    M33 = MassProfile("M33", 0)
    M33haloProfile = M33.MassEnclosed(1, test_r)
    M33diskProfile = M33.MassEnclosed(2, test_r)
    M33BulgeProfile = M33.MassEnclosed(3, test_r)
    M33totProfile = M33.MassEnclosedTotal(test_r)

    # # Plotting
    fig, ax = plt.subplots(1, 3, sharey=True, figsize = (14, 6))
    # # MW
    ax[0].semilogy(test_r, MWtotProfile, color="k", label="Total")
    ax[0].semilogy(test_r, MWhaloProfile, linestyle="dotted", color="magenta", label="Halo")
    ax[0].semilogy(test_r, MWdiskProfile, linestyle="--", color="blue", label="Disk")
    ax[0].semilogy(test_r, MWBulgeProfile, linestyle="dashdot", color="cyan", label="Bulge")
    # M31
    ax[1].semilogy(test_r, M31totProfile, color="k", label="Total")
    ax[1].semilogy(test_r, M31haloProfile, linestyle="dotted", color="magenta", label="Halo")
    ax[1].semilogy(test_r, M31diskProfile, linestyle="--", color="blue", label="Disk")
    ax[1].semilogy(test_r, M31BulgeProfile, linestyle="dashdot", color="cyan", label="Bulge")
    # M33
    ax[2].semilogy(test_r, M33totProfile, color="k", label="Total")
    ax[2].semilogy(test_r, M33haloProfile, linestyle="dotted", color="magenta", label="Halo")
    ax[2].semilogy(test_r, M33diskProfile, linestyle="--", color="blue", label="Disk")
    ax[2].semilogy(test_r, M33BulgeProfile, linestyle="dashdot", color="cyan", label="Bulge")
    # Plot labels
    ax[0].set_xlabel("r (kpc)")
    ax[1].set_xlabel("r (kpc)")
    ax[2].set_xlabel("r (kpc)")
    ax[0].set_ylabel(r"$M_{enclosed}$ ($M_{\odot}$)")
    ax[0].text(0.1, 0.9, "Milky Way", transform=ax[0].transAxes, size=14)
    ax[1].text(0.1, 0.9, "M31", transform=ax[1].transAxes, size=14)
    ax[2].text(0.1, 0.9, "M33", transform=ax[2].transAxes, size=14)
    fig.suptitle("Mass Profiles", size=18)
    ax[2].legend(loc="lower right")
    plt.show()
    

