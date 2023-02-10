from Homeworks.Homework2.ReadFile import Read
from Homeworks.Homework4.CenterOfMass_Template import CenterOfMass
import astropy.units as u
import numpy as np

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
                    returned

        '''
        COM = CenterOfMass(self.filename, prtclType)
        COMPos = COM.COM_P(0.1)
        
        # get prtcls of the correct type
        typeIndex = np.where(self.data["type"] == prtclType)
        xType = self.x[typeIndex]
        yType = self.y[typeIndex]
        zType = self.z[typeIndex]
        mType = self.m[typeIndex]
        rType = np.sqrt((xType - COMPos[0])**2 + (yType - COMPos[1])**2 + (zType - COMPos[2])**2)
        # Check the type of r
        r_vals = []
        try:
            r_len = len(r)
            r_vals = r
            if len(typeIndex) == 0:
                # Return an array of zeros with the same length as
                # r if there are no particles of the given type
                return np.array([0 for i in range(len(r))])
        except:
            r_vals.append(r)
            if len(typeIndex) == 0:
                # Return 0 if there are no prtcls of the
                # given type
                return 0
        
        masses = []
        for dist in r_vals:
            insideIndex = np.where(rType < dist)
            totMass = np.sum(mType[insideIndex]*1e10)*u.Msun
            masses.append(totMass)
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

if __name__ == "__main__":
    MW = MassProfile("MW", 0)
    test_r = np.arange(0.25, 30.5, step=1.5)
    print(MW.MassEnclosed(1, test_r))