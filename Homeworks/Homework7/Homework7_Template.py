
# # Homework 7 Template
# 
# Rixin Li & G . Besla
# 
# Make edits where instructed - look for "****", which indicates where you need to 
# add code. 

# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
from IPython.display import Latex

# **** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass2 import CenterOfMass

# **** import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass import *

# # M33AnalyticOrbit

class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, filename): # **** add inputs    
        """
        A class used to compute the analytical orbit of M33 around
        M31 and store it in a file.
        Parameters:
            filename: (string) The name of the output file
        """
        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        ### **** store the output file name
        self.filename = filename
        
        ### get the current pos/vel of M33 
        # **** create an instance of the  CenterOfMass class for M33 
        M33COM = CenterOfMass("M33_000.txt", 2)
        # **** store the position VECTOR of the M33 COM (.value to get rid of units)
        M33_Pos = M33COM.COM_P().value
        # **** store the velocity VECTOR of the M33 COM (.value to get rid of units)
        M33_v = M33COM.COM_V().value
        
        ### get the current pos/vel of M31 
        # **** create an instance of the  CenterOfMass class for M31 
        M31COM = CenterOfMass("M31_000.txt", 2)
        # **** store the position VECTOR of the M31 COM (.value to get rid of units)
        M31_Pos = M31COM.COM_P().value
        # **** store the velocity VECTOR of the M31 COM (.value to get rid of units)
        M31_v = M31COM.COM_V().value
        
        ### store the DIFFERENCE between the vectors posM33 - posM31
        # **** create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        self.r0 = M31_Pos - M33_Pos
        self.v0 = M31_v - M33_v
        ### get the mass of each component in M31 
        ### disk
        # **** self.rdisk = scale length (no units)
        self.rdisk = 5
        # **** self.Mdisk set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mdisk = ComponentMass("M31_000.txt", 2)*1e12
        ### bulge
        # **** self.rbulge = set scale length (no units)
        self.rbulge = 1
        # **** self.Mbulge  set with ComponentMass function. Remember to *1e12 to get the right units Use the right ptype
        self.Mbulge = ComponentMass("M31_000.txt", 3)*1e12
        # Halo
        # **** self.rhalo = set scale length from HW5 (no units)
        self.rhalo = 60
        # **** self.Mhalo set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mhalo = ComponentMass("M31_000.txt", 1)*1e12
     
    
    
    def HernquistAccel(self, M, r_a, r): # it is easiest if you take as an input the position VECTOR 
        """
        A function to compute the gravitational acceleration induced
        by a Hernquist profile using the usual
        a = -grad(phi)
        where phi is the Hernquist potential. This will only be able
        to compute acceleration for bulge and halo particles.
        Parameters:
            M: (float) The mass of the halo or bulge in Msun
            r_a: (float) The Hernquist scale radius in kpc
            r: (vector) The position vector to calculate the potential at.
               Each component is in kpc
        Returns:
            Hern: (vector) The acceleration vector from a Hernquist potential.
                  Each component is in kpc/Gyr^2
        """
        ### **** Store the magnitude of the position vector
        rmag = np.sqrt(np.sum([i**2 for i in r]))
        
        ### *** Store the Acceleration
        Hern =  ((-self.G*M)/(rmag*((r_a + rmag)**2))) * r #follow the formula in the HW instructions
        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        # use  -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        
        return Hern
    
    
    
    def MiyamotoNagaiAccel(self, M, r_d, r):# it is easiest if you take as an input a position VECTOR  r 
        """
        A function to compute the acceleration vector due to disk
        particles at a position r. This is done by assuming a Miyamoto-Nagai 1975 Profile
        for the disk and using the usual
        a = -grad(phi)
        with phi as the M-N potential.
        Parameters:
            M: (float) Mass of the disk in Msun
            r_d: (float) disk radius in kpc
            r: (vector) The position to compute the acceleration vector at.
               Each component is in kpc
        Returns:
            a: (vector) The acceleration vector at the given position from
               the potential. Each component is in kpc/Gyr^2
        """
        ### Acceleration **** follow the formula in the HW instructions
        # AGAIN note that we want a VECTOR to be returned  (see Hernquist instructions)
        # this can be tricky given that the z component is different than in the x or y directions. 
        # we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        # differences in the z direction:
        #  multiply the whle thing by :   np.array([1,1,ZSTUFF]) 
        # where ZSTUFF are the terms associated with the z direction
        R = np.sqrt(r[0]**2 + r[1]**2)
        z_d = self.rdisk/5
        B = r_d + np.sqrt(r[2]**2 + z_d**2)
        a = ((-self.G*M)/((R**2 + B**2)**(1.5)))*r*np.array([1, 1, ((B)/(np.sqrt(r[2]**2 + z_d**2)))])
        
        return a
        # the np.array allows for a different value for the z component of the acceleration
     
    
    def M31Accel(self, r): # input should include the position vector, r
        """
        Function to determine the total acceleration from M31
        at a given position.
        Parameters:
            r: (vector) the 3D position vector. Each component is in kpc.
        Returns:
            a: (vector) the 3D acceleration vector at this position. Each
               component is in kpc/Gyr^2
        """

        ### Call the previous functions for the halo, bulge and disk
        # **** these functions will take as inputs variable we defined in the initialization of the class like 
        # self.rdisk etc. 
        haloAccel = self.HernquistAccel(self.Mhalo, self.rhalo, r)
        bulgeAccel = self.HernquistAccel(self.Mbulge, self.rbulge, r)
        diskAccel = self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, r)
        # return the SUM of the output of the acceleration functions - this will return a VECTOR 
        return haloAccel + bulgeAccel + diskAccel
    
    
    
    def LeapFrog(self, dt, r, v): # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        """
        A function that performs one timestep of leap frog integration.
        It begins with a half timestep and determines the change in velocity
        and applies that to the position to compute a timestep.
        Parameters:
            dt: (float) The timestep of integration
            r: (vector) The initial position vector. Each component is in kpc
            v: (vector) The initial velocity vector. Each component is in kpc/Gyr
        Returns:
            rnew: (vector) The final position vector
            vnew: (vector) The final velocity vector
        """
        
        # predict the position at the next half timestep
        rhalf = r + v*(dt/2)
        
        ahalf = self.M31Accel(rhalf)
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        vnew = v + ahalf*dt
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = rhalf + vnew*(dt/2)
        
        return rnew, vnew # **** return the new position and velcoity vectors
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        """
        A function to integrate the orbit of M33 due to M31. It
        generates an output file with the orbit.
        Parameters:
            t0: (float) The starting time
            dt: (float) A time interval for integration
            tmax: (float) The time to end the integration at in Gyr
        """

        # initialize the time to the input starting time
        t = t0
        
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        orbit = np.zeros((int(tmax/dt)+2, 7))
        
        # initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        # this above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        
        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while t < tmax:  # as long as t has not exceeded the maximal time 
            
            # **** advance the time by one timestep, dt
            t += dt
            # **** store the new time in the first column of the ith row
            orbit[i, 0] = t
            
            # ***** advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like:     a,b,c = function(input)
            rnew, vnew = self.LeapFrog(dt, orbit[i-1, 1], orbit[i-1, 2])
         
    
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            orbit[i, 1:4] = rnew[0], rnew[1], rnew[2]
            orbit[i, 5:8] = vnew[0], vnew[1], vnew[2]
            # **** update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i += 1
        
        
        
        # write the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
        # there is no return function

