# Homework 6 Template
# G. Besla & R. Li

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
import os

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as mpatches

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass2 import CenterOfMass




def OrbitCOM(galaxy, start, end, n):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
          galaxy: (string) the name of the galaxy, e.g. "MW"
          start: (integer) The number of the first snapshot to be read in
          end: (integer) The number of the last snapshot to be read in
          n: (integer) The intervals over which the COM will be returned
    """
    "C:\\Users\\phart\\Desktop\\College Stuff\\Senior Year\\400B\\400B_2023_hartman\\Homeworks\\Homework6\\data\\{galaxy}_{snapName}.txt"
    # compose the filename for output
    fileout = f"C:\\Users\\phart\\Desktop\\College Stuff\\Senior Year\\400B\\400B_2023_hartman\\Homeworks\\Homework6\\output\\Orbit_{galaxy}.txt"
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    delta = 0.1
    if galaxy != "M33":
        volDec = 2
    else:
        volDec = 4
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    snap_ids = np.arange(start, end, step=n)
    if len(snap_ids) == 0:
        return "Invalid start and end snapshots"
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros([len(snap_ids), 7])
    
    # a for loop 
    for i, snap_id in enumerate(snap_ids): # loop over files
        
        # compose the data filename (be careful about the folder)
        snapName = '000' + str(snap_id)
        snapName = snapName[-3:]
        # Initialize an instance of CenterOfMass class, using disk particles
        iCOM = CenterOfMass(f"C:\\Users\\phart\\Desktop\\College Stuff\\Senior Year\\400B\\400B_2023_hartman\\Homeworks\\Homework6\\data\\{galaxy}_{snapName}.txt", 2)
        # Store the COM pos and vel. Remember that now COM_P required VolDec
        COM_P = iCOM.COM_P(delta, volDec)
        COM_v = iCOM.COM_V(COM_P[0], COM_P[1], COM_P[2])
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        time = int(iCOM.time.value) / 1000
        orbit[i] = time, *(COM_P.value), *(COM_v.value)
        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))




# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 
start = 0
end = 800
step = 5
# OrbitCOM("MW", start, end, step)
# OrbitCOM("M31", start, end, step)
# OrbitCOM("M33", start, end, step)


# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt
MWdata = np.genfromtxt(os.path.abspath("./output/Orbit_MW.txt"), dtype=None, names=True)
M31Data = np.genfromtxt(os.path.abspath("./output/Orbit_M31.txt"), dtype=None, names=True)
M33Data = np.genfromtxt(os.path.abspath("./output/Orbit_M33.txt"), dtype=None, names=True)

# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  

diff2vecs = lambda v1, v2: np.sqrt(np.sum([i**2 for i in (v1-v2)]))

# def diff2vecs(v1, v2):
#     diff = []
#     for i in range(len(v1)):
#         vecDiff = v1[i] - v2[i]
#         diff.append(np.sqrt(np.sum([i**2 for i in vecDiff])))

# Determine the magnitude of the relative position and velocities 

# of MW and M31
posVecMW = np.array([MWdata['x'], MWdata['y'], MWdata['z']])
posVecM31 = np.array([M31Data['x'], M31Data['y'], M31Data['z']])
differencePMWM31 = []
for i in range(len(MWdata['t'])):
    diffAtTimet = diff2vecs(np.array(posVecMW[:, i]), np.array(posVecM31[:, i]))
    differencePMWM31.append(diffAtTimet)
differencePMWM31 = np.array(differencePMWM31)

veloVecMW = np.array([MWdata['vx'], MWdata['vy'], MWdata['vz']])
veloVecM31 = np.array([M31Data['vx'], M31Data['vy'], M31Data['vz']])
differenceVMWM31 = []
for i in range(len(MWdata['t'])):
    diffAtTimet = diff2vecs(np.array(veloVecMW[:, i]), np.array(veloVecM31[:, i]))
    differenceVMWM31.append(diffAtTimet)
differenceVMWM31 = np.array(differenceVMWM31)

# of M33 and M31
posVecM33 = np.array([M33Data['x'], M33Data['y'], M33Data['z']])
differencePM31M33 = []
for i in range(len(MWdata['t'])):
    diffAtTimet = diff2vecs(np.array(posVecM33[:, i]), np.array(posVecM31[:, i]))
    differencePM31M33.append(diffAtTimet)
differencePM31M33 = np.array(differencePM31M33)

veloVecM33 = np.array([M33Data['vx'], M33Data['vy'], M33Data['vz']])
differenceVM31M33 = []
for i in range(len(MWdata['t'])):
    diffAtTimet =  diff2vecs(np.array(veloVecM33[:, i]), np.array(veloVecM31[:, i]))
    differenceVM31M33.append(diffAtTimet)
differenceVM31M33 = np.array(differenceVM31M33)


fig, ax = plt.subplots(2, 1, sharex=True)

# Plot the Orbit of the galaxies 
#################################
ax[0].plot(MWdata['t'], differencePMWM31, color="magenta")
ax[0].plot(MWdata['t'], differencePM31M33, color = "cyan")

# Plot the orbital velocities of the galaxies 
#################################
ax[1].plot(MWdata['t'], differenceVMWM31, color="magenta")
ax[1].plot(MWdata['t'], differenceVM31M33, color="cyan")


ax[1].set_xlabel("time (Gyr)")
ax[0].set_ylabel("Position Separation (kpc)")
ax[1].set_ylabel("Vedlocity Separation (km/s)")

magenta_patch = mpatches.Patch(color = "magenta", label="MW and M31")
cyan_patch = mpatches.Patch(color="cyan", label="M31 and M33")

ax[0].legend(handles=[magenta_patch, cyan_patch], loc="upper right")

fig.suptitle("Separation Magnitude")
plt.savefig("PositionVelocityFigures.png")