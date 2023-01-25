#!/usr/bin/python3
import numpy as np
import astropy.units as u
from ReadFile import Read

def ParticleInfo(filename, prtclType, prtclNum):
    '''
    This function gets the data of one particle in the
    given snapshot file. The particle is specified by its type
    (1=dark matter, 2=disk star, 3=bulge star) and number.
    The relevant data in this case is the velocity and position magnitudes,
    calculated from the x,y,z,vx,vy,vz values in the snapshot file.
    Parameters:
        filename: A string representing the path to the snapshot file
        prtclType: An integer representing the type of the particle
        prtclNum: An integer specifying which particle to look at
                  (prtclNum = 10 means"Get info on the 10th particle of type prtclType")
    Returns:
        A tuple of Astropy Quantites that represent the following:
        Position of the particle in kpc
        Position of the particle in lyr
        Velocity of the particle in km/s
        Mass of the particle in Msun
    '''
    # read the file into an array with column headers
    time, numprtcls, data = Read(filename) 
    # get the index where the 'prctlNum'th particle of the correct type is located
    requestedPrtclIndex = np.where(data["type"] == prtclType)[0][prtclNum - 1]    
    # pull all of the relevant quantities from the data
    prtclM = data["m"][requestedPrtclIndex] #mass
    # position values
    prtclx = data["x"][requestedPrtclIndex]
    prtcly = data["y"][requestedPrtclIndex]
    prtclz = data["z"][requestedPrtclIndex]
    # velocity values
    prtclvx = data["vx"][requestedPrtclIndex]
    prtclvy = data["vy"][requestedPrtclIndex]
    prtclvz = data["vz"][requestedPrtclIndex]
    # compute magnitudes from components, setting the proper units
    poskpc = np.sqrt(prtclx**2 + prtcly**2 + prtclz**2)*u.kpc
    veloMag = np.sqrt(prtclvx**2 + prtclvy**2 + prtclvz**2)*u.km / (1*u.s)
    mass = prtclM*1e10*u.Msun
    # convert kpc to lyr
    posly = poskpc.to(u.lyr)
    return np.around(poskpc, 3), np.around(posly, 3), np.around(veloMag, 3), np.around(mass, 3)


if __name__ == "__main__":
    prtclType = 2
    prtclNum = 100
    info = ParticleInfo("C:\\Users\\phart\\Desktop\\College Stuff\\Senior Year\\400B\\400B_2023_hartman\\Homeworks\\Homework2\MW_000.txt", prtclType, prtclNum)
    print(f"Info for Particle {prtclNum} of type {prtclType}\nParticle Distance: {info[0].value} {info[0].unit}, {info[1].value} {info[1].unit}\nParticle Velocity: {info[2].value} {info[2].unit}\nParticle Mass: {info[3].value} {info[3].unit}")