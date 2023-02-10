#!/usr/bin/python3
import numpy as np
from Homeworks.Homework2.ReadFile import Read
import astropy.units as u

def ComponentMass(filename, prtclType):
    '''
    This function reads in the total mass of a given type of particle
    from a given file.

    Parmeters:
        filename: (String) The name of the file to check
        prtclType: (Integer) The type of particle to get the mass of
    Returns:
        (float) The total mass of the particles of the given type in
        units of 10^{12} Msun
    '''
    time, numprtcl, data = Read(filename)
    reqPrtclTypeIndex = np.where(data["type"] == prtclType)[0]
    # ^ index 0 is the array, 1 is the type
    totMass = 0*u.Msun
    for index in reqPrtclTypeIndex:
        totMass += data['m'][index]*u.Msun
    totMass = totMass / 100 # gets from 1e10 to 1e12
    return np.around(totMass, 3)

def printTable(filenames):
    '''
    This function generates the rows and columns of a LaTeX table from
    given filenames. The table is considering the masses of different types of
    particles, with columns as follows:
    Galaxy Name, Halo Mass (10^12 Msun), Disk Mass (10^12 Msun), Bulge Mass (10^12 Msun), Total Mass (10^12 Msun), fbar
    with fbar as the baryon fraction.
    Parameters:
        filenames: (dict) A dictionary mapping the Galaxy Name as a string to the filepath string associated with it
    Returns:
        table: (string) The table rows in LaTeX syntax, with values filled in.
    '''
    #add table column names
    table = "Galaxy Name & Halo Mass ($10^{12}$ $M_{\odot}$) & Disk Mass ($10^{12}$ $M_{\odot}$) & Bulge Mass ($10^{12}$ $M_{\odot}$) & Total ($10^{12}$ $M_{\odot}$) & $f_{bar}$\\\\\n"
    groupDiskMass = 0
    groupHaloMass = 0
    groupBulgeMass = 0
    # Get total nasses for each row
    for galaxy in filenames:
        HaloMass = ComponentMass(filenames[galaxy], 1).value
        groupHaloMass += HaloMass
        DiskMass = ComponentMass(filenames[galaxy], 2).value
        groupDiskMass += DiskMass
        BulgeMass = ComponentMass(filenames[galaxy], 3).value
        groupBulgeMass += BulgeMass
        totalMass = BulgeMass + DiskMass + HaloMass
        fbar = np.around(((BulgeMass + DiskMass) / totalMass), 3)
        table += f"{galaxy} & {HaloMass} & {DiskMass} & {BulgeMass} & {totalMass} & {fbar}\\\\\n"
    # Round total masses
    groupBulgeMass = np.around(groupBulgeMass, 3)
    groupDiskMass = np.around(groupDiskMass, 3)
    groupHaloMass = np.around(groupHaloMass, 3)
    # compute totalMass and fbar for local group
    groupTotalMass = groupHaloMass + groupDiskMass + groupBulgeMass
    groupfbar = (groupDiskMass + groupBulgeMass) / groupTotalMass
    # Add Local Group Row
    table += f"Local Group & {groupHaloMass} & {groupDiskMass} & {groupBulgeMass} & {np.around(groupTotalMass, 3)} & {np.around(groupfbar, 3)}\\\\"
    return table

if __name__ == "__main__":
    # filenames associated with different galaxies
    MW = "C:\\Users\\phart\\Desktop\\College Stuff\\Senior Year\\400B\\400B_2023_hartman\\MW_000.txt"
    M31 = 'C:\\Users\\phart\\Desktop\\College Stuff\\Senior Year\\400B\\400B_2023_hartman\\M31_000.txt'
    M33 = 'C:\\Users\\phart\\Desktop\\College Stuff\\Senior Year\\400B\\400B_2023_hartman\\M33_000.txt'
    # Generating a LaTeX table with the galaxy names given
    print(printTable({"Milky Way": MW, "M31":M31, "M33":M33}))