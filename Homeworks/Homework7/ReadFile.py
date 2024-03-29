#!/usr/bin/python3
import numpy as np
from astropy import units as u
import os
def Read(filename):
    '''
    The Read function reads in a file with assumed format as follows:
    Time    Value
    ParticleNumber  Value
    units of each quantity
    #column headers (line starts with a # and headers are ,-separated)
    data (tab separated, with every column filled)
    Parameters:
        filename: A string representing the path to the file to be read
    Returns:
        time: The time listed in the file as an astropy quantity in Myr
        numprtcl: The number of particles listed in the file
        data: All the data in the file, organized in an array of tuples where each column
              is given the name listed on the associated column header line in
              the file
    '''
    filename = os.path.abspath(filename)
    with open(filename, "r") as f: #Don't want to forget to close files
        label, value = f.readline().split() # get first line, split to get label and time
        time = float(value)*u.Myr # ensure time is proper type and units
        # get number of particles off the second line of the file
        # file format has "Total    numprtcls", so 1 is the index with the number
        numprtcl = float(f.readline().split()[1])
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3) # The data in the file, skipping the header
    return time, numprtcl, data

if __name__ == "__main__":
    time, numprtcl, data = Read("MW_000.txt")
    print(data)
        