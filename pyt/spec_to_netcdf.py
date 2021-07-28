import glob
import os

import numpy as np
import netCDF4 as nc

from tqdm import tqdm

def read_spec(coord):

    I = np.genfromtxt('./spec/' + coord)

    return I

Nx = 512 / cube_fraction
Ny = 512 / cube_fraction
Nz = 120

I = np.zeros((Nx, Ny, Nz))

for ray_number in 

for i in tqdm(range(1, Nx + 1)):

    for j in range(1, Ny + 1):

        coord = str(i) + '.' + str(j)

        if os.path.isfile(name): I[i - 1, j - 1, :] = 
