import glob
import os
import sys

import numpy as np
import netCDF4 as nc

from tqdm import tqdm

from multiprocessing import Pool

if len(sys.argv) < 2:

    print('nproc not provided. abort.')

    sys.exit()

nproc = int(sys.argv[1])

Nx = 512
Ny = 512
Nz = 120

def read_spec(ray_number):

    i = int(ray_number / 512) + 1

    j = ray_number % 512

    if j == 0:

        i -= 1

        j = 512

    name = './spec/' + str(i) + '.' + str(j)

    intensity = np.zeros(Nz)

    if os.path.isfile(name): intensity = np.genfromtxt(name)

    return i, j, intensity

I = np.zeros((Nx, Ny, Nz))

f = open('specification.dat', 'r')

specifications = f.readlines()

rev_cube_frac = int(specifications[3])

ray_numbers = np.arange(rev_cube_frac, Nx * Ny + rev_cube_frac, rev_cube_frac)

with Pool(processes = nproc) as p:

    maximum = len(ray_numbers)

    n_chunks = 1

    with tqdm(total = maximum, position = 0) as pbar:

        results = p.imap(read_spec, ray_numbers, chunksize = 1)

        for i, result in enumerate(results):

            i, j, intensity = result

            I[i - 1, j - 1, :] = intensity

            pbar.update()

    p.close()
    p.join()

spectral_type = specifications[0].strip('\n')
magnetisation = specifications[1].strip('\n')
snapshot =      specifications[2].strip('\n')
regime =        specifications[4].strip('\n')
angle =         specifications[5].strip('\n')

np.savez('./spec/' + spectral_type + '.' 
                   + magnetisation + '.' 
                   + snapshot + '.' 
                   + str(rev_cube_frac) + '.'
                   + regime + '.' 
                   + angle, I = I)
