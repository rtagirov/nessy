import numpy as np

import sys

import os

from tqdm import tqdm

if len(sys.argv) < 2:

    print('File name not given. Abort.')

    sys.exit()

fname = sys.argv[1]

fname_list = fname.split('.')

dname = ''

for i in range(6):

    dname += fname_list[i] + '.'

dname = dname[:len(dname) - 1]

Nx = 512
Ny = 512

for i in range(Nx):

    os.makedirs(dname + '/' + str(i + 1), exist_ok = True)

I = np.load(fname)['I']

for i in tqdm(range(Nx)):

    for j in range(Ny):

        if I[i, j, 0] != 0.0:

            np.savetxt(dname + '/' + str(i + 1) + '/' + str(j + 1), I[i, j, :], fmt = ('%12.5E'))
