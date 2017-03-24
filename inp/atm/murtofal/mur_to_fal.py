import numpy as np
import sys

import importlib

import phys; importlib.reload(phys)

models = []

for arg in sys.argv:

    if arg != sys.argv[0]: models.append(arg)

for model in models:

    h = np.loadtxt(model, usecols = [0])
    T = np.loadtxt(model, usecols = [1])
    p = np.loadtxt(model, usecols = [2])

    n = p / (phys.k * T)

    h = abs(h - max(h))

    ne = np.ones(len(n))
    vt = np.ones(len(n))

    np.savetxt(model + '_fal', np.transpose((h, T, ne, n, vt)), fmt = ('%9.2f', '%9.2f', '%15.7e', '%15.7e', '%15.7e'), delimiter = '  ')
