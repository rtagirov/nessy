import numpy as np

from scipy import interpolate

atmean = 1.25650

amu = 1.66053906660e-24

boltz = 1.380649e-16

height, temp, x, ntot, x = np.loadtxt('FAL99_C_TMIN', unpack = True)

height *= 1e+5

dens = atmean * amu * ntot

pres = ntot * boltz * temp

zero = np.zeros(len(height))

coldens =  np.zeros(len(height))

heightc =  np.zeros(len(height) - 1)
coldensc = np.zeros(len(height) - 1)
dheight =  np.zeros(len(height) - 1)

for k in range(len(heightc)):

    heightc[k] = (height[k] + height[k + 1]) / 2.0
    dheight[k] = height[k + 1] - height[k]

f = interpolate.interp1d(height, dens)

densc = f(heightc)

coldensc[0] = -densc[0] * dheight[0]

for k in range(1, len(heightc)):

    coldensc[k] = coldensc[k - 1] - densc[k] * dheight[k]

f = interpolate.interp1d(heightc, coldensc, fill_value = 'extrapolate')

coldens = f(height)

if coldens[0] <= 0.0: coldens[0] = coldensc[0]

np.savetxt('FAL99_C_TMIN_kur', \
           np.column_stack([coldens, temp, pres, zero, zero, zero, zero]), \
           fmt = ('%7.5e', '%12.6f', '%7.5e', '%7.5e', '%7.5e', '%7.5e', '%7.5e'), delimiter = '  ')
