import numpy as np

h, T, ne, n, vt = np.loadtxt('FAL99_B', unpack = True)

#n /= 9.11e-1

#np.savetxt('FAL99_B_ntot',         np.transpose((h, T, ne, n, vt)), fmt = ('%6.2f', '%8.2f', '%10.5E', '%10.5E', '%4.2f'), delimiter = '   ')

T += 10.0

np.savetxt('FAL99_B_Tplus10', np.transpose((h, T, ne, n, vt)), fmt = ('%6.2f', '%8.2f', '%10.5E', '%10.5E', '%4.2f'), delimiter = '   ')
