import os

from tqdm import tqdm

f = open('fail.log', 'r')

lines = f.readlines()

for l in tqdm(lines):

    if 'convergence' in l:

       group_ray = l.split(';')[0]

       remainder = l.split(';')[1]

       ij = remainder.split(':')[0]

#       i = ij.split('.')[0]
#       j = ij.split('.')[1]

       os.system('cp ./groups/' + group_ray + '/conv.out ./noconv/' + ij)

f.close()
