import os

from tqdm import tqdm

f = open('success.log', 'r')

lines = f.readlines()

for l in tqdm(lines):

    group_ray = l.split(';')[0]

    remainder = l.split(';')[1]

    ij = remainder.split(':')[0]

    os.system('cp ./groups/' + group_ray + '/conv.out ./conv/' + ij)

f.close()
