import numpy as np

from tqdm import tqdm

from multiprocessing import Pool

import os

import sys

if len(sys.argv) < 2:

    print('nproc not provided. abort.')

    sys.exit()

nproc = int(sys.argv[1])

if len(sys.argv) < 3:

    print('regime not provided. abort.')

    sys.exit()

regime = sys.argv[2]

def exec_time(l):

    dpn = 0
    hmet = 0.0
    fiet = 0.0

    if regime == 'mr':

        if l != '\n' and len(l.split(';')) == 2:

            group_ray = l.split(';')[0]

            if len(group_ray.split(',')) == 2:

                group = group_ray.split(',')[0]
                ray =   group_ray.split(',')[1]

                if group and ray:

                    folder = './groups/' + group + '/' + ray

    if regime == 'hk':

        if l != '\n' and len(l.split(':')) == 2:

            xy = l.split(':')[0]

            if len(xy.split('.')) == 2:

                folder = './groups_hk/' + xy

    if os.path.isdir(folder):

        hmlog = open(folder + '/hminus.log', 'r')
        filog = open(folder + '/fioss.log',  'r')

        hm_lines = hmlog.readlines()
        fi_lines = filog.readlines()

        hmlog.close()
        filog.close()

        hm_lines.reverse()
        fi_lines.reverse()

        hm_user_time_line = hm_lines[21]
        fi_user_time_line = fi_lines[21]

        hmet = float(hm_user_time_line.split(':')[1].strip('\n'))
        fiet = float(fi_user_time_line.split(':')[1].strip('\n'))

        atm = open(folder + '/atm.inp', 'r')

        atm_lines = atm.readlines()

        atm.close()

        lin = []

        for l in atm_lines:

            if len(l.split(' ')) == 2:

                lin.append(l.strip('\n'))

#        print(lin)

        ray_number = int(np.loadtxt(folder + '/rn.inp'))

#        print(ray_number)

        dpn = int(lin[ray_number - 1].split(' ')[1])

#        print(dpn)

    return dpn, hmet, fiet

if regime == 'mr': f = open('success.log',    'r')
if regime == 'hk': f = open('success.hk.log', 'r')

lines = f.readlines()

f.close()

num_depth_points = []
hm_exec_time = []
fi_exec_time = []

with Pool(processes = nproc) as p:

    maximum = len(lines)

    n_chunks = 1

    with tqdm(total = maximum, position = 0) as pbar:

        results = p.imap(exec_time, lines, chunksize = 1)

        for i, result in enumerate(results):

            dpn, hmet, fiet = result

            if dpn > 0:    num_depth_points.append(dpn)
            if hmet > 0.0: hm_exec_time.append(hmet)
            if fiet > 0.0: fi_exec_time.append(fiet)

            pbar.update()

    p.close()
    p.join()

num_depth_points = np.array(num_depth_points)

hm_exec_time = np.array(hm_exec_time) / 60
fi_exec_time = np.array(fi_exec_time) / 60

np.savetxt('exec_time.out', np.transpose((num_depth_points, hm_exec_time, fi_exec_time)), fmt = ('%4i', '%6.3f', '%6.3f'), delimiter = '  ')
