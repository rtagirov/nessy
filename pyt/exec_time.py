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

if len(sys.argv) < 4:

    print('step not provided. abort.')

    sys.exit()

step = int(sys.argv[3])

def num_lines(f):

    return sum(1 for line in open(f))

def get_time(folder):

    dpn = 0
    nli = 0

    fiet = 0.0
    hmet = 0.0

    fime = 0.0
    hmme = 0.0

    if os.path.isfile(folder + '/hminus.log') and \
       os.path.isfile(folder + '/fioss.log'):

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

        hm_mem_line = hm_lines[13]
        fi_mem_line = fi_lines[13]

        hmet = float(hm_user_time_line.split(':')[1].strip('\n'))
        fiet = float(fi_user_time_line.split(':')[1].strip('\n'))

        hmme = float(hm_mem_line.split(':')[1].strip('\n'))
        fime = float(fi_mem_line.split(':')[1].strip('\n'))

        atm = open(folder + '/atm.inp', 'r')

        atm_lines = atm.readlines()

        atm.close()

        lin = []

        for l in atm_lines:

            if len(l.split(' ')) == 2:

                lin.append(l.strip('\n'))

        ray_number = int(np.loadtxt(folder + '/rn.inp'))

        dpn = int(lin[ray_number - 1].split(' ')[1])

        nli = num_lines(folder + '/conv.out')

    return dpn, nli, hmet, fiet, hmme, fime

def exec_time(l):

    group_ray = ''
    xy = ''

    dpn = 0
    nli = 0

    hmet = 0.0
    fiet = 0.0

    hmme = 0.0
    fime = 0.0

    if regime == 'mr':

        if l != '\n' and len(l.split(';')) == 2:

            group_ray = l.split(';')[0]

            if len(group_ray.split('/')) == 2:

                group = group_ray.split('/')[0]
                ray =   group_ray.split('/')[1]

                if group and ray:

                    folder = './groups/' + group + '/' + ray

                    if os.path.isdir(folder):

                        dpn, nli, hmet, fiet, hmme, fime = get_time(folder)

    if regime == 'hk':

        if l != '\n' and len(l.split(':')) == 2:

            xy = l.split(':')[0]

            if len(xy.split('.')) == 2:

                folder = './groups_hk/' + xy

                if os.path.isdir(folder):

                    dpn, nli, hmet, fiet, hmme, fime = get_time(folder)

    if regime == 'mr': return group_ray, dpn, nli, hmet, fiet, hmme, fime
    if regime == 'hk': return xy,        dpn, nli, hmet, fiet, hmme, fime

if regime == 'mr': f = open('success.log',    'r')
if regime == 'hk': f = open('success.hk.log', 'r')

lines = f.readlines()[::step]

f.close()

rid0 = []
rid1 = []
num_depth_points = []
num_lambda_iter = []
hm_exec_time = []
fi_exec_time = []
hm_mem = []
fi_mem = []

with Pool(processes = nproc) as p:

    maximum = len(lines)

    n_chunks = 1

    with tqdm(total = maximum, position = 0) as pbar:

        results = p.imap(exec_time, lines, chunksize = 1)

        for i, result in enumerate(results):

            rid, dpn, nli, hmet, fiet, hmme, fime = result

            if dpn > 0 and hmet > 0.0 and fiet > 0.0:

                num_depth_points.append(dpn)
                num_lambda_iter.append(nli)
                hm_exec_time.append(hmet)
                fi_exec_time.append(fiet)
                hm_mem.append(hmme)
                fi_mem.append(fime)

                if regime == 'mr':

                    rid0.append(int(rid.split('/')[0]))
                    rid1.append(int(rid.split('/')[1]))

                if regime == 'hk':

                    rid0.append(int(rid.split('.')[0]))
                    rid1.append(int(rid.split('.')[1]))

            pbar.update()

    p.close()
    p.join()

rid0 = np.array(rid0)
rid1 = np.array(rid1)

num_depth_points = np.array(num_depth_points)
num_lambda_iter = np.array(num_lambda_iter)

hm_exec_time = np.array(hm_exec_time) / 60
fi_exec_time = np.array(fi_exec_time) / 60

hm_mem = np.array(hm_mem) / 1000
fi_mem = np.array(fi_mem) / 1000

np.savetxt('exec_time.out', np.transpose((rid0, rid1, num_depth_points, num_lambda_iter, hm_exec_time, fi_exec_time, hm_mem, fi_mem)), fmt = ('%4i', '%4i', '%4i', '%4i', '%6.3f', '%6.3f', '%6.3f', '%6.3f'), delimiter = '  ')
