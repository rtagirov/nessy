import numpy as np

from tqdm import tqdm

from multiprocessing import Pool

import os

def exec_time(l):

    hmet = 0.0
    fiet = 0.0

    if l != '\n' and len(l.split(';')) == 2:

        group_ray = l.split(';')[0]

        if len(group_ray.split(',')) == 2:

            group = group_ray.split(',')[0]
            ray =   group_ray.split(',')[1]

            if group and ray:

                folder = './groups/' + group + '/' + ray

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

    return hmet, fiet

f = open('success.log', 'r')
#f = open('success_test.log', 'r')

lines = f.readlines()

f.close()

hm_exec_time = []
fi_exec_time = []

with Pool(processes = 16) as p:

    maximum = len(lines)

    n_chunks = 1

    with tqdm(total = maximum, position = 0) as pbar:

        results = p.imap(exec_time, lines, chunksize = 1)

        for i, result in enumerate(results):

            hmet, fiet = result

            if hmet > 0.0: hm_exec_time.append(hmet)
            if fiet > 0.0: fi_exec_time.append(fiet)

            pbar.update()

    p.close()
    p.join()

hm_exec_time = np.array(hm_exec_time) / 60
fi_exec_time = np.array(fi_exec_time) / 60

np.savetxt('exec_time.out', np.transpose((hm_exec_time, fi_exec_time)), fmt = ('%6.3f', '%6.3f'), delimiter = '  ')
