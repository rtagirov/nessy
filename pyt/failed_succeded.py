import os

from tqdm import tqdm

def job_id(l):

    group_ray = l.split(';')[0]

    fioss_log = open('./groups/' + group_ray + '/fioss.log', 'r')

    fioss_log_lines = fioss_log.readlines()

    fioss_log.close()

    job_line = fioss_log_lines[0]

    slash_split = job_line.split('/')

    dot_split = slash_split[2].split('.')

    return dot_split[0]

fail_file = open('fail.log', 'r')
success_file = open('success.log', 'r')

fail_lines = fail_file.readlines()
success_lines = success_file.readlines()

fail_file.close()
success_file.close()

fail_id = open('fail.id', 'w')
success_id = open('success.id', 'w')

for l in tqdm(fail_lines):

    if 'unknown' in l:

        fail_id.write(job_id(l) + '\n')

for l in tqdm(success_lines):

    success_id.write(job_id(l) + '\n')

fail_id.close()
success_id.close()
