#!/usr/bin/env python3.5
import subprocess
import os
import fnmatch
__author__ = 'jonathan'


def jobs_running() -> dict:
    proc = subprocess.Popen(['bjobs', '-u',  'all'], stdout=subprocess.PIPE)
    bjobs = proc.stdout.read()
    runners, penders = [], []
    for l in str(bjobs).split('\\n'):
        s = l.split()
        if len(s) > 6:
            if s[1] == 'jonatha':
                if s[2] == 'PEND':
                    penders.append(s[-4].split('.')[1])
                elif s[2] == 'RUN':
                    runners.append(s[-4].split('.')[1])
    return {'runs': runners, 'pend': penders}


def find_files(directory='src', pattern='job.*'):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename


if __name__ == '__main__':
    run_pend = jobs_running()
    dirs_running = []
    dirs_pending = []
    for file_name in find_files('/home/labs/fleishman/jonathaw/', 'job.*'):
        folder = file_name.split('/home/labs/fleishman/jonathaw/')[1].split('job.')[0]
        if folder in dirs_running or folder in dirs_pending:
            continue
        job_name = file_name.split('job.')[1]
        if job_name in run_pend['runs']:
            dirs_running.append(folder)
        if job_name in run_pend['pend']:
            dirs_pending.append(folder)
    print('Found these fodlers to be running:\n%s' % '\n'.join(dirs_running))
    print('Found these fodlers to be pending:\n%s' % '\n'.join(dirs_pending))