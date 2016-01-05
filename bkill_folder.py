import argparse
import os
from LSFManager import get_my_running_jobs
import subprocess
__author__ = 'jonathan'


def bkill_folder(folder: str) -> None:
    """
    bkills all jobs that are from folder
    :param folder: a folder address
    :return: None
    """
    global log
    pwd = os.getcwd()
    os.chdir(folder)
    folder_jobs = [a for a in os.listdir(folder) if a[:4] == 'job.']
    running_jobs = get_my_running_jobs()
    if not running_jobs:
        log += 'found NO running jobs...\n'
        return
    for folder_job in folder_jobs:
        if folder_job[4:] in running_jobs.values():
            os.system('bkill %s' % [k for k, v in running_jobs.items() if v == folder_job[4:]][0])
    os.chdir(pwd)



parser = argparse.ArgumentParser()
parser.add_argument('-folder')
args = vars(parser.parse_args())

bkill_folder(args['folder'])