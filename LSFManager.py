#!/usr/bin/env python3.5
"""
A script to manage LSF runs.
"""
import os
import sys
import time
import smtplib
import subprocess
from result_processor import process_folder
from DoCohResultProcessor import generate_run_filters, all_who_pass_run_filters
from RosettaFilter import score2dict

__author__ = 'jonathan'
global log


def main():
    global log
    original_pwd = os.getcwd()
    running = get_running_folders()
    pending = get_pending_folders()
    run_filters = generate_run_filters(args={'ddg': 25.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6, 'buried_2': 3})
    for runner in running:
        os.chdir(runner)
        if is_folder_finished(runner):
            print('processing %s' % runner)
            log += 'folder is finished, processing %s\n' % runner.split('/')[-2]
            process_folder({'folder': runner, 'force_process': False, 'remove_pdbs': False})
            move_to_processed(runner, running, 0)
        else:
            score_dict = folder_scores(runner)
            passed, failed = all_who_pass_run_filters({}, score_dict, run_filters)
            log += 'passed %i, failed %i' % (len(passed), len(failed))
            if len(passed) >= 50:
                print('found enough passed, stopping folder %s' % runner.split('/')[-1])
                log += 'found enough passed, stopping folder\n'
                bkill_folder(runner)
                process_folder({'folder': runner, 'force_process': True, 'remove_pdbs': False})
                move_to_processed(runner, running, len(passed))
            else:
                print('not enough finished, letting him be %s, found %i passed and %i failed' %
                      (runner.split('/')[-2], len(passed), len(failed)))
                log += 'not enough finished, letting him be %s, found %i passed and %i failed\n' % \
                       (runner.split('/')[-2], len(passed), len(failed))
        os.chdir(original_pwd)

    # for pender in pending:
    #     os.chdir(pender)
    #     lsf_status, pends = how_many_queue()
    #     if lsf_status['fleishman'] < 12000:
    #         print('found %i jobs in fleishman, submitting %s' % (lsf_status['fleishman'], pender))
    #         log += 'found %i jobs in fleishman, submitting %s\n' % (lsf_status['fleishman'], pender)
    #         submit_folder(pender)
    #         move_pender_to_runner(pending, pender)
    #     os.chdir(original_pwd)

    # lsf_status, pends = how_many_queue()
    # if lsf_status['new-all.q'] <= 2000:
    #     bswitch_pends(pends, 2000-lsf_status['new-all.q'])
    os.chdir(original_pwd)


def bswitch_pends(pends: list, num: int) -> None:
    """
    :param pends: a list of pending jobs
    :param num: number of jobs to switch
    :return: None
    """
    global log
    log += 'switching %i jobs to new-all.q\n' % num
    for pnd in pends[:num]:
        os.system('bmod -q new-all.q %s' % pnd)


def how_many_queue() -> (dict, list):
    """
    checks how many jobs are running in each queue
    :return: {fleishman: #jobs, new-all.q: #jobs} and a list of job id that are pending in fleishman
    """
    proc = subprocess.Popen(['bjobs', '-u',  'all'], stdout=subprocess.PIPE)
    bjobs = proc.stdout.read()
    results = {'fleishman': 0, 'new-all.q': 0}
    pends = []
    for l in str(bjobs).split('\\n'):
        s = l.split()
        if len(s) > 1:
            if s[3] in ['fleishman', 'new-all.q']:
                results[s[3]] += 1
            if s[1] == 'jonatha' and s[2] == 'PEND' and s[3] == 'fleishman':
                pends.append(s[0])
    return results, pends


def move_to_processed(folder: str, running: list, passed: int) -> None:
    """
    deletes folder form the running_folders list, and adds it to the processed list
    :param folder: a folder address that finished being processed
    :param running: a list of running folders
    :return: None
    """
    with open('/home/labs/fleishman/jonathaw/general_lists/processed_folders.txt', 'a') as fout:
        fout.write('%s\t%i\n' % (folder, passed))
    with open('/home/labs/fleishman/jonathaw/general_lists/running_folders.txt', 'w+') as fout:
        for runner in running:
            if runner != folder:
                fout.write('%s\n' % runner)


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
    killed = 0
    for folder_job in folder_jobs:
        if folder_job[4:] in running_jobs.values():
            os.system('bkill %s 1>2>/dev/null' % [k for k, v in running_jobs.items() if v == folder_job[4:]][0])
            killed += 1
    log += 'in %s KILLED %i jobs' % (folder.split('/')[-1], killed)
    os.chdir(pwd)


def get_my_running_jobs() -> dict:
    """
    returns a dict of all my running jobs
    :return: {job_id: job_name}
    """

    results = {}
    proc = subprocess.Popen(['bjobs'], stdout=subprocess.PIPE)
    bjobs = proc.stdout.read()
    t = 0
    for l in str(bjobs).split('\\n'):
        s = l.split()
        if len(s) < 4:
            continue
        if s[1] == 'jonatha' and s[0] != 'wexac':
            try:
                results[s[0]] = s[-4].split('.')[1]
                t += 1
            except:
                pass
    return results


def folder_scores(folder: str) -> dict:
    """
    concatenates all the score files on the folder to one score dict
    :param folder: a folder address
    :return: {name: {filter: grade}} a score dict for the entire folder
    """
    results = {}
    score_files = [a for a in os.listdir(folder) if a[-3:] == '.sc']
    for score in score_files:
        results.update(score2dict(score))
    return results


def is_folder_finished(folder: str) -> bool:
    """
    checks if #outs == #jobs, if it is, checks if thei'r all finished
    :param folder: a folder address
    :return: True iff the folder finished running all jobs, and all of them finished
    """
    global log
    job_files = [a for a in os.listdir(folder) if a[:4] == 'job.']
    out_files = [a for a in os.listdir(folder) if a[:4] == 'out.']
    if len(job_files) != len(out_files):
        log += 'in %s found %i outs, and %i jobs, not finished\n' % (folder.split('/')[-2], len(out_files),
                                                                     len(job_files))
        return False
    for out in out_files:
        with open(out, 'r') as fin:
            cont = fin.read().split('\n')
        if not any(['protocols.jd2.JobDistributor: no more batches to process...' in a for a in cont]):
            log += '%s has enough outs, but NOT FINISHED\n' % folder.split('/')[-2]
            return False
    return True


def submit_folder(folder: str) -> None:
    """
    submits jobs from folder
    :param folder: a folder address
    :return: None
    """
    try:
        os.system('sh %scommand' % folder)
    except:
        pass


def move_pender_to_runner(pending: list, pender: str) -> None:
    """
    removes pender from the pending list, and places it in the running list
    :param pending: list of folders pending for submission
    :param pender: a pender that is being submitted
    :return: None
    """
    with open('/home/labs/fleishman/jonathaw/general_lists/pending_folders.txt', 'w') as fout:
        for p in pending:
            if p != pender:
                fout.write('%s\n' % p)
    with open('/home/labs/fleishman/jonathaw/general_lists/running_folders.txt', 'a') as fout:
        fout.write('%s\n' % pender)


def get_running_folders() -> list:
    """
    :return: list of folders that are currently running
    """
    with open('/home/labs/fleishman/jonathaw/general_lists/running_folders.txt', 'r') as fin:
        cont = fin.read().split('\n')
    resutls = []
    for l in cont:
        if len(l) != 0:
            resutls.append(l if l[-1] == '/' else l + '/')
    return resutls


def get_pending_folders() -> list:
    """
    :return: a list of folders that need to run
    """
    with open('/home/labs/fleishman/jonathaw/general_lists/pending_folders.txt', 'r') as fin:
        cont = fin.read().split('\n')
    resutls = []
    for l in cont:
        if len(l) != 0:
            resutls.append(l if l[-1] == '/' else l + '/')
    return resutls


def lists_status() -> None:
    """
    prints how many folders are in which list
    :return: None
    """
    global log
    for kind in ['pending', 'processed', 'running']:
        num_lines = sum(1 for line in open('/home/labs/fleishman/jonathaw/general_lists/%s_folders.txt' % kind, 'r'))
        log += 'found %i lines in %s\n' % (num_lines, kind)
    num_lines = sum(1 for line in open('/home/labs/fleishman/jonathaw/general_lists/switched_jobs.txt', 'r'))
    log += 'found %i lines in switched_jobs\n' % num_lines


def am_i_running() -> bool:
    """
    :return: whether the .cron_script_running file exists or not
    """
    return os.path.isfile('/home/labs/fleishman/jonathaw/.cron_script_running')


def set_as_running() -> None:
    """
    :return: None. creates  the .cron_script_running file
    """
    os.mknod('/home/labs/fleishman/jonathaw/.cron_script_running')


def set_as_not_running() -> None:
    """
    :return: None. removes the .cron_script_running
    """
    os.remove('/home/labs/fleishman/jonathaw/.cron_script_running')


def should_i_run() -> bool:
    return os.path.isfile('/home/labs/fleishman/jonathaw/.run_lsf_manager')


if __name__ == '__main__':
    global log
    sender = 'jonathan.weinstein2012@gmail.com'
    smtpObj = smtplib.SMTP('localhost')
    while should_i_run():
        set_as_running()
        log = ("From: LSFManager <LSF@manager.com>\n"
           "To: Me <jonathan.weinstein2012@gmail.com>\n"
           "Subject: LSFManager Report\n")
        print('starting run!!!', time.ctime())
        main()
        print('finished run, sending email')
        smtpObj.sendmail(sender, [sender], log)
        set_as_not_running()
        for i in range(61):
            sys.stdout.write('\r')
            sys.stdout.write("[%-60s] %d%%" % ('='*i, 100./60.*float(i)))
            sys.stdout.flush()
            time.sleep(10)
        print('\n')
    print('found I should not run anymore...')
    smtpObj.sendmail(sender, [sender], 'found I should not run anymore...\n')
    # if not am_i_running():
    #     print('cron is not running, so now it is')
    #     set_as_running()
    #     main()
    #     set_as_not_running()
    # else:
    #     print('cron is still running...')
