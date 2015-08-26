def main():
    import timeit
    timer_start = timeit.default_timer()
    PURPLES_THRESHOLD = 15
    with open('/home/labs/fleishman/jonathaw/general_lists/dirs_running_jobs', 'r') as dir_list:
        new_file = open('/home/labs/fleishman/jonathaw/general_lists/dirs_running_jobs_temp', 'wr+')
        for dir_o in dir_list:
            dir_o = dir_o.rstrip()
            print dir_o
            if dir_o[-7:] == 'command':
                dir_o = dir_o[:-7]
            os.chdir(dir_o)
            if is_design(dir_o):
                # if folder_jobs_are_finished(os.getcwd()):
                if is_folder_finished_running(os.getcwd()):
                    subprocess.call(['sh', '/home/labs/fleishman/jonathaw/bin/result_processor.sh'])
            else:
                purples_num = how_many_purples(dir_o)
                if purples_num >= PURPLES_THRESHOLD:
                    print 'was found to have ', purples_num, 'purples, and was stopped'
                    stop_jobs_from_dir(dir_o)
                    subprocess.call(['sh', '/home/labs/fleishman/jonathaw/bin/result_processor.sh'])
                    add_result_to_table(purples_num, dir_o)
                # elif folder_jobs_are_finished(os.getcwd()):
                elif is_folder_finished_running(os.getcwd()):
                    print 'has FINISHED and was processed', purples_num, 'purples were found'
                    subprocess.call(['sh', '/home/labs/fleishman/jonathaw/bin/result_processor.sh'])
                    add_result_to_table(purples_num, dir_o)
                else:
                    print 'had only', purples_num, 'and was not stopped'
                    new_file.write(dir_o + '\n')
    new_file.close()
    shutil.move('/home/labs/fleishman/jonathaw/general_lists/dirs_running_jobs_temp',
                '/home/labs/fleishman/jonathaw/general_lists/dirs_running_jobs')
    import job_submitter
    import running_status
    running_status.my_lsf_status(True)
    job_submitter.job_submitter_main()
    # subprocess.call(['sh', '/home/labs/fleishman/jonathaw/bin/lsfStatus.sh'])
    print 'process took', (timeit.default_timer() - timer_start) / 60, 'minutes'


def score_passes_thresholds(score_line, thresholds, filter_fields):
    score_line = score_line.split()
    if float(score_line[filter_fields['ddg']]) < thresholds['a_ddg']\
            and float(score_line[filter_fields['sasa']]) > thresholds['a_sasa'] \
            and float(score_line[filter_fields['pack']]) > thresholds['a_pack']\
            and float(score_line[filter_fields['shape']]) > thresholds['a_shape']\
            and float(score_line[filter_fields['buried']]) <= thresholds['a_burr']:
        return True
    else:
        return False


def stop_jobs_from_dir(dire):
    os.chdir(dire)
    jobs_in_dir = [f.split('.')[1] for f in os.listdir('.') if re.match('job.*', f)]
    jobs_running = get_my_running_jobs()
    for job in jobs_running.keys():
        if job in jobs_in_dir:
            subprocess.call(['bkill', jobs_running[job]], stdout=open('./killing.log', 'w'),
                            stderr=open('./killing_err.log', 'w'))


def get_my_running_jobs():
    results = {}
    proc = subprocess.Popen(['bjobs'], stdout=subprocess.PIPE)
    bjobs = proc.stdout.read()
    bjobs_split = bjobs.split('\n')
    for bjobs_line in bjobs_split:
        if len(bjobs_line.split()) >= 6:
            try:
                results[(bjobs_line.split()[6].split('.')[1])] = bjobs_line.split()[0]
            except:
                continue
    return results


def folder_jobs_are_finished(cwd):
    os.chdir(cwd)
    running_jobs = get_my_running_jobs()
    jobs_in_dir = [f.split('.')[1] for f in os.listdir('.') if re.match('job.*', f)]
    intercept = [x for x in running_jobs.keys() if x in jobs_in_dir]
    return True if len(intercept) == 0 else False


def add_result_to_table(purple_num, name):
    with open('/home/labs/fleishman/jonathaw/general_lists/score_lists', 'a') as f:
        f.write(name + '\t' + str(purple_num) + '\n')


def how_many_purples(dir_i):
    import os
    import re
    cwd = os.getcwd()
    os.chdir(dir_i)
    thresholds = {'a_ddg': -20.0, 'a_sasa': 1300.0, 'a_pack': 0.6, 'a_shape': 0.5, 'a_burr': 1.0}
    purples_file = open('./purples_script.score', 'w')
    out_files = [f for f in os.listdir('.') if re.match('.*\.out', f)]
    filter_fields = {}
    purples_num = 0
    add_header = True
    for file_name in out_files:
        score_file = open(file_name, 'r')
        for file_line in score_file:
            line_split = file_line.split()
            if 'SCORE:' in line_split:
                if line_split[1] == 'score':
                    if add_header:
                        purples_file.write(file_line.rstrip() + '\n')
                        add_header = False
                    filter_fields['ddg'] = line_split.index('a_ddg')
                    filter_fields['sasa'] = line_split.index('a_sasa')
                    filter_fields['pack'] = line_split.index('a_packstat')
                    filter_fields['shape'] = line_split.index('a_shape')
                    filter_fields['buried'] = line_split.index('a_buried_2')
                else:
                    if score_passes_thresholds(file_line, thresholds, filter_fields):
                        purples_file.write(file_line.rstrip() + '\n')
                        purples_num += 1
        score_file.close()
    purples_file.close()
    os.chdir(cwd)
    return purples_num


def is_design(dir_i):
    with open('/home/labs/fleishman/jonathaw/general_lists/design_folders', 'r') as f:
        for line in f:
            if dir_i == line:
                return True
    return False


def is_folder_finished_running(folder):
    import re
    import os
    finished_out_num = 0
    cwd = os.getcwd()
    os.chdir(folder)
    out_file_list = [x for x in os.listdir(folder)
                       if re.match('out.*', x)]
    for out_file in out_file_list:
        for line in open(out_file):
            pass
        if line.split()[0] == 'protocols.jd2.JobDistributor:' and \
                        line.split()[2] == 'jobs' and line.split()[3] == 'considered,':
            finished_out_num += 1
    job_file_list = [x for x in os.listdir(folder) if re.match('job.*', x)]
    job_num = len(job_file_list)
    print 'found %i job files and %i finished out files' % (job_num, finished_out_num)
    print 'declared %s' % ('finished' if job_num == finished_out_num else 'NOT finished')
    os.chdir(cwd)
    return True if job_num == finished_out_num else False

if __name__ == '__main__':
    import os
    import re
    import subprocess
    import shutil
    main()
