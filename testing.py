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
    print job_num, finished_out_num
    os.chdir(cwd)
    return True if job_num == finished_out_num else False

if __name__ == '__main__':
    import os
    print is_folder_finished_running(os.getcwd())
