def my_lsf_status(email):
    import subprocess
    import os.path
    from whose_running_2 import lsf_user
    from my_colors import my_colors
    my_colors_ = my_colors()
    global QUEUES, MODES
    QUEUES = ['fleishman', 'centos', 'new-all.q']
    MODES = ['RUN', 'PEND', 'SSUSP', 'UNKWN', 'PSUSP', 'USUSP']
    proc = subprocess.Popen(['bjobs', '-u', 'all'], stdout=subprocess.PIPE)
    bjobs = proc.stdout.read()
    bjobs_split = bjobs.split('\n')
    jonathaw = lsf_user('jonatha')
    for line in bjobs_split:
        line_split = line.split()
        if len(line_split) > 3 and line_split[3] in QUEUES and line_split[1] == 'jonatha':
            jonathaw.add_to_count(line_split[3], line_split[2])
    print '%-9s' % '' + ''.join("%-30s" % queue for queue in QUEUES)
    print '%-9s' % 'NAME' + ''.join(''.join("%-6s"  % mode for mode in MODES) for queue in QUEUES)
    if email:
        print jonathaw.print_results()
    else:
        print my_colors_['red'].format(jonathaw.print_results())
    dirs = sum(1 for line in open('/home/labs/fleishman/jonathaw/general_lists/dirs_running_jobs', 'r'))
    coms = sum(1 for line in open('/home/labs/fleishman/jonathaw/bin/command_list', 'r'))
    cron = os.path.isfile('/home/labs/fleishman/jonathaw/.cron_script_running')
    if email:
        print "dirs to check %s\t\tcoms to submit %s\tcron is running %s" % (dirs, coms, cron)
    else:
        print my_colors_['green'].format("dirs to check %s\t\tcoms to submit %s\tcron is running %s" % (dirs, coms, cron))


if __name__ == '__main__':
    my_lsf_status(False)
