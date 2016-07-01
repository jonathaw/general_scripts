def job_submitter_main():
    import subprocess
    import os
    running = jobs_parser()
    if running['fleishman'] <= 7000:
        cmd = get_command_name()
        print 'submitting %s to fleishman where %i jobs are running' % (cmd, running['fleishman'])
        subprocess.call(['sh', cmd], stdout=open(os.devnull, 'w'))
        place_cmd_in_list(cmd)
    # if running['centos'] <= 200:
    #     cmd = get_command_name().lstrip()
    #     print 'submitting %s to centos where %i jobs are running' % (cmd, running['centos'])
    #     convert_command_to_centos(cmd)
    #     subprocess.call(['sh', cmd], stdout=open(os.devnull, 'w'))
    #     place_cmd_in_list(cmd)
    #     with open('/home/labs/fleishman/jonathaw/general_lists/commands_converted_to_centos', 'a') as t:
    #         t.write(cmd)


def jobs_parser():
    import subprocess
    proc = subprocess.Popen(['bjobs'], stdout=subprocess.PIPE)
    bjobs = proc.stdout.read()
    bjobs_split = bjobs.split('\n')
    results = {'fleishman': 0, 'centos': 0}
    for bjobs_line in bjobs_split:
        bjobs_line_split = bjobs_line.split()
        if 'fleishman' in bjobs_line_split or 'centos' in bjobs_line_split:
            results[bjobs_line_split[3]] += 1
    return results


def convert_command_to_centos(cmd_i):
    import re
    import shutil
    with open(cmd_i.lstrip(), 'r') as cmd_in:
        with open('/home/labs/fleishman/jonathaw/bin/temp_command_file', 'w') as cmd_out:
            for cmd_line in cmd_in:
                # print cmd_line
                cmd_line = re.sub(r' fleishman ', r' centos ', cmd_line)
                cmd_out.write(cmd_line)
    shutil.move('/home/labs/fleishman/jonathaw/bin/temp_command_file', cmd_i)


def get_command_name():
    import shutil
    import sys
    first = True
    result = ''
    with open('/home/labs/fleishman/jonathaw/bin/command_list', 'r') as cmd_in:
        with open('/home/labs/fleishman/jonathaw/bin/command_list_temp', 'w') as cmd_out:
    # with open('/home/labs/fleishman/jonathaw/temp_lsf_manager/temp_dirs', 'r') as cmd_in:
    #     with open('/home/labs/fleishman/jonathaw/temp_lsf_manager/temp_dirs_temp', 'w') as cmd_out:
            for line in cmd_in:
                if first:
                    result = line.rstrip()
                    first = False
                else:
                    cmd_out.write(line)
    if len(result) == 0:
        sys.exit()
    shutil.move('/home/labs/fleishman/jonathaw/bin/command_list_temp',
                '/home/labs/fleishman/jonathaw/bin/command_list')
    # shutil.move('/home/labs/fleishman/jonathaw/temp_lsf_manager/temp_dirs_temp',
    #             '/home/labs/fleishman/jonathaw/temp_lsf_manager/temp_dirs')
    return result


def place_cmd_in_list(cmd_i):
    with open('/home/labs/fleishman/jonathaw/general_lists/dirs_running_jobs', 'a') as f:
        f.write(cmd_i.rstrip() + '\n')


if __name__ == '__main__':
    job_submitter_main()