import os
import re
import subprocess
import shutil

with_burried = True
PURPLES_THRESHOLD = 15
thresholds = {'a_ddg': -20, 'a_sasa': 1300, 'a_pack': 0.6, 'a_shape': 0.5, 'a_burr': 1.0}
# thresholds = {'a_ddg': -0, 'a_sasa': 10, 'a_pack': 0, 'a_shape': 0, 'a_burr': 1}


def score_passes_thresholds(score_line):
    score_line = score_line.split()
    if float(score_line[ddg_field]) <= thresholds['a_ddg']\
            and float(score_line[sasa_field]) >= thresholds['a_sasa'] \
            and float(score_line[pack_field]) >= thresholds['a_pack']\
            and float(score_line[shape_field]) >= thresholds['a_shape']:
        if with_burried and float(score_line[buried_field]) > thresholds['a_burr']:
            return False
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


with open('/home/labs/fleishman/jonathaw/bin/dirs_running_jobs', 'r') as dir_list:
    new_file = open('/home/labs/fleishman/jonathaw/bin/dirs_running_jobs_temp', 'wr+')
    for line in dir_list:
        os.chdir(line.rstrip())
        purples_file = open('./purples_script.score', 'w')
        out_files = [f for f in os.listdir('.') if re.match('.*\.out', f)]
        purples_num = 0
        add_header = True
        for file_name in out_files:
            score_file = open(file_name, 'r')
            for file_line in score_file:
                line_split = file_line.split()
                if 'SCORE:' in line_split:
                    if line_split[1] == 'score':
                        if add_header:
                            purples_file.write(file_line + '\n')
                        ddg_field = line_split.index('a_ddg')
                        sasa_field = line_split.index('a_sasa')
                        pack_field = line_split.index('a_packstat')
                        shape_field = line_split.index('a_shape')
                        if with_burried:
                            buried_field = line_split.index('a_buried_2')
                    else:
		    	#print 'in', line, 'and', score_passes_thresholds(file_line)
                        if score_passes_thresholds(file_line):
                            purples_file.write(file_line)
                            purples_num += 1
			    #print 'now has', purples_num
            score_file.close()
        if purples_num > PURPLES_THRESHOLD:
            stop_jobs_from_dir(line.rstrip())
            subprocess.call(['sh', '/home/labs/fleishman/jonathaw/bin/result_processor.sh'])
            print line, 'was found to have ', purples_num, 'purples, and was stopped'
        else:
            print line.rstrip(), 'had only', purples_num, 'and was not stopped'
            new_file.write(line)

        purples_file.close()
new_file.close()
shutil.move('/home/labs/fleishman/jonathaw/bin/dirs_running_jobs_temp', '/home/labs/fleishman/jonathaw/bin/dirs_running_jobs')
