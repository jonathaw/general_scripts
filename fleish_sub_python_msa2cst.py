#!/usr/bin/env python

def main():
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument('-name')
    parser.add_argument('-msa_percentile')
    parser.add_argument('-hp_threshold')
    parser.add_argument('-z_0')
    parser.add_argument('-w')
    args = vars(parser.parse_args())
    args['name'] = args['name'].lower()
    fout = open('job.'+args['name'], 'wr+')
    command = open('command', 'a+')
    fout.write('#!/bin/bash\n')
    fout.write('. /usr/share/lsf/conf/profile.lsf\n')
    fout.write('cd %s\n' % os.getcwd())
    fout.write('/apps/RH6U4/python/2.7.6/bin/python /home/labs/fleishman/jonathaw/membrane_prediciton/TMpredict_WinGrade.py '
               '-name %s -with_msa True -msa_percentile %s -w %s -z_0 %s -hp_threshold %s \n' % (args['name'], args['msa_percentile'],
                                                                               args['w'], args['z_0'],args['hp_threshold']))
    fout.write('/apps/RH6U4/python/2.7.6/bin/python /home/labs/fleishman/jonathaw/membrane_prediciton/TMConstraint.py '
               '-mode pred2cst -name %s\n' % args['name'])
    fout.write('/apps/RH6U4/python/2.7.6/bin/python /home/labs/fleishman/jonathaw/membrane_prediciton/TMpredict_WinGrade.py '
               '-name %s -with_cst True -w %s -z_0 %s -hp_threshold 6 \n' % (args['name'], args['w'], args['z_0']))

    command.write('bsub -u /dev/null -G fleishman-wx-grp-lsf -q new-all.q -o %s/out.%s -e %s/err.%s %s/job.%s /apps/RH6U4/blcr/0.8.5/bin/cr_run\n' % (os.getcwd(), args['name'], os.getcwd(), args['name'], os.getcwd(), args['name']))
    os.system('chmod +x %s/job.%s' % (os.getcwd(), args['name']))#-R rusage[mem=4096]
    fout.close()
# -M 4194304
# -R rusage"[mem=2048]"
if __name__ == '__main__':
    main()