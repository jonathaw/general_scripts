#!/usr/bin/env python3.5
"""
process rosetta run results
"""
import os
import sys
import argparse
import subprocess
import RosettaFilter as RF
import DoCohResultProcessor as DCRP


def main():
    import time
    parser = argparse.ArgumentParser()
    parser.add_argument('-save_pass', type=bool, default=False)
    args = vars(parser.parse_args())
    pwd = os.getcwd()+'/'
    name = pwd.split('/')[-2]
    time = time.strftime("%d.%0-m")

    if os.path.isfile('%sall_%s_%s.err' % (pwd, name, time)) or os.path.isfile('%sall_%s_%s.score' % (pwd, name, time)):
        print('found %sall_%s_%s.err / score, STOPPING' % (pwd, name, time))
        sys.exit()

    sc_files = [a for a in os.listdir(pwd) if a[-3:] == '.sc']
    pdb_files = [a for a in os.listdir(pwd) if a[-4:] == '.pdb']
    err_files = [a for a in os.listdir(pwd) if a[:4] == 'err.']
    job_files = [a for a in os.listdir(pwd) if a[:4] == 'job.']
    print('found a total of %i job files, %i err files, %i pdbs and %i scores' % (len(job_files), len(err_files),
                                                                                  len(pdb_files), len(sc_files)))
    if len(sc_files) == 0:
        print('found 0 .sc files, so EXITING')
        sys.exit()
    combine_scores('%sall_%s_%s.score' % (pwd, name, time), sc_files)
    non_triv_errs = process_errors('%sall_%s_%s.err' % (pwd, name, time), err_files)

    if non_triv_errs == 0:
        print('removing out.* and job.*')
        [os.remove(out) for out in os.listdir(pwd) if out[:4] == 'out.']
        [os.remove(job) for job in job_files]
        os.remove('./command')

    run_filters = DCRP.generate_run_filters()
    score_dict = RF.score2dict('%sall_%s_%s.score' % (pwd, name, time))
    passed, failed = DCRP.all_who_pass_run_filters({}, score_dict, run_filters)
    if len(passed) != 0:
        print('there are %i passed scores, so choosing from there' % len(list(passed.keys())))
        best_structs = DCRP.best_n_structures({'filter': 'ddg', 'n': min([10, len(list(passed.keys()))])}, passed)
    else:
        print('there are no passed, so choosing from the failed')
        best_structs = DCRP.best_n_structures({'filter': 'ddg', 'n': min([10, len(list(failed.keys()))])}, failed)
    best_names = [a['description'] for a in best_structs]
    print('the best:', best_names)
    print('removing all other pdbs')
    if not args['save_pass']:
        [os.remove(pdb) for pdb in pdb_files if pdb[:-4] not in best_names]
    else:
        pass_names = list(passed.keys())
        [os.remove(pdb) for pdb in pdb_files if pdb[:-4] not in pass_names]


def process_errors(file_name, err_list):
    non_triv = []
    with open(file_name, 'w+') as err_out:
        for err in err_list:
            only_non_triv = True
            with open(err, 'r') as err_in:
                for l in err_in.read().split('\n'):
                    if 'load in l':
                        continue
                    else:
                        err_out.write(l+'\n')
                        only_non_triv = False
                        non_triv.append(l)
            if only_non_triv:
                os.remove(err)
    print('found %i non trivial errors, in %i files' % (len(non_triv), len(err_list)))
    return len(non_triv)


def combine_scores(file_name, sc_list):
    header = subprocess.check_output('grep description %s' % sc_list[0], shell=True).decode()
    cont = subprocess.check_output('grep SCORE: *.sc | grep -v description', shell=True)
    with open(file_name, 'w+') as fout:
        fout.write(str(header))
        for l in str(cont).split('\\n'):
            fout.write(str(l)+'\n')
    if len(cont) > 100:
        print('found many scores, removing *sc')
        [os.remove(sc) for sc in sc_list]


if __name__ == '__main__':
    main()