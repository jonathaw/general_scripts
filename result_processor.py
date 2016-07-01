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


def process_folder(args):
    import time
    pwd = os.getcwd()+'/'
    os.chdir(args['folder'])

    name = pwd.split('/')[-2]
    time = time.strftime("%d.%0-m")

    if not args['force_process']:
        if os.path.isfile('%sall_%s_%s.err' % (pwd, name, time)):
            print('found %sall_%s_%s.err, STOPPING' % (pwd, name, time))
            if args['remove_pdbs']:
                remove_pdbs_only(pwd, name, time)
                sys.exit()
            return 'not finished'
        if not is_folder_finished(pwd):
            return 'not finished'
    sc_files = [a for a in os.listdir(pwd) if a[-3:] == '.sc']
    pdb_files = [a for a in os.listdir(pwd) if a[-7:] == '.pdb.gz']
    err_files = [a for a in os.listdir(pwd) if a[:4] == 'err.']
    job_files = [a for a in os.listdir(pwd) if a[:4] == 'job.']
    print('found a total of %i job files, %i err files, %i pdbs and %i scores' % (len(job_files), len(err_files),
                                                                                  len(pdb_files), len(sc_files)))
    if len(sc_files) == 0:
        return 'no scores'
    combine_scores('%sall_%s_%s.score' % (pwd, name, time), sc_files)
    non_triv_errs = process_errors('%sall_%s_%s.err' % (pwd, name, time), err_files)

    if non_triv_errs == 0:
        print('removing out.* and job.*')
        [os.remove(out) for out in os.listdir(pwd) if out[:4] == 'out.']
        [os.remove(job) for job in job_files]
        try:
            os.remove('./command')
        except:
            pass

    run_filters = DCRP.generate_run_filters(args={'ddg': 24.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6,
                                                  'buried_2': 3, 'hbonds': -10.})
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
    [os.remove(pdb) for pdb in pdb_files if pdb[:-7] not in best_names]
    os.chdir(pwd)


def remove_pdbs_only(pwd, name, time):
    pdb_files = [a for a in os.listdir(pwd) if a[-7:] == '.pdb.gz']
    run_filters = DCRP.generate_run_filters(args={'ddg': 24.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6,
                                                  'buried_2': 3, 'hbonds': -10.})
    score_dict = RF.score2dict('%sall_%s_%s.score' % (pwd, name, time))
    passed, failed = DCRP.all_who_pass_run_filters({}, score_dict, run_filters)
    if len(passed) != 0:
        print('there are %i passed scores, so choosing from there' % len(list(passed.keys())))
        best_structs = DCRP.best_n_structures({'filter': 'ddg', 'n': min([10, len(list(passed.keys()))])}, passed)
    else:
        print('there are no passed, so choosing from the failed')
        best_structs = DCRP.best_n_structures({'filter': 'ddg', 'n': min([5, len(list(failed.keys()))])}, failed)
    best_names = [a['description'] for a in best_structs]
    [os.remove(pdb) for pdb in pdb_files if pdb[:-7] not in best_names]


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


def is_folder_finished(folder):
    out_files = [a for a in os.listdir(folder) if a[:4] == 'out.']
    for out in out_files:
        with open(out, 'r') as fin:
            cont = fin.read().split('\n')
        if not any(['protocols.jd2.JobDistributor: no more batches to process...' in a for a in cont]):
            print('NOT FINISHED', out)
            return False
    return True

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-save_pass', type=bool, default=False)
    parser.add_argument('-remove_pdbs', type=bool, default=False, help='if script ran over folder but did not erase the '
                                                                                                        'pdbs, use this')
    parser.add_argument('-force_process', type=bool, default=False)
    parser.add_argument('-folder', default=os.getcwd()+'/')
    args = vars(parser.parse_args())

    process_folder(args)