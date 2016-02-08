#!/usr/bin/env python3.5
"""
a bundle of scripts to create and run Elazar membrane energy function benchmarks
"""
import argparse
import os
import sys
import time
import shutil
import subprocess

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from collections import OrderedDict
from Equation import Expression

from RosettaFilter import score2dict, score_dict2df, df2boxplots

lsf_username = 'jonatha'
pwd = os.getcwd()+'/'
all_configurations_path = pwd+'all_configs/'

rosetta_executables_path = '/home/labs/fleishman/jonathaw/Rosetta/main/source/bin/'
rosetta_script_exec_path = 'rosetta_scripts.default.linuxgccrelease'
rosetta_database_path = '/home/labs/fleishman/jonathaw/Rosetta/main/database/'
protocols_path = '/home/labs/fleishman/jonathaw/RMP_temp/protocols/'
polyval_file_path = protocols_path+'polyval_21_5_15.txt'
MPMutateRelax = 'MPMutateRelax.xml'

num_As = 26
name = 'polyA'
nstruct = 5
membrane_half_depth = 15
lazaridis_poly_deg = 4

aas_names = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG',
             'SER', 'THR', 'VAL', 'TRP', 'TYR']
aas_3_1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K',
           'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
           'TRP': 'W', 'TYR': 'Y'}
aas_1_3 = {v: k for k, v in aas_3_1.items()}

###
polyA_total_mean = -45.51
###

def main():
    global polyA_total_mean
    create_polyA_fasta()
    sequence_to_idealized_helix()
    create_spanfile()
    trunctate_2nd_mem_res()
    polyA_bench_sc = run_polyA_score_benchmark()
    polyA_total_mean = np.mean(polyA_bench_sc['score'].tolist())
    polyA_total_std = np.std(polyA_bench_sc['score'].tolist())
    print('average polyA total score %.3f with std %.3f' % (polyA_total_mean, polyA_total_std))
    run_all_configurations()


def run_all_configurations():
    """
    :return:
    """
    os.mkdir(pwd+'all_configs')
    os.chdir(pwd+'all_configs')
    make_all_configurations_jobs()
    os.system('sh command')
    sleep_until_jobs_finish()
    create_profiles()
    os.chdir(pwd)


def create_profiles():
    profiles = {aa: {pos: None for pos in range(1, num_As+1)} for aa in aas_names}
    sc_files = [a for a in os.listdir(all_configurations_path) if '.sc' in a]
    for sc_file in sc_files:
        sc = score2dict(all_configurations_path+sc_file)
        df = score_dict2df(sc)
        aa = [a for a in aas_names if a in sc_file][0]
        pos = int(sc_file.split('_')[2].split('.')[0])
        profiles[aa][pos] = np.mean(df['score']) - polyA_total_mean
    membrane_position = np.linspace(-membrane_half_depth, membrane_half_depth, endpoint=True, num=num_As+1)
    mm_profiles = {aa: OrderedDict((membrane_position[pos], profiles[aa][pos]) for pos in range(1, num_As+1))
                   for aa in aas_names}
    mm_polys = {aa: np.polyfit(list(mm_profiles[aa].keys()), list(mm_profiles[aa].values()), lazaridis_poly_deg)
                for aa in aas_names}
    mm_eqs = {aa: Expression('+'.join('%f*z^%i' % (n, i) for i, n in enumerate(mm_polys[aa][::-1])),
                             argorder=['z']) for aa, val in mm_polys.items()}
    draw_profiles(mm_profiles, mm_eqs)


def get_elazar_scale():
    resutls = {}
    for res, val in MakeHydrophobicityGrade().items():
        resutls[res] = Expression("%f*z^4+%f*z^3+%f*z^2+%f*z+%f" % (val[0], val[1], val[2], val[3], val[4]),
                                  argorder=['z'])
    return resutls


def MakeHydrophobicityGrade():
    """
    :return: returns a dictionary of the polynom values for each residue
    """
    global hydrophobicity_polyval
    hydrophobicity_grade = open(polyval_file_path, 'r')
    hydrophobicity_polyval = {}
    for line in hydrophobicity_grade:
        split = line.split()
        hydrophobicity_polyval[aas_1_3[split[0]].upper()] = [float(n) for n in split[1:6]]
    hydrophobicity_grade.close()
    return hydrophobicity_polyval


def draw_profiles(mm_profiles: dict, mm_eqs: dict):
    fig = plt.figure()
    elazar_profiles = get_elazar_scale()
    z = np.linspace(-membrane_half_depth, membrane_half_depth, endpoint=True, num=num_As+1)
    for i, (aa, val) in enumerate(mm_profiles.items()):
        ax = plt.subplot(5, 4, i+1)
        # Lazaridis - line
        plt.plot([k for k in val.keys()], [v for v in val.values()], label='Lazaridis', color='k')

        # Lazaridis - polyfit
        plt.plot(z, [mm_eqs[aa](z_) for z_ in z], label='Lazaridis-fit', color='r')

        # Elazar fit
        plt.plot(z, [elazar_profiles[aa](z_) for z_ in z], label='Elazar', color='b')

        # plot attributes
        plt.title(aa)
        plt.xlim([-(membrane_half_depth+0.5), membrane_half_depth+0.5])
        plt.ylim([-3., 3.])#-5., 6.
    ax.legend(bbox_to_anchor=(1.05, 0), loc='lower left', borderaxespad=0.)
    plt.show()


def sleep_until_jobs_finish() -> None:
    while True:
        bjobs = str(subprocess.Popen(['bjobs'], stderr=subprocess.PIPE, stdout=subprocess.PIPE).stderr.read())
        if 'No unfinished job found' in bjobs:
            print('all jobs finished')
            break
        time.sleep(5)


def make_all_configurations_jobs():
    shutil.copy(pwd+name+'.span', './')
    shutil.copy(pwd+name+'.pdb', './')
    for i in range(1, num_As+1):
        for aa in aas_names:
            job_args = {
                    '-parser:protocol': '%s%s' % (protocols_path, MPMutateRelax),
                    '-database': '/home/labs/fleishman/jonathaw/Rosetta/main/database/',
                    '-overwrite': None,
                    '-s': '%s.pdb' % name,
                    '-out:suffix': '_%s_%i.pdb' % (aa, i),
                    '-mp:setup:spanfiles': name + '.span',
                    '-nstruct': nstruct,
                    '-mp:scoring:hbond': None,
                }
            script_vars = {'mut_pos=': i, 'mut_iden=': aa}
            make_jobs(aa, i, job_args, script_vars)
    print('created all jobs to run for all mutations in all %i positions' % num_As)


def run_polyA_score_benchmark() -> pd.DataFrame:
    os.system('%s%s -database %s -parser:protocol %s -s %s -mp:setup:spanfiles %s -nstruct %i -overwrite '
              '-mp:scoring:hbond '
              '-parser:script_vars mut_pos=1 mut_iden=ALA' %
              (rosetta_executables_path, rosetta_script_exec_path, rosetta_database_path,
               protocols_path+MPMutateRelax, pwd+name+'.pdb', pwd+name+'.span', nstruct))
    shutil.move('%sscore.sc' % pwd, '%spolyA_benchmark.score' % pwd)
    print('ran a FastRelax benchmark and placed the score in polyA_benchmark.score')
    sc_dict = score2dict('%spolyA_benchmark.score' % pwd)
    return score_dict2df(sc_dict)


def create_polyA_fasta() -> None:
    """
    creates a fasta file with a num_As polyA in it
    """
    with open(pwd+name+'.fa', 'w+') as fout:
        fout.write('>polyA\n%s\n' % ''.join(['A']*num_As))
        print('created polyA file at %s.fa with %i As' % (pwd+name, num_As))


def make_jobs(aa: str, pos: int, job_args: dict, script_vars: dict=None) -> None:
    """
    :param aa: aa type
    :param pos: helix sequence position
    :param job_args: arguments for Rosetta
    :param script_vars: script vars (mut_iden and mut_pos)
    :return: creates the jobs that will be the benchmark
    """
    name = 'cnfg_%i_%s' % (pos, aa)
    jobname = '%sjob.%s' % (all_configurations_path, name)
    outname = '%sout.%s' % (all_configurations_path, name)
    errname = '%serr.%s' % (all_configurations_path, name)
    cmdname = '%scommand' % all_configurations_path
    with open(jobname, 'w+') as job:
        job.write('#!/bin/bash\n')
        job.write('. /usr/share/lsf/conf/profile.lsf\n')
        job.write('cd ' + all_configurations_path + '\n')
        if args['queue'] == 'new-all.q':
            job.write(
                '/apps/RH6U4/blcr/0.8.5/bin/cr_run '
                '/home/labs/fleishman/jonathaw/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease ')
        else:
            job.write('/home/labs/fleishman/jonathaw/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease ')
        for k, v in job_args.items():
            if v is None:
                job.write('%s ' % k)
            else:
                job.write('%s %s ' % (k, v))
        if script_vars is not None:
            job.write('-parser:script_vars ')
            for k, v in script_vars.items():
                job.write('%s%s ' % (k, v))
        job.write('\n')
    with open(cmdname, 'a+') as cmd:
        if args['queue'] == 'new-all.q':
            cmd.write(str(
                'bsub -C 1024 -u /dev/null -N -u /dev/null -R rusage[mem=1024] -L /bin/bash -G fleishman-wx-grp-lsf -q ' +
                args['queue'] + ' -o ' + outname + ' -e ' + errname + ' /apps/RH6U4/blcr/0.8.5/bin/cr_run ' +
                jobname + '\n'))
        else:
            cmd.write(str('bsub -L /bin/bash -N -u /dev/null -G fleishman-wx-grp-lsf -q ' + args['queue'] + ' -o ' +
                          outname + ' -e ' + errname + ' ' + jobname + '\n'))
    os.system('chmod +x %s' % jobname)


def sequence_to_idealized_helix() -> None:
    """
    calls the Rosetta application that turns a sequence into a membrane embedded helix pdb
    """
    os.system('%shelix_from_sequence.default.linuxgccrelease -in:file:fasta %s '
                    '-mp:setup:transform_into_membrane 1' % (rosetta_executables_path, pwd+name+'.fa'))
    shutil.move(pwd+'helix_from_sequence.pdb', pwd+name+'.pdb')
    print('created an idealised helix from %s and put it in %s' % (name, name+'.pdb'))


def create_spanfile() -> None:
    """
    creates a simple spanfile for the polyA pdb. basically all 1-num_As residues are in the membrane
    """
    with open('%s%s.span' % (pwd, name), 'w+') as fout:
        fout.write('Rosetta-generated spanfile from SpanningTopology object\n1 %i\nantiparallel\nn2c\n\t1	%i\n'
                   % (num_As, num_As))


def trunctate_2nd_mem_res() -> None:
    """
    truncate the 2nd membrane residue from the helix pdb. for some reason Rosetta doesn't like it
    """
    with open('%s%s.pdb' % (pwd, name), 'r') as fin:
        with open('%s%s.tmp' % (pwd, name), 'w+') as fout:
            for l in fin.read().split('\n'):
                if 'XXXX' not in l:
                    fout.write(l+'\n')
    shutil.move('%s%s.tmp' % (pwd, name), '%s%s.pdb' % (pwd, name))
    print('removed the second membrane residue')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode')
    parser.add_argument('-queue', default='new-all.q')

    args = vars(parser.parse_args())

    if args['mode'] == 'main':
        main()

    elif args['mode'] == 'create_polyA':
        create_polyA_fasta()

    elif args['mode'] == 'sequence_to_idealized_helix':
        sequence_to_idealized_helix()

    elif args['mode'] == 'create_spanfile':
        create_spanfile()

    elif args['mode'] == 'trunctate_2nd_mem_res':
        trunctate_2nd_mem_res()

    elif args['mode'] == 'profiles':
        create_profiles()

    elif args['mode'] == 'sleep':
        sleep_until_jobs_finish()

    elif args['mode'] == 'test':
        create_profiles()

    else:
        print('unknown mode')