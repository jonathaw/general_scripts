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
import matplotlib as mpl
import pandas as pd
import numpy as np

from collections import OrderedDict
from Equation import Expression


from RosettaFilter import score2dict, score_dict2df, df2boxplots
from Logger import Logger

# mpl.rc_context(fname="/home/labs/fleishman/jonathaw/.matplotlib/publishable_matplotlibrc")
lsf_username = 'jonatha'
pwd = os.getcwd()+'/'
all_configurations_path = pwd+'all_configs/'

rosetta_executables_path = '/home/labs/fleishman/jonathaw/Rosetta/main/source/bin/'
rosetta_script_exec_path = 'rosetta_scripts.default.linuxgccrelease'
rosetta_database_path = '/home/labs/fleishman/jonathaw/Rosetta/main/database/'
rosetta_membranepotential_path = rosetta_database_path + 'scoring/score_functions/MembranePotential/'
protocols_path = '/home/labs/fleishman/jonathaw/elazaridis/protocols/'
polyval_file_path = '/home/labs/fleishman/jonathaw/membrane_prediciton/mother_fucker_MET_LUE_VAL_sym_k_neg.txt' #'''polyval_21_5_15.txt'
MPMutateRelax = 'MPMutateRelax.xml'
MPFilterScan = 'MPFilterScan_ELazaridis.xml'

num_As = 26
name = 'polyA'
nstruct = 1 #0
membrane_half_depth = 15
lazaridis_poly_deg = 4
flank_size = 6
Z = np.linspace(-membrane_half_depth, membrane_half_depth, num=num_As)
pos_z_dict = {i+1: z for i, z in enumerate(Z)}

aas_names = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG',
             'SER', 'THR', 'VAL', 'TRP', 'TYR']
aas_3_1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K',
           'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
           'TRP': 'W', 'TYR': 'Y'}
aas_1_3 = {v: k for k, v in aas_3_1.items()}
AAs = list('ACDEFGHIKLMNPQRSTVWY')

color_map = {'full': 'blue', 'no_Menv': 'grey', 'elazar': 'red', 'ResSolv': 'purple', 'diff_ips': 'orange'}
rmsd_iter = {aa: [] for aa in AAs}


class PoZEnergy:
    def __init__(self, pos: int, z: float, energy: float):
        self.pos = pos
        self.z = z
        self.energy = energy

    def __repr__(self):
        return '~%i %.1f %.1f~' % (self.pos, self.z, self.energy)

    def __sub__(self, other):
        """
        changes self.energy to self.energy - other.energy
        """
        self.energy -= other.energy


class InsertionProfile:
    def __init__(self, aa: str, pos_score: dict, membrane_half_depth: int=membrane_half_depth, residue_num: int=num_As):
        """

        """
        self.AA = aa
        self.membrane_half_depth = membrane_half_depth
        self.residue_num = residue_num
        self.pos_score = pos_score
        self.poz_energy = pos_energy_dict_to_PoZEnergy_list(pos_score)
        self.polynom = np.polyfit(Z, [x.energy for x in self.poz_energy], 4)

    def polynom_at_z(self, z: float) -> float:
        """
        returns the polynoms value at z
        """
        return np.polyval(self.polynom, z)

    def format_polyval(self):
        """
        print the polyval values for table
        """
        return ' '.join([str(a) for a in self.polynom])

    def rmsd_ips(self, other) -> float:
        """
        retruns the difference between the IPs calculated as by RMSD over Z
        """
        res = 0.0
        for z in Z:
            res += (self.polynom_at_z(z) - other.polynom_at_z(z))**2
        return np.sqrt(np.mean(res))


def pos_energy_dict_to_PoZEnergy_list(pos_energy_dict: dict) -> list():
    """
    creates an ordered list of PoZEnergy instances corresponding to their positions
    """
    result = []
    for pos in range(1, num_As+1):
        result.append(PoZEnergy(pos, pos_z_dict[pos], pos_energy_dict[pos]))
    return result


def subtract_IP_from_IP(ip1: InsertionProfile, ip2: InsertionProfile) -> InsertionProfile:
    new_pos_score = {k: v - ip2.pos_score[k] for k, v in ip1.pos_score.items()}
    return InsertionProfile(ip1.AA, new_pos_score)


def main():
    # create files for running benchmark
    if args['full']:
        create_polyA_fasta()
        sequence_to_idealized_helix()
        create_spanfile()
        trunctate_2nd_mem_res()

    elazar_ips = create_elazar_ips()

    # empty IPs for later
    MPResSolv_remade_ips = {aa: InsertionProfile(aa, pos_score={x: 0 for x in range(1, num_As+1)}) for aa in AAs}

    # first FilterScan run. using null ResSolv
    all_cnfgs_df = filterscan_analysis()
    full_ips = create_insertion_profiles(all_cnfgs_df, 'full_normed')
    noMenv_ips = create_insertion_profiles(all_cnfgs_df, 'noMenv_normed')
    MPResSolv_ips = create_insertion_profiles(all_cnfgs_df, 'ResSolv_normed')

    # calc the difference InsertionProfiles between Elazar and Rosetta. assign them as the polynom table
    diff_ips = {k: subtract_IP_from_IP(elazar_ips[k], noMenv_ips[k]) for k in AAs}
    create_polyval_table(diff_ips, 'ELazaridis_0.txt')

    # analyse Rosetta again.
    all_cnfgs_df = filterscan_analysis()
    MPResSolv_current_ips = create_insertion_profiles(all_cnfgs_df, 'ResSolv_normed')


    # if args['show_fig']:
    #     draw_filterscan_profiles(OrderedDict({'full': full_ips, 'no_Menv': noMenv_ips, 'ResSolv': MPResSolv_ips,
    #                                           'diff_ips': diff_ips, 'elazar': elazar_ips}))

    for aa in AAs:
        rmsd = elazar_ips[aa].rmsd_ips(MPResSolv_current_ips[aa])
        rmsd_iter[aa].append(rmsd)
        logger.log('the RMSD between Elazar and ResSolv for %s is %.2f' % (aa, rmsd))

    # as long as one residue's RMSD is higher than 0.5, keep iterating
    iteration_num = 1
    while any([elazar_ips[aa].rmsd_ips(MPResSolv_ips[aa]) > 0.5 for aa in AAs]):
        a_profiles = {aa: subtract_IP_from_IP(elazar_ips[aa], MPResSolv_current_ips[aa]) for aa in AAs}
        new_profiles = {aa: subtract_IP_from_IP(MPResSolv_current_ips[aa], a_profiles[aa]) for aa in AAs}
        logger.log('created new profiles by subtracting the last MPResSolv from Elazar:')
        for aa in AAs:
            logger.log('%s %s' % (aa, new_profiles[aa].format_polyval()))
        create_polyval_table(new_profiles, 'ELazaridis_%i.txt' % iteration_num)

        logger.log('rerunning FilterScan analysis...')
        all_cnfgs_df = filterscan_analysis()
        MPResSolv_current_ips = create_insertion_profiles(all_cnfgs_df, 'ResSolv_normed')

        for aa in AAs:
            rmsd = elazar_ips[aa].rmsd_ips(MPResSolv_current_ips[aa])
            rmsd_iter[aa].append(rmsd)
            logger.log('the RMSD between Elazar and ResSolv for %s is %.2f' % (aa, rmsd))

        logger.log('finished iteration number %i' % iteration_num)

        iteration_num += 1
        if iteration_num == 3:
            break

    logger.log('throught the run, the RMSD for the AAs have changes as such:')
    for aa in AAs:
        logger.log('%s %r' % (aa, rmsd_iter[aa]))

    draw_rmsd_plots()
    # if args['show_fig']:
    #     draw_filterscan_profiles(OrderedDict({'full': full_ips, 'no_Menv': noMenv_ips, 'ResSolv': MPResSolv_ips,
    #                                           'diff_ips': diff_ips, 'elazar': elazar_ips}))


def create_polyval_table(diff_ips: dict, file_name: str) -> None:
    """
    creates a file_name with polynom values for all residues:
    A float float float float (highest power first)
    """
    with open(pwd+file_name, 'w+') as fout:
        for aa in AAs:
            fout.write('%s %s\n' % (aa, diff_ips[aa].format_polyval()))
    logger.log('created table at %s' % pwd+file_name)
    if args['change_rosetta_table']:
        shutil.copy(pwd+file_name, rosetta_membranepotential_path + args['polynom_table_file_name'])
        logger.log('copied table to %s' % rosetta_membranepotential_path + args['polynom_table_file_name'])


def draw_filterscan_profiles(ips_dict: OrderedDict) -> None:
    """
    draws all profiles
    """
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.15, hspace=0.45)
    for i, aa in enumerate(AAs):
        plt.subplot(5, 4, 1+i)
        for name, ips in ips_dict.items():
            plt.plot(Z, [x.energy for x in ips[aa].poz_energy], color=color_map[name], label=name)
        # plt.plot(Z, [x.energy for x in full_ips[aa].poz_energy], color='blue', label='full')
        # plt.plot(Z, [x.energy for x in noMenv_ips[aa].poz_energy], color='red', linestyle='dashed', label='no M env')
        # plt.plot(Z, [x.energy-y.energy for x, y in zip(full_ips[aa].poz_energy, noMenv_ips[aa].poz_energy)], color='green')
        # plt.plot(Z, [x.energy for x in elazar_ips[aa].poz_energy], color='black', linestyle='dashed', label='Elazar')
        # plt.plot(Z, [x.energy for x in diff_ips[aa].poz_energy], color='purple', linestyle='dashed', label='Elazar-no_Menv')
        # plt.plot(Z, [np.polyval(diff_ips[aa].polynom, z) for z in Z], color='green', label='diff fit')
        # plt.plot(Z, [np.polyval(MPResSolv_ips[aa].polynom, z) for z in Z], color='orange', label='diff fit')
        # plt.plot(Z, [np.polyval(MPResSolv_remade_ips[aa].polynom, z) for z in Z], color='cyan', label='remade ResSolv')
        plt.title(aa.upper())
        if aa == 'P':
            plt.ylim([-12, 5])
        else:
            plt.ylim([-5, 5])
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    file_location = '%sprofile_comparison.png' % pwd
    plt.savefig(file_location, dpi=600)
    logger.log('saving profile comparison figure to %s' % file_location)
    plt.show()


def draw_sigmoidal_profiles(elazar_ips: dict) -> None:
    from mpl_toolkits.mplot3d import axes3d
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x_neighbors = np.array(range(0, 41))
    z_Z = Z
    y_energy = [[sigmoid(n, 0.5, -4) * pze.energy for pze in elazar_ips['W'].poz_energy] for n in x_neighbors]
    print(x_neighbors)
    print(z_Z)
    print(y_energy)

    ax.plot_wireframe(z_Z, x_neighbors, y_energy, rstride=10, cstride=10)


def sigmoid(x, slope: float, offset: float) -> float:
    return 1. / (1 + np.exp(slope*x+offset))


def create_elazar_ips() -> dict:
    """
    creates {AA: ip} for the Elazar scale
    """
    logger.log('creating InsertionProfiles for Elazar')
    elazar_polyval = MakeHydrophobicityGrade()
    result = {}
    for aa in AAs:
        pos_score = {i+1: a for i, a in enumerate([np.polyval(elazar_polyval[aas_1_3[aa]], z) for z in Z])}
        # adjust kcal/mol to REUs according to kcal/mol=0.57REU.
        # suggested by "Role of conformational sampling in computing mutation-induced changes in protein structure
        # and stability."
        pos_score = {k: v/0.57 for k, v in pos_score.items()}
        ip = InsertionProfile(aa, pos_score=pos_score)
        result[aa] = ip
    logger.log('adjusting kcal/mol to REU by REU=kcal/mol / 0.57')
    return result


def create_insertion_profiles(df: pd.DataFrame, column: str) -> dict:
    """
    creates {AA: InsertionProfile} using energy column in df
    """
    logger.log('creating InsertionProfiles for %s' % column)
    result = {}
    for aa in AAs:
        pos_score = {i: df[((df['aa'] == aa) & (df['pos'] == i))][column].values[0] for i in range(1, num_As+1)}
        ip = InsertionProfile(aa, pos_score=pos_score)
        result[aa] = ip
    return result


def filterscan_analysis() -> pd.DataFrame:
    """
    run the FilteScan protocol on the ployA and return results in data frame
    """
    # run FilterScan protocol to create sclog files for both full and no M env score functions
    if args['full']: # -unmute core.scoring.membrane.MPResSolvEnergy
        command = '%s%s -database %s -parser:protocol %s -s %s -mp:setup:spanfiles %s -nstruct %i -overwrite ' \
                  '-mp:scoring:hbond -mute all' % (rosetta_executables_path, rosetta_script_exec_path,
                                                   rosetta_database_path, protocols_path + MPFilterScan,
                                                   pwd + name + '.pdb', pwd + name + '.span', nstruct)
        logger.log('running FilterScan protocol, issuing command:\n%s' % command)
        os.system(command)

        logger.log('ran FilterScan. finished at %s' % time.strftime("%H:%M_%d%b"))

    # parse both sclog files to {pos: {AA: score}}
    full_fs_log = parse_filterscan_log('memb_hires_full.sclog')
    noMenv_fs_log = parse_filterscan_log('memb_hires_no_Menv.sclog')
    MPResSolv_fs_log = parse_filterscan_log('memb_hires_ResSolv.sclog')

    # create DataFrame
    df = pd.DataFrame()
    for aa in AAs:
        for pos in range(7, 33+1):
            df = df.append({'pos': pos-6, 'aa': aa,
                            'full': full_fs_log[pos][aa],
                            'noMenv': noMenv_fs_log[pos][aa],
                            'res_solv': MPResSolv_fs_log[pos][aa]},
                           ignore_index=True)

    # calculate polyA scores in both score funcitons
    full_A_mean = df[df['aa'] == 'A']['full'].mean()
    noMenv_A_mean = df[df['aa'] == 'A']['noMenv'].mean()
    MPResSolv_A_mean = df[df['aa'] == 'A']['res_solv'].mean()

    # calcualte delta of every mutant and the polyA for both score functions
    df['full_normed'] = df['full']-full_A_mean
    df['noMenv_normed'] = df['noMenv'] - noMenv_A_mean
    df['ResSolv_normed'] = df['res_solv'] - MPResSolv_A_mean

    logger.log('mean of A at full score is %f' % full_A_mean)
    logger.log('mean of A at no M env is %f' % noMenv_A_mean)
    logger.log('mean of A at MPResSolv env is %f' % MPResSolv_A_mean)
    return df


def draw_rmsd_plots() -> None:
    """
    draw the rmsd over iteration plots
    """
    i = 0
    for aa, rmsd_list in rmsd_iter.items():
        plt.subplot(5, 4, 1 + i)
        plt.plot(range(len(rmsd_list)), rmsd_list)
        plt.title(aa)

        i += 1
    plt.show()

def parse_filterscan_log(file_name) -> dict:
    """
    parse a FilterScan run log file, returns dict {position: {AA: score}}
    """
    result = {i: {aa: None for aa in AAs} for i in range(1, num_As+flank_size*2+1)}
    for l in open(file_name, 'r'):
        s = l.split()
        if len(s) >= 4:
            result[int(s[0])][s[2]] = float(s[3])
    return result


def analyse_all_configuration():
    """
    goes over all configuration score files, returns DF of scores
    """
    sc_files = [a for a in os.listdir(pwd+'all_configs') if '.sc' in a]
    for sc_file in sc_files:
        sc_df = score_dict2df(score2dict(pwd+'all_configs/'+sc_file))
        print(sc_df)


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
    profiles = {aa: {pos: None for pos in range(flank_size+1, flank_size+num_As+1)} for aa in aas_names}
    sc_files = [a for a in os.listdir(all_configurations_path) if '.sc' in a]
    for sc_file in sc_files:
        sc = score2dict(all_configurations_path+sc_file)
        df = score_dict2df(sc)
        aa = [a for a in aas_names if a in sc_file][0]
        pos = int(sc_file.split('_')[2].split('.')[0])
        profiles[aa][pos] = np.mean(df['score']) - polyA_total_mean
        print('for res %s in pos %i found %.3f mean with %.3f std' % (aa, pos, np.mean(df['score']), np.std(df['score'])))
    membrane_position = np.linspace(-membrane_half_depth, membrane_half_depth, endpoint=True, num=num_As)
    print('AAAA', membrane_position)
    mm_profiles = {aa: OrderedDict((membrane_position[pos-flank_size-1], profiles[aa][pos])
                                   for pos in range(flank_size+1, flank_size+num_As+1)) for aa in aas_names}
    print(mm_profiles)
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
    try:
        hydrophobicity_grade = open(polyval_file_path, 'r')
    except:
        hydrophobicity_grade = open('/Volumes/labs/fleishman/jonathaw/membrane_prediciton/mother_fucker_MET_LUE_VAL_sym_k_neg.txt', 'r')
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
    for i in range(flank_size+1, flank_size+num_As+1):
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


def create_polyA_fasta() -> None:
    """
    creates a fasta file with a num_As polyA in it
    """
    with open(pwd+name+'.fa', 'w+') as fout:
        fout.write('>polyA\n%s\n' % ''.join(['A']*(num_As+2*flank_size)))
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
                    '-mp:setup:transform_into_membrane 1 -mute all' % (rosetta_executables_path, pwd+name+'.fa'))
    shutil.move(pwd+'helix_from_sequence.pdb', pwd+name+'.pdb')
    print('created an idealised helix from %s and put it in %s' % (name, name+'.pdb'))


def create_spanfile() -> None:
    """
    creates a simple spanfile for the polyA pdb. basically all 1-num_As residues are in the membrane
    """
    with open('%s%s.span' % (pwd, name), 'w+') as fout:
        fout.write('Rosetta-generated spanfile from SpanningTopology object\n%i %i\nantiparallel\nn2c\n\t%i	%i\n'
                   % (flank_size+1, flank_size+num_As, flank_size+1, flank_size+num_As))


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
    global polyA_total_mean, logger, args

    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', default='', type=str)
    parser.add_argument('-queue', default='new-all.q', type=str)
    parser.add_argument('-full', default=False, type=bool)
    parser.add_argument('-show_fig', default=False, type=bool)
    parser.add_argument('-polynom_table_file_name', default='ELazaridis_polynom_table.txt', type=str)
    parser.add_argument('-change_rosetta_table', default=False, type=bool)

    logger = Logger('elazaridis_%s.log' % time.strftime("%H_%M_%d_%b"))

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
        elazar_ips = create_elazar_ips()
        draw_sigmoidal_profiles(elazar_ips)


    else:
        print('unknown mode')

    logger.close()