#!/usr/bin/env python3.5
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available
# (c) under license.
# (c) The Rosetta software is developed by the contributing members of the
# (c) Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington UW
# (c) TechTransfer, email: license@u.washington.edu.
"""Brief:   this script creates a polyval table to read by MPResSolvEnergy. it optimises the fit of the Rosetta energy
            function for membranes to the dsTBL experimental hydrophobicity scale

Params:  use --help to get all params

Example: CreateEMPBenchmark.py -mode main -change_rosetta_table True -full Tru

Remarks: Blah blah blah, blah, blah.

Author:  Jonathan Weinstein

"""
# general imports
import argparse
import os
import sys
import time
import shutil
import subprocess

# 3rd party imports
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
from collections import OrderedDict
from Equation import Expression

# specific imports
from RosettaFilter import score2dict, score_dict2df, df2boxplots
from Logger import Logger

mpl.rc_context(fname="/home/labs/fleishman/jonathaw/.matplotlib/publishable_matplotlibrc")
LSF_USERNAME = 'jonatha'
PWD = os.getcwd() + '/'
all_configurations_path = PWD + 'all_configs/'

ROSETTA_EXECUTABLES_PATH = '/home/labs/fleishman/jonathaw/Rosetta/main/source/bin/'
ROSETTA_SCRIPTS_EXEC_PATH = 'rosetta_scripts.default.linuxgccrelease'
ROSETTA_DATABASE_PATH = '/home/labs/fleishman/jonathaw/Rosetta/main/database/'
ROSETTA_MEM_POTENTIAL_PATH = ROSETTA_DATABASE_PATH + 'scoring/score_functions/MembranePotential/'
PROTOCOLS_PATH = '/home/labs/fleishman/jonathaw/elazaridis/protocols/'
ELAZAR_POLYVAL_PATH = '/home/labs/fleishman/jonathaw/membrane_prediciton/mother_fucker_MET_LUE_VAL_sym_k_neg.txt'
MPMutateRelax_XML = 'MPMutateRelax.xml'
MPFilterScan_XML = 'MPFilterScan_ELazaridis.xml'
MPFilterScanDifferentSFAAs = 'MPFilterScanDifferentSFAAs.xml'

RMSD_THRESHOLD = 0.5
NUM_AAS = 26
POLY_A_NAME = 'polyA'
NSTRUCT = 1  # 0
MEMBRANE_HALF_WIDTH = 15
LAZARIDIS_POLY_DEG = 4
FLANK_SIZE = 6 # 30  # 6
Z = np.linspace(-MEMBRANE_HALF_WIDTH, MEMBRANE_HALF_WIDTH, num=NUM_AAS)
POS_Z_DICT = {i + 1: z for i, z in enumerate(Z)}

aas_names = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG',
             'SER', 'THR', 'VAL', 'TRP', 'TYR']
aas_3_1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K',
           'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
           'TRP': 'W', 'TYR': 'Y'}
aas_1_3 = {v: k for k, v in aas_3_1.items()}
AAs = list('ACDEFGHIKLMNPQRSTVWY')

COLOR_MAP = {'full': 'blue', 'no_Menv': 'grey', 'elazar': 'red', 'ResSolv': 'purple', 'ResSolvCEN': 'pink',
             'fullCEN': 'yellow', 'no_MenvCEN': 'black', 'diff_ips': 'orange', 'fa_intra_rep': 'green',
             'fa_mpsolv': 'pink', 'fa_rep': 'black', 'p_aa_pp': 'blue', 'rama': 'brown'}
RMSD_ITERATIONS = {aa: [] for aa in AAs}


class PoZEnergy:
    """
    container for sequence position, membrane depth and energy. simplifies code.
    """
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
    """
    a single residue insertion profile described by either a sorted list of PoZEnergy instances (poz_energy) or by
    polynum
    """
    def __init__(self, aa: str, pos_score: dict, membrane_half_depth: int=MEMBRANE_HALF_WIDTH,
                 residue_num: int=NUM_AAS):
        """

        """
        self.AA = aa
        self.membrane_half_depth = membrane_half_depth
        self.residue_num = residue_num
        self.pos_score = pos_score
        self.poz_energy = pos_energy_dict_to_PoZEnergy_list(pos_score)
        self.polynom = np.polyfit(Z, [x.energy for x in self.poz_energy], LAZARIDIS_POLY_DEG)

    def __repr__(self):
        res = '<InsertionProfile for %s>\n' % self.AA
        res += '\t<poz_energy>\n'
        for a in self.poz_energy:
            res += '\t\t<%s/>\n' % a
        res += '\t<poz_energy/>\n'
        res += '<polynom: %.2f*z^4 + %.2f*z^3 %.2f*z^2 %.2f*z + %.2f />' % self.polynom
        res += '<InsertionProfile for %s/>\n' % self.AA
        return res

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
    for pos in range(1, NUM_AAS+1):
        result.append(PoZEnergy(pos, POS_Z_DICT[pos], pos_energy_dict[pos]))
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

    # first FilterScan run. using null ResSolv
    full_ips = filterscan_analysis_energy_func('full', residues_to_test=AAs)
    fullCEN_ips = filterscan_analysis_energy_func('fullCEN', residues_to_test=AAs)
    noMenv_ips = filterscan_analysis_energy_func('noMenv', residues_to_test=AAs)
    noMenvCEN_ips = filterscan_analysis_energy_func('noMenvCEN', residues_to_test=AAs)

    # calc the difference InsertionProfiles between Elazar and Rosetta. assign them as the polynom table
    diff_ips = {k: subtract_IP_from_IP(elazar_ips[k], noMenv_ips[k]) for k in AAs}
    diff_ips_CEN = {k: subtract_IP_from_IP(elazar_ips[k], noMenvCEN_ips[k]) for k in AAs}
    create_polyval_table(diff_ips, 'ELazaridis_0.txt')
    create_polyval_table(diff_ips_CEN, 'ELazaridis_CEN_0.txt')

    # analyse Rosetta again.
    MPResSolv_current_ips = filterscan_analysis_energy_func('ResSolv', residues_to_test=AAs)
    MPResSolvCEN_current_ips = filterscan_analysis_energy_func('ResSolvCEN', residues_to_test=AAs)

    draw_filterscan_profiles(OrderedDict({'full': full_ips, 'no_Menv': noMenv_ips, 'ResSolv': MPResSolv_current_ips,
                                          'diff_ips': diff_ips, 'elazar': elazar_ips,
                                          'fullCEN': fullCEN_ips, 'no_MenvCEN': noMenvCEN_ips,
                                          'ResSolvCEN': MPResSolvCEN_current_ips}))

    for aa in AAs:
        rmsd = elazar_ips[aa].rmsd_ips(MPResSolv_current_ips[aa])
        RMSD_ITERATIONS[aa].append(rmsd)
        logger.log('the RMSD between Elazar and ResSolv for %s is %.2f' % (aa, rmsd))

    # as long as one residue's RMSD is higher than threshold, keep iterating
    if args['improve']:
        iteration_num = 1
        while any([elazar_ips[aa].rmsd_ips(MPResSolv_current_ips[aa]) > RMSD_THRESHOLD for aa in AAs]):
            # polynom arithmetics to calibrate the polynoms, only for AAs with RMSD higher than threshold
            a_profiles = {aa: subtract_IP_from_IP(elazar_ips[aa], MPResSolv_current_ips[aa]) if elazar_ips[aa].rmsd_ips(
                MPResSolv_current_ips[aa]) > RMSD_THRESHOLD else MPResSolv_current_ips[aa] for aa in AAs}
            new_profiles = {aa: subtract_IP_from_IP(MPResSolv_current_ips[aa], a_profiles[aa]) if elazar_ips[aa].rmsd_ips(
                MPResSolv_current_ips[aa]) > RMSD_THRESHOLD else MPResSolv_current_ips[aa] for aa in AAs}

            # report the status
            logger.log('created new profiles by subtracting the last MPResSolv from Elazar:')
            for aa in AAs:
                logger.log('%s %s' % (aa, new_profiles[aa].format_polyval()))

            create_polyval_table(new_profiles, 'ELazaridis_%i.txt' % iteration_num)

            logger.log('rerunning FilterScan analysis...')
            MPResSolv_current_ips = filterscan_analysis_energy_func('ResSolv')

            for aa in AAs:
                rmsd = elazar_ips[aa].rmsd_ips(MPResSolv_current_ips[aa])
                RMSD_ITERATIONS[aa].append(rmsd)
                logger.log('the RMSD between Elazar and ResSolv for %s is %.2f' % (aa, rmsd))

            logger.log('finished iteration number %i' % iteration_num)

            iteration_num += 1
            if iteration_num == 3:
                break

        logger.log('throught the run, the RMSD for the AAs have changes as such:')
        for aa in AAs:
            logger.log('%s %r' % (aa, RMSD_ITERATIONS[aa]))

        draw_rmsd_plots()
    # draw_filterscan_profiles(OrderedDict({'full': full_ips, 'no_Menv': noMenv_ips, 'ResSolv': MPResSolv_current_ips,
    #                                       'diff_ips': diff_ips, 'elazar': elazar_ips}))
    if args['show_fig']:
        plt.show()


def draw_rosetta_profiles(args):
    ### create files for running benchmark
    if args['full']:
        create_polyA_fasta()
        sequence_to_idealized_helix()
        create_spanfile()
        trunctate_2nd_mem_res()
    ###
    elazar_ips = create_elazar_ips()
    full_ips = filterscan_analysis_energy_func('full', residues_to_test=AAs)

    MPResSolv_current_ips = filterscan_analysis_energy_func('ResSolv', residues_to_test=AAs, to_dump_pdbs=True)
    MPResSolvCEN_current_ips = filterscan_analysis_energy_func('ResSolv_cen', residues_to_test=AAs, to_dump_pdbs=True)


    # e_term_ips = create_e_term_specific_profiles('./', args['e_terms'])

    dct = OrderedDict({'ResSolv': MPResSolv_current_ips, 'ResSolvCEN': MPResSolvCEN_current_ips,
                       'elazar': elazar_ips, 'full': full_ips})

    # for e_term in args['e_terms']:
    #     dct[e_term] = e_term_ips[e_term]

    draw_filterscan_profiles(dct)
    plt.show()


def create_e_term_specific_profiles(path_pdbs: str, e_terms: list) -> dict:
    """
    creates residue specific insertion profiles for every e_term in the list
    """
    # get energy terms for all res / position combinations
    all_results = get_all_e_terms(path_pdbs)

    # create DataFrame
    df = pd.DataFrame()
    for aa in AAs:
        for pos in range(FLANK_SIZE + 1, NUM_AAS + FLANK_SIZE + 1 + 1):
            dct = {'pos': pos - FLANK_SIZE, 'aa': aa}
            for e_term in e_terms:
                # print('a', all_results[aa][pos].keys())
                dct[e_term] = all_results[aa][pos][e_term]
            df = df.append(dct, ignore_index=True)
    # each res / position combination has a row with values of all e_terms

    e_term_ips = {}
    for e_term in e_terms:
        # find Ala mean
        A_mean = df[df['aa'] == 'A'][e_term].mean()

        # normalise by Ala mean
        df['%s_normed' % e_term] = df[e_term] - A_mean

        # create insertion profiles for e_term
        e_term_ips[e_term] = create_insertion_profiles(df, '%s_normed' % e_term)

    return e_term_ips


def get_all_e_terms(path_pdbs: str) -> dict:
    """
    goes over all residue/position combiantions, and gets the total energy terms out of their files
    """
    result = {a: {} for a in AAs}
    for res in sorted(aas_3_1.keys()):
        for ind in range(1, 2*FLANK_SIZE+NUM_AAS):
            for l in open('polyA.pdbALA%i%s.pdb' % (ind, res), 'r'):
                if 'TOTAL_WTD' in l:
                    s = l.split()
                    result[aas_3_1[res]][ind] = {s[i].replace(':', ''): float(s[i+1]) for i in range(1, len(s), 2)}
    return result


def create_polyval_table(diff_ips: dict, file_name: str) -> None:
    """
    creates a file_name with polynom values for all residues:
    A float float float float (highest power first)
    """
    with open(PWD+file_name, 'w+') as fout:
        for aa in AAs:
            fout.write('%s %s\n' % (aa, diff_ips[aa].format_polyval()))
    logger.log('created table at %s' % PWD + file_name)
    if args['change_rosetta_table']:
        shutil.copy(PWD + file_name, ROSETTA_MEM_POTENTIAL_PATH + args['polynom_table_file_name'])
        logger.log('copied table to %s' % ROSETTA_MEM_POTENTIAL_PATH + args['polynom_table_file_name'])


def draw_filterscan_profiles(ips_dict: OrderedDict) -> None:
    """
    draws all profiles
    """
    plt.figure()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.15, hspace=0.45)
    for i, aa in enumerate(AAs):
        plt.subplot(5, 4, 1+i)
        for name, ips in ips_dict.items():
            plt.plot(Z, [x.energy for x in ips[aa].poz_energy], color=COLOR_MAP[name], label=name)
        plt.title(aa.upper())
        if aa == 'P':
            plt.ylim([-12, 5])
        else:
            plt.ylim([-5, 5])
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    file_location = '%sprofile_comparison.png' % PWD
    plt.savefig(file_location, dpi=600)
    logger.log('saving profile comparison figure to %s' % file_location)


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
        pos_score = {i: df[((df['aa'] == aa) & (df['pos'] == i))][column].values[0] for i in range(1, NUM_AAS + 1)}
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
                  '-mp:scoring:hbond -mute all' % (ROSETTA_EXECUTABLES_PATH, ROSETTA_SCRIPTS_EXEC_PATH,
                                                   ROSETTA_DATABASE_PATH, PROTOCOLS_PATH + MPFilterScan_XML,
                                                   PWD + POLY_A_NAME + '.pdb', PWD + POLY_A_NAME + '.span', NSTRUCT)
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


def filterscan_analysis_energy_func(energy_function: str, residues_to_test: list=AAs, to_dump_pdbs: bool=False) -> dict:
    """
    run the FilteScan protocol on the ployA and return InsertionProfiles dict for energy_function
    """
    # run FilterScan protocol to create sclog files for both full and no M env score functions
    if args['full']: # -unmute core.scoring.membrane.MPResSolvEnergy
        command = '%s%s -database %s -parser:protocol %s -s %s -mp:setup:spanfiles %s -nstruct %i -overwrite ' \
                  '-mp:scoring:hbond -mute all -parser:script_vars energy_function=%s residues_to_test=%s to_dump=%i' % \
                  (ROSETTA_EXECUTABLES_PATH, ROSETTA_SCRIPTS_EXEC_PATH, ROSETTA_DATABASE_PATH,
                   PROTOCOLS_PATH + MPFilterScanDifferentSFAAs, PWD + POLY_A_NAME + '.pdb', PWD + POLY_A_NAME + '.span', NSTRUCT,
                   energy_function, ''.join(residues_to_test), 1 if to_dump_pdbs else 0)
        logger.log('running FilterScan for the energy function %s protocol, issuing command:\n%s' %
                   (energy_function, command))
        os.system(command)

        logger.log('ran FilterScan for energy function %s. finished at %s' %
                   (energy_function, time.strftime("%H:%M_%d%b")))

    # parse both sclog files to {pos: {AA: score}}
    if args['full']:
        shutil.move('temp.sclog', '%s.sclog' % energy_function)
    temp_fs_log = parse_filterscan_log('%s.sclog' % energy_function)
    logger.log('saved the FilterScan sclog to %s.sclog' % energy_function)

    # create DataFrame
    df = pd.DataFrame()
    for aa in AAs:
        for pos in range(FLANK_SIZE+1, NUM_AAS+FLANK_SIZE+1+1):
            df = df.append({'pos': pos-FLANK_SIZE, 'aa': aa, energy_function: temp_fs_log[pos][aa]}, ignore_index=True)

    # calculate polyA scores in both score funcitons
    temp_A_mean = df[df['aa'] == 'A'][energy_function].mean()

    # calcualte delta of every mutant and the polyA for both score functions
    df['%s_normed' % energy_function] = df[energy_function]-temp_A_mean

    logger.log('mean of A for %s is %f' % (energy_function, temp_A_mean))

    ips = create_insertion_profiles(df, '%s_normed' % energy_function)
    return ips


def draw_rmsd_plots() -> None:
    """
    draw the rmsd over iteration plots
    """
    i = 0
    plt.figure()
    for aa, rmsd_list in RMSD_ITERATIONS.items():
        plt.subplot(5, 4, 1 + i)
        plt.plot(range(len(rmsd_list)), rmsd_list)
        plt.title(aa)

        i += 1
    plt.savefig('rmsd_plt.png')


def parse_filterscan_log(file_name) -> dict:
    """
    parse a FilterScan run log file, returns dict {position: {AA: score}}
    """
    result = {i: {aa: None for aa in AAs} for i in range(1, NUM_AAS + FLANK_SIZE * 2 + 1)}
    for l in open(file_name, 'r'):
        s = l.split()
        if len(s) >= 4:
            result[int(s[0])][s[2]] = float(s[3])
    return result


def analyse_all_configuration():
    """
    goes over all configuration score files, returns DF of scores
    """
    sc_files = [a for a in os.listdir(PWD + 'all_configs') if '.sc' in a]
    for sc_file in sc_files:
        sc_df = score_dict2df(score2dict(PWD + 'all_configs/' + sc_file))
        print(sc_df)


def run_all_configurations():
    """
    :return:
    """
    os.mkdir(PWD + 'all_configs')
    os.chdir(PWD + 'all_configs')
    make_all_configurations_jobs()
    os.system('sh command')
    sleep_until_jobs_finish()
    create_profiles()
    os.chdir(PWD)


def create_profiles():
    profiles = {aa: {pos: None for pos in range(FLANK_SIZE + 1, FLANK_SIZE + NUM_AAS + 1)} for aa in aas_names}
    sc_files = [a for a in os.listdir(all_configurations_path) if '.sc' in a]
    for sc_file in sc_files:
        sc = score2dict(all_configurations_path+sc_file)
        df = score_dict2df(sc)
        aa = [a for a in aas_names if a in sc_file][0]
        pos = int(sc_file.split('_')[2].split('.')[0])
        profiles[aa][pos] = np.mean(df['score']) - polyA_total_mean
        print('for res %s in pos %i found %.3f mean with %.3f std' % (aa, pos, np.mean(df['score']), np.std(df['score'])))
    membrane_position = np.linspace(-MEMBRANE_HALF_WIDTH, MEMBRANE_HALF_WIDTH, endpoint=True, num=NUM_AAS)
    print('AAAA', membrane_position)
    mm_profiles = {aa: OrderedDict((membrane_position[pos - FLANK_SIZE - 1], profiles[aa][pos])
                                   for pos in range(FLANK_SIZE + 1, FLANK_SIZE + NUM_AAS + 1)) for aa in aas_names}
    print(mm_profiles)
    mm_polys = {aa: np.polyfit(list(mm_profiles[aa].keys()), list(mm_profiles[aa].values()), LAZARIDIS_POLY_DEG)
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
        hydrophobicity_grade = open(ELAZAR_POLYVAL_PATH, 'r')
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
    z = np.linspace(-MEMBRANE_HALF_WIDTH, MEMBRANE_HALF_WIDTH, endpoint=True, num=NUM_AAS + 1)
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
        plt.xlim([-(MEMBRANE_HALF_WIDTH + 0.5), MEMBRANE_HALF_WIDTH + 0.5])
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
    shutil.copy(PWD + POLY_A_NAME + '.span', './')
    shutil.copy(PWD + POLY_A_NAME + '.pdb', './')
    for i in range(FLANK_SIZE+1, FLANK_SIZE+NUM_AAS+1):
        for aa in aas_names:
            job_args = {
                    '-parser:protocol': '%s%s' % (PROTOCOLS_PATH, MPMutateRelax_XML),
                    '-database': '/home/labs/fleishman/jonathaw/Rosetta/main/database/',
                    '-overwrite': None,
                    '-s': '%s.pdb' % POLY_A_NAME,
                    '-out:suffix': '_%s_%i.pdb' % (aa, i),
                    '-mp:setup:spanfiles': POLY_A_NAME + '.span',
                    '-nstruct': NSTRUCT,
                    '-mp:scoring:hbond': None,
                }
            script_vars = {'mut_pos=': i, 'mut_iden=': aa}
            make_jobs(aa, i, job_args, script_vars)
    print('created all jobs to run for all mutations in all %i positions' % NUM_AAS)


def run_polyA_score_benchmark() -> pd.DataFrame:
    os.system('%s%s -database %s -parser:protocol %s -s %s -mp:setup:spanfiles %s -nstruct %i -overwrite '
              '-mp:scoring:hbond '
              '-parser:script_vars mut_pos=1 mut_iden=ALA' %
              (ROSETTA_EXECUTABLES_PATH, ROSETTA_SCRIPTS_EXEC_PATH, ROSETTA_DATABASE_PATH,
               PROTOCOLS_PATH + MPMutateRelax_XML, PWD + POLY_A_NAME + '.pdb', PWD + POLY_A_NAME + '.span', NSTRUCT))
    shutil.move('%sscore.sc' % PWD, '%spolyA_benchmark.score' % PWD)
    print('ran a FastRelax benchmark and placed the score in polyA_benchmark.score')


def create_polyA_fasta() -> None:
    """
    creates a fasta file with a num_As polyA in it
    """
    with open(PWD+POLY_A_NAME+ '.fa', 'w+') as fout:
        fout.write('>polyA\n%s\n' % ''.join(['A'] * (NUM_AAS + 2 * FLANK_SIZE)))
        logger.log('created polyA file at %s.fa with %i As' % (PWD + POLY_A_NAME, NUM_AAS))


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
                    '-mp:setup:transform_into_membrane 1 -mute all' % (ROSETTA_EXECUTABLES_PATH, PWD + POLY_A_NAME + '.fa'))
    shutil.move(PWD + 'helix_from_sequence.pdb', PWD + POLY_A_NAME + '.pdb')
    logger.log('created an idealised helix from %s and put it in %s' % (POLY_A_NAME, POLY_A_NAME + '.pdb'))


def create_spanfile() -> None:
    """
    creates a simple spanfile for the polyA pdb. basically all 1-num_As residues are in the membrane
    """
    with open('%s%s.span' % (PWD, POLY_A_NAME), 'w+') as fout:
        fout.write('Rosetta-generated spanfile from SpanningTopology object\n%i %i\nantiparallel\nn2c\n\t%i	%i\n'
                   % (FLANK_SIZE + 1, FLANK_SIZE + NUM_AAS, FLANK_SIZE + 1, FLANK_SIZE + NUM_AAS))


def trunctate_2nd_mem_res() -> None:
    """
    truncate the 2nd membrane residue from the helix pdb. for some reason Rosetta doesn't like it
    """
    with open('%s%s.pdb' % (PWD, POLY_A_NAME), 'r') as fin:
        with open('%s%s.tmp' % (PWD, POLY_A_NAME), 'w+') as fout:
            for l in fin.read().split('\n'):
                if 'XXXX' not in l:
                    fout.write(l+'\n')
    shutil.move('%s%s.tmp' % (PWD, POLY_A_NAME), '%s%s.pdb' % (PWD, POLY_A_NAME))
    logger.log('removed the second membrane residue')

if __name__ == '__main__':
    global polyA_total_mean, logger, args

    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', default='', type=str)
    parser.add_argument('-queue', default='new-all.q', type=str)
    parser.add_argument('-full', default=False, type=bool)
    parser.add_argument('-show_fig', default=False, type=bool)
    parser.add_argument('-polynom_table_file_name', default='ELazaridis_polynom_table.txt', type=str)
    parser.add_argument('-change_rosetta_table', default=False, type=bool)
    parser.add_argument('-improve', default=False, type=bool)
    parser.add_argument('-e_terms', default=['fa_intra_rep', 'fa_mpsolv', 'fa_rep', 'p_aa_pp', 'rama'])

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

    elif args['mode'] == 'draw_profiles':
        draw_rosetta_profiles(args)

    elif args['mode'] == 'test':
        create_e_term_specific_profiles('./', ['fa_rep', 'fa_mpsolv'])

    else:
        print('unknown mode')

    logger.close()