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
import copy
import shutil
import subprocess
from scipy.signal import argrelextrema

# 3rd party imports
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from collections import OrderedDict
# from Equation import Expression
from scipy import interpolate


# specific imports
from RosettaFilter import score2dict, score_dict2df, df2boxplots
from Logger import Logger

mpl.rc_context(fname="/home/labs/fleishman/jonathaw/.matplotlib/publishable_matplotlibrc")
LSF_USERNAME = 'jonatha'
# PWD = os.getcwd() + '/'
# all_configurations_path = PWD + 'all_configs/'

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
TOTAL_HALF_WIDTH = 134.5/2 # 50
LAZARIDIS_POLY_DEG = 4
FLANK_SIZE = 30  # 6
TOTAL_AAS = NUM_AAS + (FLANK_SIZE * 2)
Z = np.linspace(-MEMBRANE_HALF_WIDTH, MEMBRANE_HALF_WIDTH, num=NUM_AAS)
Z_total = np.linspace(-TOTAL_HALF_WIDTH, TOTAL_HALF_WIDTH, num=TOTAL_AAS)
POS_Z_DICT_total = {i + 1: z for i, z in enumerate(Z_total)}
POS_Z_DICT = {i + 1: z for i, z in enumerate(Z)}
POS_RANGE = range(1, TOTAL_AAS+1)
aas_names = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG',
             'SER', 'THR', 'VAL', 'TRP', 'TYR']
aas_3_1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K',
           'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
           'TRP': 'W', 'TYR': 'Y'}
aas_1_3 = {v: k for k, v in aas_3_1.items()}
AAs = list('ACDEFGHIKLMNPQRSTVWY')
skip_aas = ['P'] #, 'T', 'S']
COLOR_MAP = {'full': 'blue', 'no_Menv': 'grey', 'ResSolv': 'purple',
             'fullCEN': 'blue', 'no_MenvCEN': 'grey', 'ResSolvCEN': 'purple',
             'elazar': 'red', 'diff_ips': 'orange', 'diff_ips_CEN': 'orange',
             'fa_intra_rep': 'green', 'fa_mpsolv': 'pink', 'fa_rep': 'black',
             'p_aa_pp': 'blue', 'rama': 'brown', 'no_res_solv': 'blue'}
RMSD_ITERATIONS = {aa: [] for aa in AAs}

z_range_aa = {aa: [-20, +20] if aa not in ['R', 'K', 'H'] else [-23, +20] for
              aa in AAs}

SPLINE_SMOOTHNESS = 0
SPLINE_LIM = 25 # 30 1Nov

class InsertionProfile:
    """
    a single residue insertion profile described by either a sorted list of PoZEnergy instances (poz_energy) or by
    polynum
    """
    def __init__(self, aa: str, pos_score: dict, membrane_half_depth: int=MEMBRANE_HALF_WIDTH,
                 residue_num: int=NUM_AAS, adjust_extra_membranal=True, poly_edges: list=[]):
        """

        """
        self.AA = aa
        self.membrane_half_depth = membrane_half_depth
        self.residue_num = residue_num
        # self.extra_membrane_adjusted = adjust_extra_membranal
        # if adjust_extra_membranal:
        #     for pos in pos_score.keys():
        #         if -20 > POS_Z_DICT_total[pos] or +20 < POS_Z_DICT_total[pos]:
        #             pos_score[pos] = 0.0
        self.pos_score = pos_score
        # self.poz_energy = pos_energy_dict_to_PoZEnergy_list(pos_score)
        # self.polynom = np.polyfit(Z_total, [x.energy for x in self.poz_energy], LAZARIDIS_POLY_DEG)
        self.extramembrane_adjusted = False
        self.poly_edges = poly_edges

    def __repr__(self):
        res = '<InsertionProfile for %s>\n' % self.AA
        res += '\t<poz_energy>\n'
        for a in self.poz_energy:
            res += '\t\t<%s/>\n' % a
        res += '\t<poz_energy/>\n'
        print('ddd', self.polynom)
        try:
            res += '<polynom: %.2f*z^4 + %.2f*z^3 %.2f*z^2 %.2f*z + %.2f />' % self.polynom
        except:
            pass
        res += '<InsertionProfile for %s/>\n' % self.AA
        return res

    def polynom_at_z(self, z: float) -> float:
        """
        returns the polynoms value at z
        :type z: float
        """
        return np.polyval(self.polynom, z)

    def format_polyval(self):
        """
        print the polyval values for table
        """
        return ' '.join([str(a) for a in self.polynom])

    def format_spline_energies(self):
        """
        :return: string of all energies separated by spaces
        """
        # if self.AA not in ['R', 'K', 'H', 'P']:
            # return ' '.join(
                # str(self.pos_score[pos] if -15 <= POS_Z_DICT_total[pos] <= 15 else 0.0) for pos in POS_RANGE)
        # elif self.AA == 'P':
            # return ' '.join("0.0" for pos in POS_RANGE)
        # else:
            # return ' '.join(
                # str(self.pos_score[pos] if -23 <= POS_Z_DICT_total[pos] <= 15 else 0.0) for pos in POS_RANGE)
        if self.AA not in skip_aas:
            # return ' '.join(str(self.pos_score[pos] if z_range_aa[self.AA][0]
                                # <= POS_Z_DICT_total[pos] <=
                                # z_range_aa[self.AA][1] else 0.0) for pos in POS_RANGE)
            return ' '.join(str(self.pos_score[pos] if -SPLINE_LIM <= POS_Z_DICT_total[pos] <= SPLINE_LIM else 0.0) for pos in POS_RANGE)
        else:
            logger.log("creating %s spline as 0.0" % self.AA)
            return ' '.join("0.0" for pos in POS_RANGE)

    def rmsd_ips(self, other) -> float:
        """
        retruns the difference between the IPs calculated as by RMSD over Z
        """
        res = 0.0
        for pos in range(1, TOTAL_AAS+1):
            if self.poly_edges[0] <= POS_Z_DICT_total[pos] <= self.poly_edges[1]:
                res += (self.pos_score[pos] - other.pos_score[pos])**2
        return np.sqrt(np.mean(res))

    def adjust_exta_membrane(self):
        """
        set all positions outside [-20, 20] to 0. use only for setting splines in Rosetta
        :return:
        """
        self.extramembrane_adjusted = True
        print('ADJUSTING !!! !STOP ME !!!!')
        sys.exit()
        for pos in POS_RANGE:
            if -15 > POS_Z_DICT_total[pos] or +15 < POS_Z_DICT_total[pos]:
                self.pos_score[pos] = 0


def pos_energy_dict_to_PoZEnergy_list(pos_energy_dict: dict) -> list():
    """
    creates an ordered list of PoZEnergy instances corresponding to their positions
    """
    result = []
    for pos in range(1, TOTAL_AAS+1):
        result.append(PoZEnergy(pos, POS_Z_DICT_total[pos], pos_energy_dict[pos]))
    return result


def subtract_IP_from_IP(ip1: InsertionProfile, ip2: InsertionProfile, verbose: bool = False, smooth: bool=True) -> InsertionProfile:
    """
    """
    new_pos_score = {}
    if not smooth:
        for pos in POS_RANGE:
            new_pos_score[pos] = ip1.pos_score[pos] - ip2.pos_score[pos]
    else:
        # smooth the transition from water to membrane between +/-15A to +/-25A for resulting splines
        y, x = [], []
        for pos in POS_RANGE:
            if -SPLINE_LIM > POS_Z_DICT_total[pos] or POS_Z_DICT_total[pos] > SPLINE_LIM:
                y.append(0.0)
                x.append(pos)
            elif ip1.poly_edges[0] <= POS_Z_DICT_total[pos] <= ip1.poly_edges[1]:
                y.append(ip1.pos_score[pos] - ip2.pos_score[pos])
                x.append(pos)
        tck = interpolate.splrep(x, y, s=SPLINE_SMOOTHNESS)
        new_pos_score = {pos: interpolate.splev(pos, tck) 
                         if -SPLINE_LIM <= POS_Z_DICT_total[pos] <= +SPLINE_LIM else 0.0 
                         for pos in POS_RANGE}

    return InsertionProfile(ip1.AA, new_pos_score)#, adjust_extra_membranal=ip1.extra_membrane_adjusted)


def add_IP_to_IP(ip1: InsertionProfile, ip2: InsertionProfile, verbose: bool = False) -> InsertionProfile:
    """
    """
    new_pos_score = {}
    for pos in POS_RANGE:
        new_pos_score = ip1.pos_score[pos] = ip2.pos_score[pos]
        if verbose:
            print(pos, ip1.pos_score[pos], ip2.pos_score[pos], ip1.pos_score[pos]+ip2.pos_score[pos])
    return InsertionProfile(ip1.AA, new_pos_score)#, adjust_extra_membranal=ip1.extra_membrane_adjusted)


def calibrate_energy_functions(args):
    """
    
    :param args:
    :return:
    """
    global PWD

    # fa_cen_for_scores = OrderedDict(
        # dict(score0='centroid', score1='centroid', score2='centroid', score3='centroid', score5='centroid',
             # talaris2014='fa_standard'))
    # fa_cen_for_scores = OrderedDict(dict(score5='centroid', talaris2014='fa_standard'))
    fa_cen_for_scores = OrderedDict(dict(score5='centroid', beta_nov15='fa_standard'))
    score_funcs_to_calibrate = list(fa_cen_for_scores.keys())
    original_dir = os.getcwd()
    logger.log("will calibrate the score functions %r" % score_funcs_to_calibrate)
    for en_func in score_funcs_to_calibrate:
        if en_func != 'beta_nov15':
            continue
        # if en_func != 'talaris2014': continue
        os.mkdir('%s/%s' % (original_dir, en_func))
        os.chdir('%s/%s' % (original_dir, en_func))
        logger.create_header("calibrating %s" % en_func)
        PWD = '%s/%s/' % (original_dir, en_func)

        calibrate_function(en_func + '_elazaridis', fa_cen=fa_cen_for_scores[en_func])

        os.chdir('%s/' % original_dir)


def calibrate_function(score_func='talaris2014_elazaridis', fa_cen='fa_standard'):
    # create files for running benchmark
    if args['full']:
        if not args['use_made_pdb']:
            create_polyA_fasta()
            sequence_to_idealized_helix()
        else:
            copy_path = '/home/labs/fleishman/jonathaw/elazaridis/file_safe/polyA_inMemb.pdb'
            logger.log('USING THE PREMADE SAVE PDB !!! from %s' % copy_path)
            shutil.copy(copy_path, 'polyA.pdb')
        # create_spanfile()
        # trunctate_2nd_mem_res()
        # os.system('/home/labs/fleishman/jonathaw/elazaridis/protocols/make_csts.sh %s > %s' % (PWD + POLY_A_NAME + '.pdb', PWD+POLY_A_NAME+'.cst'))

    # elazar_ips = create_elazar_ips()

    # first FilterScan run. using null ResSolv
    full_ips = filterscan_analysis_energy_func(score_func, res_solv_weight=0.0, fa_cen=fa_cen, residues_to_test=AAs,
                                               print_xml=True, adjust_extra_membranal=False)

    logger.create_header('creating and adjusting elazar profiles')
    elazar_ips = create_elazar_ips()
    for aa in ['V', 'M', 'G', 'T', 'S']:
        if aa not in ['S']:
            rhs_avg = np.mean([full_ips[aa].pos_score[pos] for pos in range(6, 17)])
        elif aa == 'S':
            rhs_avg = 0.5 # local maxima in elazar S profile
        logger.log('for res %s found mean to be %.2f' % (aa, rhs_avg))
        for pos in POS_RANGE:
            if aa in ['V', 'M']:
                if elazar_ips[aa].pos_score[pos] > rhs_avg:
                    elazar_ips[aa].pos_score[pos] = rhs_avg
            if aa in ['G', 'T']:
                if elazar_ips[aa].pos_score[pos] < rhs_avg:
                    elazar_ips[aa].pos_score[pos] = rhs_avg
            if aa in ['S']:
                    elazar_ips[aa].pos_score[pos] = rhs_avg

    # calc the difference InsertionProfiles between Elazar and Rosetta. assign them as the polynom table
    diff_ips = {0: {k: subtract_IP_from_IP(elazar_ips[k], full_ips[k]) for k in AAs}}

    create_spline_table(diff_ips[0], 'spline_%s_fa.txt' % score_func, 'spline_test_%s.txt' % ('fa' if fa_cen == 'fa_standard' else 'cen'), args['note'])

    # analyse Rosetta again.
    MPResSolv_current_ips = {
        0: filterscan_analysis_energy_func(score_func, res_solv_weight=1.0, fa_cen=fa_cen, residues_to_test=AAs,
                                           adjust_extra_membranal=False, to_dump_pdbs=False)}
    for aa in AAs:
        rmsd = elazar_ips[aa].rmsd_ips(MPResSolv_current_ips[0][aa])
        RMSD_ITERATIONS[aa].append(rmsd)
        logger.log('the RMSD between Elazar and ResSolv for %s is %.2f' % (aa, rmsd))

    # as long as one residue's RMSD is higher than threshold, keep iterating
    if args['improve']:
        args['full'] = True
        rmsds = {aa: elazar_ips[aa].rmsd_ips(MPResSolv_current_ips[0][aa]) for
                 aa in AAs if aa not in skip_aas}
        fixed_ips = {}
        iter_num = 1
        while any([rmsd > RMSD_THRESHOLD for rmsd in rmsds.values()]) and iter_num < 10:
            aas_improve = [aa for aa in AAs if aa not in skip_aas if rmsds[aa] > RMSD_THRESHOLD]
            if aas_improve == []:
                logger.create_header('finished improving')
                break
            logger.log('starting round %i for AAs %s' % (iter_num, aas_improve))
            diff_ips[iter_num] = {'A': diff_ips[iter_num-1]['A']}

            # check which residues are "good enough" by RMSD, and fix them.
            for aa in AAs:
                if aa in skip_aas:
                    continue
                RMSD_ITERATIONS[aa].append(rmsds[aa])
                if aa not in fixed_ips.keys() and aa not in aas_improve:
                    logger.log('fixing %s on round %i' % (aa, iter_num))
                    fixed_ips[aa] = diff_ips[iter_num-1][aa]
                if aa not in aas_improve:
                    diff_ips[iter_num][aa] = fixed_ips[aa]

            for aa in aas_improve:
                if aa in skip_aas:
                    continue
                logger.log('improve %s at %.2f round %i' % (aa, rmsds[aa], iter_num))

                # create a spline that describes the required profile. in the elazar range (-23, 15 or -15, 15)
                # it will be what is required to get to elazar within the membrane.
                # in (inf, -25) and (+25, inf) it is 0.
                y, x = [], []
                for pos in POS_RANGE:
                    if -SPLINE_LIM > POS_Z_DICT_total[pos] or POS_Z_DICT_total[pos] > SPLINE_LIM:
                        y.append(0.0)
                        x.append(pos)
                    elif elazar_ips[aa].poly_edges[0] <= POS_Z_DICT_total[pos] <= elazar_ips[aa].poly_edges[1]:
                        # train the spline on the difference between the profile from the previous iteration
                        # and what is reuqired to pull it closer to the Elazar 
                        y.append(diff_ips[iter_num-1][aa].pos_score[pos] + elazar_ips[aa].pos_score[pos] 
                                 - MPResSolv_current_ips[iter_num-1][aa].pos_score[pos])
                        x.append(pos)
                tck = interpolate.splrep(x, y, s=SPLINE_SMOOTHNESS)
                diff_ips[iter_num][aa] = InsertionProfile(aa, {pos: interpolate.splev(pos, tck) 
                                                               if -SPLINE_LIM <= POS_Z_DICT_total[pos] <= +SPLINE_LIM else 0.0 
                                                               for pos in POS_RANGE})
            create_spline_table(diff_ips[iter_num], 'spline_%s_fa_%i.txt' % (score_func, iter_num),
                                'spline_test_%s.txt' % ('fa' if fa_cen == 'fa_standard' else 'cen'), args['note'])
            MPResSolv_current_ips[iter_num] = filterscan_analysis_energy_func(score_func, res_solv_weight=1.0,
                                                                              fa_cen=fa_cen,
                                                                              residues_to_test=AAs,
                                                                              adjust_extra_membranal=False,
                                                                              to_dump_pdbs=False)

            rmsds = {aa:
                     elazar_ips[aa].rmsd_ips(MPResSolv_current_ips[iter_num][aa])
                     for aa in AAs if aa != 'P'}
            iter_num += 1

        draw_rmsd_plots()
        draw_filterscan_profiles(OrderedDict({'ResSolv': MPResSolv_current_ips[iter_num-1],
                                              'elazar': elazar_ips}), cen_fa='%s_%s' % (score_func, fa_cen))
        logger.log('finished calibrating %s %s' % (score_func, fa_cen))
        logger.log('got these RMSDs:')
        for aa in AAs:
            if aa not in skip_aas:
                logger.log('for %s got %.2f after %i rounds' % (aa, rmsds[aa], iter_num))


def follow_ip(aa, ips, title):
    points = [43] # [10, 35, 43, 53, 76]
    print('flw %s %s, %s' % (title, aa, ', '.join('%i: %.2f' % (i, ips[aa].pos_score[i]) for i in points)))


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
    fullCEN_ips = filterscan_analysis_energy_func('fullCEN', residues_to_test=AAs)

    MPResSolv_current_ips = filterscan_analysis_energy_func('talaris2014_elazaridis', residues_to_test=AAs, to_dump_pdbs=True)
    MPResSolvCEN_current_ips = filterscan_analysis_energy_func('ResSolvCEN', residues_to_test=AAs, to_dump_pdbs=True)


    # e_term_ips = create_e_term_specific_profiles('./', args['e_terms'])

    dct = OrderedDict({'ResSolv': MPResSolv_current_ips, 'ResSolvCEN': MPResSolvCEN_current_ips,
                       'elazar': elazar_ips, 'full': full_ips})

    # for e_term in args['e_terms']:
    #     dct[e_term] = e_term_ips[e_term]

    draw_filterscan_profiles(dct)
    plt.show()


def draw_rosetta_profiles_fa_cen(args):
    global PWD
    PWD = os.getcwd()+'/'
    ### create files for running benchmark
    if args['full']:
        create_polyA_fasta()
        sequence_to_idealized_helix()
        create_spanfile()
        trunctate_2nd_mem_res()
    ###
    elazar_ips = create_elazar_ips()

    # caluclate and draw full-atom level terms
    # fa_terms = ['fa_atr', 'fa_rep', 'fa_intra_rep', 'fa_pair', 'fa_dun', 'hbond_lr_bb', 'hbond_sr_bb', 'hbond_sc', 'p_aa_pp'] # canceled fa_sol, fa_plane, weight=0, 'fa_mpsolv', 'fa_mpenv','fa_mpenv_smooth' , 'ref, 'hbond_bb_sc'
    # full_ips = filterscan_analysis_energy_func(args['energy_func_fa'], residues_to_test=AAs, to_dump_pdbs=True)
    fa_e_term_ips, fa_terms = create_e_term_specific_profiles(args, './', args['energy_func_fa'])

    print('GOT HESES TERMS', fa_terms)
    for term in fa_terms:
        plt.figure()
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.15, hspace=0.45)
        for i, aa in enumerate(AAs):
            plt.subplot(5, 4, 1 + i)
            # if aa != 'P': continue
            # plt.plot(Z_total, [x.energy for x in fa_e_term_ips[term][aa].poz_energy], color='k', label=term)
            plt.plot(Z_total, [fa_e_term_ips[term][aa].pos_score[pos] for pos in POS_RANGE], color='k', label=term)
            plt.title(aa.upper())
            # if aa == 'P':
                # plt.ylim([-12, 5])
            # else:
                # plt.ylim([-5, 5])
        plt.suptitle(term)
        plt.savefig('fa_%s.png' % term)
        plt.close()


def create_e_term_specific_profiles(args, path_pdbs: str, energy_func, res_solv_weight: float=0.0) -> dict:
    """
    creates residue specific insertion profiles for every e_term in the list
    """
    global PWD
    PWD = os.getcwd()+'/'
    if args['full']:
        create_polyA_fasta()
        sequence_to_idealized_helix()
        create_spanfile()
        trunctate_2nd_mem_res()
        filterscan_analysis_energy_func('talaris2014_elazaridis',
                                        res_solv_weight=0.0,
                                        residues_to_test=AAs,
                                        to_dump_pdbs=True,
                                        fa_cen='fa_standard',
                                       )
    # get energy terms for all res / position combinations
    all_results, e_terms = get_all_e_terms(path_pdbs)

    # create DataFrame
    df_ = pd.DataFrame()
    for aa in AAs:
        for pos in range(1, TOTAL_AAS + 1):
            dct = {'pos': pos, 'aa': aa}
            for e_term in e_terms:
                # print('a', all_results[aa][pos].keys())
                dct[e_term] = all_results[aa][pos][e_term]
            df_ = df_.append(dct, ignore_index=True)
    # each res / position combination has a row with values of all e_terms

    e_term_ips = {}
    for e_term in e_terms:
        # find Ala mean
        A_mean = df_[df_['aa'] == 'A'][e_term].mean()

        # normalise by Ala mean
        df_['%s_normed' % e_term] = df_[e_term] - A_mean

        # create insertion profiles for e_term
        e_term_ips[e_term] = create_insertion_profiles(df_, '%s_normed' % e_term, adjust_extra_membranal=False, smooth=res_solv_weight == 1.0)

    return e_term_ips, e_terms


def get_all_e_terms(path_pdbs: str) -> (dict, list):
    """
    goes over all residue/position combiantions, and gets the total energy terms out of their files
    """
    result = {a: {} for a in AAs}
    all_terms = []
    for res in sorted(aas_3_1.keys()):
        for ind in range(1, TOTAL_AAS+1):
            for l in open('polyA.pdbALA%i%s.pdb' % (ind, res), 'r'):
                if 'TOTAL_WTD' in l:
                    s = l.split()
                    result[aas_3_1[res]][ind] = {s[i].replace(':', ''): float(s[i+1]) for i in range(1, len(s), 2)}
                    if not all_terms:
                        all_terms = [s[i].replace(':', '') for i in range(1, len(s), 2)]
    return result, all_terms


def create_polyval_table(diff_ips: dict, file_name: str, rosetta_table_name) -> None:
    """
    creates a file_name with polynom values for all residues:
    A float float float float (highest power first)
    """
    with open(PWD+file_name, 'w+') as fout:
        for aa in AAs:
            fout.write('%s %s\n' % (aa, diff_ips[aa].format_polyval()))
    logger.log('created table at %s' % PWD + file_name)
    if args['change_rosetta_table']:
        shutil.copy(PWD + file_name, ROSETTA_MEM_POTENTIAL_PATH + rosetta_table_name)
        logger.log('copied table to %s' % ROSETTA_MEM_POTENTIAL_PATH + rosetta_table_name)


def create_spline_table(diff_ips_: dict, file_name: str, rosetta_spline_table_name, note=None) -> None:
    """
    :param diff_ips_: {AA: IP}
    :param file_name: local file in which to create the table
    :param rosetta_spline_table_name: destination to which place the table if "change_rosetta_table"
    :return: None
    """
    with open(PWD+file_name, 'w+') as fout:
        fout.write('# splines generated on %s\n' % time.strftime("%H_%M_%d_%b"))
        if note is not None:
            fout.write('# NOTE: %s\n' % note)
        for aa in AAs:
            if aa in skip_aas:
                fout.write('%s %s\n' % (aa, InsertionProfile(aa, {}).format_spline_energies()))
            else:
                fout.write('%s %s\n' % (aa, diff_ips_[aa].format_spline_energies()))
    logger.log('created table at %s' % PWD + file_name)
    if args['change_rosetta_table']:
        shutil.copy(PWD + file_name, ROSETTA_MEM_POTENTIAL_PATH + rosetta_spline_table_name)
        logger.log('copied table to %s' % ROSETTA_MEM_POTENTIAL_PATH + rosetta_spline_table_name)


def just_draw_current_profiles():
    global PWD
    args['full'] = True
    PWD = os.getcwd()+'/'
    create_polyA_fasta()
    sequence_to_idealized_helix()
    create_spanfile()
    trunctate_2nd_mem_res()
    elazar_ips = create_elazar_ips()
    MPResSolv_current_ips = filterscan_analysis_energy_func('talaris2014_elazaridis', 
                                                            0.0, 
                                                            'fa_standard', 
                                                            residues_to_test=AAs,
                                                            adjust_extra_membranal=False, 
                                                            to_dump_pdbs=False)
    with_ressolv = filterscan_analysis_energy_func('talaris2014_elazaridis', 
                                                   1.0, 
                                                   'fa_standard', 
                                                   residues_to_test=AAs, 
                                                   adjust_extra_membranal=False, 
                                                   to_dump_pdbs=False)
    dct = {'elazar': elazar_ips, 'ResSolv': with_ressolv, 'no_res_solv': MPResSolv_current_ips}
    draw_filterscan_profiles(dct, show=True)


def draw_filterscan_profiles(ips_dict: OrderedDict, cen_fa='fa', show: bool = False) -> None:
    """
    draws all profiles
    """
    plt.figure()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.15, hspace=0.45)
    for i, aa in enumerate(AAs):
        plt.subplot(5, 4, 1+i)
        for name, ips in ips_dict.items():
            plt.plot(Z_total, [ips[aa].pos_score[pos] for pos in POS_RANGE], color=COLOR_MAP[name], label=name)
        plt.title(aa.upper())
        if aa != 'P':
            plt.ylim([-5, 5])
        plt.xlim([-50, 50])
        plt.axvline(z_range_aa[aa][0], color='grey')
        plt.axvline(z_range_aa[aa][1], color='grey')
        plt.axvline(-SPLINE_LIM, color='blue')
        plt.axvline(+SPLINE_LIM, color='blue')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    file_location = '%sprofile_comparison_%s.png' % (PWD, cen_fa)
    plt.savefig(file_location, dpi=600)
    logger.log('saving profile comparison figure to %s' % file_location)
    if show:
        plt.show()
    else:
        plt.close()


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
        pos_score = {i+1: a for i, a in enumerate([np.polyval(elazar_polyval[aas_1_3[aa]], z)
                                                  if -20 <= z <= +20 else 0.0 for z in Z_total])}
        # pos_score = {i+1: a for i, a in enumerate([np.polyval(elazar_polyval[aas_1_3[aa]], z)
                                                   # if -MEMBRANE_HALF_WIDTH <= z <= MEMBRANE_HALF_WIDTH else 0.0
                                                   # for z in Z_total])}
        # extend positive inside rule into the IN side
        if aa in ['R', 'K', 'H']:
            edge_score = np.polyval(elazar_polyval[aas_1_3[aa]], -20)
            logger.log('adjusting %s -23 <= z <= -20 to be %.2f' % (aa, edge_score))
            for pos in POS_RANGE:
                if -23 <= POS_Z_DICT_total[pos] <= -20:
                    pos_score[pos] = edge_score

        # force DEQN profiles to be linear, equivalent to the membrane core energy
        if aa in ['D', 'E', 'Q', 'N']:
            core_avg = np.mean([pos_score[pos] for pos in POS_RANGE if -10 <= POS_Z_DICT_total[pos] <= 10])
            for pos in POS_RANGE:
                if -20 <= POS_Z_DICT_total[pos] <= 20:
                    pos_score[pos] = core_avg

        # force H to be linear in (0, +20)
        if aa == 'H':
            zero_to_twenty_avg = np.max([pos_score[pos] for pos in POS_RANGE if -MEMBRANE_HALF_WIDTH <= POS_Z_DICT_total[pos] <= MEMBRANE_HALF_WIDTH])
            for pos in POS_RANGE:
                if 0 <= POS_Z_DICT_total[pos] <= 20:
                    pos_score[pos] = zero_to_twenty_avg

        # for all other AAs, stop polynom influence at furthest max/min points. this
        # prevents the tiny troffs created by the polynom min points
        if aa not in ['D', 'Q', 'N', 'E', 'H', 'A', 'T', 'G', 'K', 'R']:
            x = np.array([pos_score[pos] for pos in POS_RANGE])
            max_pnts = [a+1 for a in argrelextrema(x, np.greater)[0]]
            min_pnts = [a+1 for a in argrelextrema(x, np.less)[0]]
            # for I, there's a minimum at position 56 which is WRONG, skip it...
            if aa == 'I':
                min_pnts = [min_pnts[0]]
            edge_pnts = [POS_Z_DICT_total[np.min(max_pnts+min_pnts)], POS_Z_DICT_total[np.max(max_pnts+min_pnts)]]
            if not any([a > 10 for a in edge_pnts]):
                edge_pnts[1] = z_range_aa[aa][1]
            if not any([a < -10 for a in edge_pnts]):
                edge_pnts[0] = z_range_aa[aa][0]
        else:
            edge_pnts = [z_range_aa[aa][0], z_range_aa[aa][1]]

        logger.log('for %s setting the poly edges at %.2f, %.2f' % (aa, edge_pnts[0], edge_pnts[1]))

        # adjust kcal/mol to REUs according to kcal/mol=0.57REU.
        # suggested by "Role of conformational sampling in computing mutation-induced changes in protein structure
        # and stability."
        pos_score = {k: v/0.57 for k, v in pos_score.items()}
        ip = InsertionProfile(aa, pos_score=pos_score, poly_edges=edge_pnts)
        result[aa] = ip
    logger.log('adjusting kcal/mol to REU by REU=kcal/mol / 0.57')
    return result


def create_insertion_profiles(df: pd.DataFrame, column: str, adjust_extra_membranal: bool=True, smooth: bool=False) -> dict:
    """
    creates {AA: InsertionProfile} using energy column in df
    """
    logger.log('creating InsertionProfiles for %s' % column)
    result = {}
    for aa in AAs:
        # pos_score = {i: df[((df['aa'] == aa) & (df['pos'] == i))][column].values[0] for i in range(1, NUM_AAS + 1)}
        pos_score = {i: df[((df['aa'] == aa) & (df['pos'] == i))][column].values[0] for i in POS_RANGE} # switched for -50 to 50

        ip = InsertionProfile(aa, pos_score=pos_score, adjust_extra_membranal=adjust_extra_membranal)
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


def filterscan_analysis_energy_func(energy_function: str, res_solv_weight: float, fa_cen: str,
                                    residues_to_test: list=AAs, to_dump_pdbs: bool=False,
                                    adjust_extra_membranal: bool=True, print_xml: bool=False) -> dict:
    """
    energy_function = which energy function to calibrate
    run the FilteScan protocol on the ployA and return InsertionProfiles dict for energy_function
    """
    # run FilterScan protocol to create sclog files for both full and no M env score functions
    if args['full']: # -unmute core.scoring.membrane.MPResSolvEnergy
        if energy_function is None:
            print('no energy function provided!')
            sys.exit() # -corrections::beta_nov15 -score::elec_memb_sig_die
        command = '%s%s -database %s -parser:protocol %s -s %s -nstruct %i -overwrite ' \
                  '-mp:scoring:hbond -mute all -ex1 -ex2 -ex3 -ex4 -parser:script_vars energy_function=%s ' \
                'residues_to_test=%s to_dump=%i res_solv_weight=%.2f fa_or_cen=%s' % \
                  (ROSETTA_EXECUTABLES_PATH, ROSETTA_SCRIPTS_EXEC_PATH, ROSETTA_DATABASE_PATH,
                   PROTOCOLS_PATH + MPFilterScanDifferentSFAAs, PWD + POLY_A_NAME + '.pdb', NSTRUCT,
                   energy_function, ''.join(residues_to_test), 1 if to_dump_pdbs else 0, res_solv_weight, fa_cen)
        # -mp:setup:spanfiles %s , PWD + POLY_A_NAME + '.span'
        if 'beta' in energy_function:
            command += ' -corrections::beta_nov15 '
            logger.log('ADDING CORRECTIONS FOR beta')
            if args['elec_memb_sig_die'] and 'beta' in energy_function:
                command += ' -score::elec_memb_sig_die '
                logger.log('ADD ELEC_MEMB_SIG_DIE')
        logger.log('running FilterScan for the energy function %s protocol, issuing command:\n%s' %
                   (energy_function, command))
        logger.log_text_file('%s' % str(PROTOCOLS_PATH + MPFilterScanDifferentSFAAs), to_print=print_xml)
        os.system(command)

        logger.log('ran FilterScan for energy function %s. finished at %s' %
                   (energy_function, time.strftime("%H:%M_%d%b")))

    # parse both sclog files to {pos: {AA: score}}
    if args['full']:
        shutil.move('temp.sclog', '%s.sclog' % energy_function)
        logger.log('saved the FilterScan sclog to %s.sclog' % energy_function)
    temp_fs_log = parse_filterscan_log('%s.sclog' % energy_function)

    # create DataFrame
    df = pd.DataFrame()
    for aa in AAs:
        # for pos in range(FLANK_SIZE+1, NUM_AAS+FLANK_SIZE+1+1):
        #     df = df.append({'pos': pos-FLANK_SIZE, 'aa': aa, energy_function: temp_fs_log[pos][aa]}, ignore_index=True) # switched for -50 to 50
        for pos in range(1, TOTAL_AAS+1):
            df = df.append({'pos': pos, 'aa': aa, energy_function: temp_fs_log[pos][aa]}, ignore_index=True)

    # calculate polyA scores in both score funcitons
    temp_A_mean = df[df['aa'] == 'A'][energy_function].mean()

    # calcualte delta of every mutant and the polyA for both score functions
    df['%s_normed' % energy_function] = df[energy_function]-temp_A_mean

    logger.log('mean of A for %s is %f' % (energy_function, temp_A_mean))

    ips = create_insertion_profiles(df, '%s_normed' % energy_function, adjust_extra_membranal, smooth=adjust_extra_membranal)
    return ips


def filterscan_parallel(energy_function: str, res_solv_weight: float, fa_cen: str, residues_to_test:
                        list=AAs, to_dump_pdbs: bool=False, adjust_extra_membranal: bool=True, 
                        print_xml: bool=False, iteration=0) -> dict:
    PWD = os.getcwd() + '/'

    # run FilterScan protocol to create sclog files for both full and no M env score functions
    if args['full']: # -unmute core.scoring.membrane.MPResSolvEnergy
        if energy_function is None:
            print('no energy function provided!')
            sys.exit()
    cmd = open('command', 'w+')
    if 'A' in residues_to_test:
        residues_to_test = [aa for aa in residues_to_test if aa != 'A']
    logger.log('creating jobs for residues %s, %i AAs' % (''.join(residues_to_test), len(residues_to_test)))
    num = 1
    iter_sclog_files = []
    for pos in POS_RANGE:
        for aa_ in residues_to_test:
            aa = aa_ if aa_ != residues_to_test[0] else residues_to_test[0] + 'A'
            name = '%s%i' % (aa, pos)
            command = '%s%s -database %s -parser:protocol %s -s %s -nstruct %i -overwrite -mp:scoring:hbond -mute all -ex1 -ex2 -ex3 -ex4 -parser:script_vars energy_function=%s residues_to_test=%s to_dump=%i res_solv_weight=%.2f fa_or_cen=%s res_num=%i sclog=%s cst_file=%s' % (ROSETTA_EXECUTABLES_PATH, ROSETTA_SCRIPTS_EXEC_PATH, ROSETTA_DATABASE_PATH, PROTOCOLS_PATH + 'MPFilterScanDifferentSFAAsPos.xml', PWD + POLY_A_NAME + '.pdb', NSTRUCT, energy_function, ''.join(aa), 1 if to_dump_pdbs else 0, res_solv_weight, fa_cen, pos, './pos_%s_%i.sclog' % (aa, pos), PWD+POLY_A_NAME+'.cst')
            iter_sclog_files.append('./pos_%s_%i.sclog' % (aa, pos))
            with open('job.%s' % name, 'w+') as job:
                job.write('#!/bin/bash\n. /usr/share/lsf/conf/profile.lsf\ncd %s\n'
                        % PWD)
                os.system('chmod +x job.%s' % name)
                job.write('%s\n' % command)
            res = subprocess.check_output('bsub -L /bin/bash -N -u /dev/null -G fleishman-wx-grp-lsf -q fleishman -o /dev/null 2>&1 -e /dev/null 2>&1 %sjob.%s' % (PWD, name), shell=True)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            # print('[%s%s]\r' % ('-'*num, ' '*(len(POS_RANGE)*len(residues_to_test)-num)))
            num += 1
            cmd.write('bsub -L /bin/bash -N -u /dev/null -G fleishman-wx-grp-lsf -q fleishman -o /dev/null 2>&1 -e /dev/null 2>&1 %sjob.%s\n' % (PWD, name))
    cmd.close()
    logger.log('finished submitting jobs')
    logger.log('running FilterScan for the energy function %s protocol, issuing command:' % energy_function)
    logger.log_text_file('%s' % str(PROTOCOLS_PATH + MPFilterScanDifferentSFAAs), to_print=print_xml)
    sclog_files = [a for a in os.listdir('./') if '.sclog' in a and 'pos' in a]
    while len(sclog_files) < len(POS_RANGE)*len(residues_to_test):
        logger.log('not enough finished %i' % len(sclog_files))
        time.sleep(10)
        sclog_files = [a for a in os.listdir('./') if '.sclog' in a and 'pos' in a]
    for sclog in iter_sclog_files:
        os.system('cat %s >> %s_iter_%i.sclog' % (sclog, energy_function, iteration))

    temp_fs_log = parse_filterscan_log('%s_iter_%i.sclog' % (energy_function, iteration))

    # create DataFrame
    df = pd.DataFrame()
    for aa in AAs:
        for pos in range(1, TOTAL_AAS+1):
            df = df.append({'pos': pos, 'aa': aa, energy_function: temp_fs_log[pos][aa]}, ignore_index=True)
    # calculate polyA scores in both score funcitons
    temp_A_mean = df[df['aa'] == 'A'][energy_function].mean()

    # calcualte delta of every mutant and the polyA for both score functions
    df['%s_normed' % energy_function] = df[energy_function]-temp_A_mean

    logger.log('mean of A for %s is %f' % (energy_function, temp_A_mean))

    ips = create_insertion_profiles(df, '%s_normed' % energy_function, adjust_extra_membranal)
    logger.log('erasing err.* out.* job.* command and pos_*sclog')
    os.system('rm -rf job.* err.* out.* command pos_*sclog')
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
    plt.close()


def parse_filterscan_log(file_name) -> dict:
    """
    parse a FilterScan run log file, returns dict {position: {AA: score}}
    """
    result = {i: {aa: None for aa in AAs} for i in range(1, TOTAL_AAS + 1)}
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
    # profiles = {aa: {pos: None for pos in range(FLANK_SIZE + 1, FLANK_SIZE + NUM_AAS + 1)} for aa in aas_names}
    profiles = {aa: {pos: None for pos in range(1, TOTAL_AAS + 1)} for aa in aas_names} # switched to this for -50 to 50
    sc_files = [a for a in os.listdir(all_configurations_path) if '.sc' in a]
    for sc_file in sc_files:
        sc = score2dict(all_configurations_path+sc_file)
        df = score_dict2df(sc)
        aa = [a for a in aas_names if a in sc_file][0]
        pos = int(sc_file.split('_')[2].split('.')[0])
        profiles[aa][pos] = np.mean(df['score']) - polyA_total_mean
        print('for res %s in pos %i found %.3f mean with %.3f std' % (aa, pos, np.mean(df['score']), np.std(df['score'])))
    membrane_position = np.linspace(-TOTAL_HALF_WIDTH, TOTAL_HALF_WIDTH, endpoint=True, num=NUM_AAS)
    mm_profiles = {aa: OrderedDict((membrane_position[pos - FLANK_SIZE - 1], profiles[aa][pos])
                                   for pos in range(FLANK_SIZE + 1, FLANK_SIZE + NUM_AAS + 1)) for aa in aas_names}
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
    logger.log('making T to be G with max at 1')
    hydrophobicity_polyval['THR'] = hydrophobicity_polyval['GLY'].copy()
    hydrophobicity_polyval['THR'][-1] = 1.0
    hydrophobicity_grade.close()
    return hydrophobicity_polyval


def draw_profiles(mm_profiles: dict, mm_eqs: dict):
    fig = plt.figure()
    elazar_profiles = get_elazar_scale()
    z = np.linspace(-TOTAL_HALF_WIDTH, TOTAL_HALF_WIDTH, endpoint=True, num=NUM_AAS + 1)
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
        plt.xlim([-(TOTAL_HALF_WIDTH + 0.5), TOTAL_HALF_WIDTH + 0.5])
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
        logger.log('created polyA file at %s.fa with %i As and %i flank size' %
                   (PWD + POLY_A_NAME, NUM_AAS, FLANK_SIZE))


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
    # cmd = '%shelix_from_sequence.default.linuxgccrelease -in:file:fasta %s -mp:setup:transform_into_membrane 1 -mute all' % (ROSETTA_EXECUTABLES_PATH, PWD + POLY_A_NAME + '.fa')
    cmd = '%shelix_from_sequence.default.linuxgccrelease -in:file:fasta %s -mute all' % (ROSETTA_EXECUTABLES_PATH, PWD + POLY_A_NAME + '.fa')
    logger.log('issuing command\n%s' % cmd)
    os.system(cmd)
    try:
        shutil.move(PWD + 'helix_from_sequence.pdb', PWD + POLY_A_NAME + '.pdb')
    except:
        shutil.move(PWD + 'S_0001.pdb', PWD + POLY_A_NAME + '.pdb')
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


def draw_elazar_splines(args) -> None:
    """
    make and draw Elazar splines
    :return:
    """
    profiles = MakeHydrophobicityGrade()
    x = np.arange(-50, +50, 0.1)
    Z = np.arange(-15, 15, 0.1)
    plt.figure()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.15, hspace=0.45)
    for i, aa in enumerate(list('ACDEFGHIKLMNPQRSTVWY')):
        y = [0] * (35 * 10) + [np.polyval(profiles[aas_1_3[aa]], z) for z in Z] + [0] * (35 * 10)
        if aa not in ['V', 'M', 'H']:
            tck = interpolate.splrep(x, y, s=25)
        else:
            tck = interpolate.splrep(x, y, s=5)
        xnew = np.arange(-50, 50, 0.05)
        ynew = interpolate.splev(xnew, tck, der=0)
        plt.subplot(5, 4, 1 + i)
        plt.plot(x, y, c='k')
        plt.plot(xnew, ynew, c='r')
        plt.scatter(x, [bspleval(x_, tck[0], tck[1], tck[2]) for x_ in x], c='g', marker='+')
        plt.vlines(-15, -2, 3, color='grey', linestyles='dashed')
        plt.vlines(15, -2, 3, color='grey', linestyles='dashed')
        plt.title(aa.upper())
        plt.ylim([-2, 3])
    plt.show()


def bspleval(x, knots, coeffs, order, debug=False):
    """
    adopted from http://scipy.github.io/old-wiki/pages/Numpy_Example_List_With_Doc.html

    Evaluate a B-spline at a set of points.

    Parameters
    ----------
    x : list or ndarray
        The set of points at which to evaluate the spline.
    knots : list or ndarray
        The set of knots used to define the spline.
    coeffs : list of ndarray
        The set of spline coefficients.
    order : int
        The order of the spline.

    Returns
    -------
    y : ndarray
        The value of the spline at each point in x.
    """

    k = order
    t = knots
    m = np.alen(t)
    npts = np.alen(x)
    B = np.zeros((m-1,k+1,npts))

    if debug:
        print('k=%i, m=%i, npts=%i' % (k, m, npts))
        print('t=', t)
        print('coeffs=', coeffs)

    ## Create the zero-order B-spline basis functions.
    for i in range(m-1):
        B[i,0,:] = np.float64(np.logical_and(x >= t[i], x < t[i+1]))

    if (k == 0):
        B[m-2,0,-1] = 1.0

    ## Next iteratively define the higher-order basis functions, working from lower order to higher.
    for j in range(1,k+1):
        for i in range(m-j-1):
            if (t[i+j] - t[i] == 0.0):
                first_term = 0.0
            else:
                first_term = ((x - t[i]) / (t[i+j] - t[i])) * B[i,j-1,:]

            if (t[i+j+1] - t[i+1] == 0.0):
                second_term = 0.0
            else:
                second_term = ((t[i+j+1] - x) / (t[i+j+1] - t[i+1])) * B[i+1,j-1,:]

            B[i,j,:] = first_term + second_term
        B[m-j-2,j,-1] = 1.0

    if debug:
        plt.figure()
        for i in range(m-1):
            plt.plot(x, B[i,k,:])
        plt.title('B-spline basis functions')

    ## Evaluate the spline by multiplying the coefficients with the highest-order basis functions.
    y = np.zeros(npts)
    for i in range(m-k-1):
        y += coeffs[i] * B[i,k,:]

    if debug:
        plt.figure()
        plt.plot(x, y)
        plt.title('spline curve')
        plt.show()

    return(y)


def understanding_spline_improvements():
    import random
    obj_f = InsertionProfile('A', pos_score={pos: pos*2 for pos in POS_RANGE}, adjust_extra_membranal=False)
    int_f = InsertionProfile('A', pos_score={pos: np.cosh(pos) for pos in POS_RANGE}, adjust_extra_membranal=False)

    for iter_num in range(10):
        for pos in POS_RANGE:
            int_f.pos_score[pos] += (obj_f.pos_score[pos] - int_f.pos_score[pos]) * random.randrange(-10, 10, 2)*0.001

    for pos in POS_RANGE:
        print(pos, obj_f.pos_score[pos], int_f.pos_score[pos], obj_f.pos_score[pos]-int_f.pos_score[pos])


def compare_multiple_splines():
    """
    compare splines from multiple energy funcs
    :return:
    """
    score_funcs = ['score0', 'score1', 'score2', 'score3', 'score5']
    score_splines = {k: dict() for k in score_funcs}
    for score in score_funcs:
        spline_file = sorted([a for a in os.listdir(score) if 'spline' in a])[0]
        for l in open(score+'/'+spline_file, 'r'):
            s = l.split()
            score_splines[score][s[0]] = [float(a) for a in s[1:]]

    elazar_ips = create_elazar_ips()

    plt.figure()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.15, hspace=0.45)
    for i, aa in enumerate(AAs):
        plt.subplot(5, 4, 1 + i)
        plt.plot(Z_total, [elazar_ips[aa].pos_score[pos] for pos in POS_RANGE], label='elazar')
        for sc in score_funcs:
            plt.plot(Z_total, score_splines[sc][aa], label=sc)
        if aa != 'P':
            plt.ylim([-5, 5])
        plt.title(aa)
    plt.legend()
    plt.savefig('spline_comparisons.png')
    plt.show()


def create_resfile(file_name):
    with open(file_name, 'w+') as fout:
        fout.write('PIKAA\nstart\n')
        for pos in POS_RANGE:
            fout.write('%i\tA\tPIKAA\t ACDEFGHIKLMNPQRSTVWY\n' % pos)


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
    parser.add_argument('-energy_func_fa', default='talaris2014_elazaridis')
    parser.add_argument('-energy_func_cen')
    parser.add_argument('-note', default=None, type=str, help='add note to spline files')
    parser.add_argument('-use_made_pdb', default=False)
    parser.add_argument('-elec_memb_sig_die', default=False)

    logger = Logger('elazaridis_%s.log' % time.strftime("%H_%M_%d_%b"))

    args = vars(parser.parse_args())
    
    if args['mode'] == 'main':
        global PWD
        PWD = os.getcwd() + '/'
        calibrate_function()

    elif args['mode'] == 'calibrate_multiple':
        calibrate_energy_functions(args)

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

    elif args['mode'] == 'draw_rosetta_profiles_fa_cen':
        draw_rosetta_profiles_fa_cen(args)

    elif args['mode'] == 'draw_elazar_splines':
        draw_elazar_splines(args)

    elif args['mode'] == 'create_e_term_specific_profiles':
        create_e_term_specific_profiles(args, './', ['fa_rep', 'fa_mpsolv'])

    elif args['mode'] == 'just':
        just_draw_current_profiles()

    elif args['mode'] == 'compare_multiple_splines':
        compare_multiple_splines()

    elif args['mode'] == 'test':
        create_resfile('test.txt')

    else:
        print('unknown mode')

    logger.close()
