#!/usr/bin/env python3.5
"""
"""
import os
import sys
import copy
import time
import argparse
import numpy as np
import scipy as sc
import pandas as pd
from matplotlib import cm
import matplotlib.pyplot as plt

import MyPDB as mp
import MyPDB_funcs as mpf
import RosettaFilter as rf
from Logger import Logger
from draw_dielectric_membrane_mesh import create_xml_flags
from RosettaData import rosetta_atom_radii
from draw_dielectric_membrane_mesh import rosetta_dz_model, parse_rosetta_log


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode')
    parser.add_argument('-stage')
    parser.add_argument('-sc', type=str, default='scores/all_scores.score')
    parser.add_argument('-z', default=0, type=float)
    parser.add_argument('-in_membrane', type=bool, default=True)
    parser.add_argument('-memb_sig', type=bool, default=True)
    parser.add_argument('-log')
    parser.add_argument('-plot_type', default='1ds')
    parser.add_argument('-dslk_scores')
    args = vars(parser.parse_args())

    args['logger'] = Logger('logeer_%s_%s.log' % (args['mode'], time.strftime("%d.%0-m")))
    args['res_maker'] = ResMaker()
    ds = np.append( np.arange(0.1, 5, 0.1), np.arange(5, 50, 5) )
    zs = np.append( np.arange(0.0, 2.5, 0.1), np.arange(2.5, 25, 2.5) )

    ds = np.arange( 0.1, 50, 0.1 )
    zs = np.arange( 0.0, 25, 1 )
    # zs = [ 0.0 ]

    if args['mode'] == 'setup':
        prepare_scan( args, ds, zs )

    elif args['mode'] == 'analyse_z0':
        analyse_scan_z0( args, ds, zs )

    elif args['mode'] == 'analyse_dz':
        analyse_dz( args, ds, zs )

    elif args['mode'] == 'test':
        create_pdb_AA_AA_d_z('D', 'R', 1.4, 0, ResMaker())

    elif args['mode'] == 'caca':
        CaCa_d_by_z(args)

    elif args['mode'] == 'just_model':
        draw_just_the_model( args )

    else:
        print('no mode given')


def draw_just_the_model( args: dict ):
    ds = np.arange(0.1, 15, 0.1)
    zs = np.arange(0, 25, 0.5)

    # D, Z, E = [], [], []
    df = pd.DataFrame(columns=['d', 'z', 'e'])
    for d in ds:
        for z in zs:
            # D.append( d )
            # Z.append( z )
            # E.append( rosetta_dz_model( d, z ))
            df = df.append({'d': d, 'z': z, 'e': rosetta_dz_model( d, z )}, ignore_index=1)
    fig = plt.figure()
    fig.subplots_adjust(hspace=.5)
    i = 0
    for z in np.arange(0, 20, 1):
        plt.subplot(4, 5, 1+i)
        plt.plot( df[ df['z' ] == z]['d'], df[ df['z'] == z ]['e'] )
        plt.title('z=%i' % z)
        plt.xlim([0, 16])
        i += 1
    plt.show()


def analyse_dz( args: dict, ds: np.arange, zs: np.arange ) -> None:
    df = read_all_scores( 'scores/all_scores.score', args )
    df = df[ df['d'] < 10 ]

    fig = plt.figure()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.15, hspace=0.45)
    i = 0
    for aa1 in ['D', 'S']:
        for aa2 in ['K', 'R']:
            plt.subplot(2, 2, 1+i)
            mini_df = df[ (df['aa1'] == aa1) & (df['aa2'] == aa2) ]
            mini_df.sort_values('d', inplace=1)
            args['logger'].create_header('analysing %s - %s' % ( aa1, aa2 ))
            for n in ['d', 'z', 'fa_elec']:
                args['logger'].log('for %s there are %i values, mean %.2f, min %.2f, max %.2f' %
                                   (n, len(mini_df[n]), mini_df[n].mean(), mini_df[n].min(), mini_df[n].max()))
                args['logger'].log('the set for %s is %r' % (n, np.round(sorted(list(set(mini_df[n].values))), 2)))

            D, Z, E = [], [], []
            for d in ds:
                if d >= 10:
                    continue
                for z in zs:
                    # if z > 5: continue
                    D.append( d )
                    Z.append( z )
                    E.append( mini_df[ (mini_df['d'] >= d-0.01) &
                                       (mini_df['d'] <= d+0.01) &
                                       (mini_df['z'] == z) ]['fa_elec'].values[0] )
            args['logger'].log('finished transforming dataframe to lists')

            D_, Z_ = np.meshgrid( D, Z )

            levels = np.linspace( np.min(E), np.max(E), 500 )
            E_ = sc.interpolate.griddata( (D, Z), E, (D_, Z_), method='cubic' )
            cs = plt.contourf( D_, Z_, E_, levels=levels, cmap=cm.coolwarm)
            plt.xlabel('d (A)')
            plt.ylabel('z (A)')
            plt.title('%s vs. %s' % ( aa1, aa2 ))
            fig.colorbar(cs, format="%.2f")
            i += 1
            # break
        # break
    plt.show()


def CaCa_d_by_z( args: dict ) -> None:
    # reading the Ca-Ca reuslts, from args['sc']
    df = rf.score_file2df( args['sc'] )
    desc = df['description'].str.split('_')
    df['d'] = desc.str.get( 0 ).astype( np.float64 )
    df['z'] = desc.str.get( 1 ).astype( np.float64 )

    # get data from spline log
    spline_log = parse_rosetta_log( args )
    spline_log_df = pd.DataFrame({'z': spline_log[0], 'd': spline_log[1], 'e': spline_log[2]})

    # get DSRK data as well
    dslk_df = read_all_scores( args['dslk_scores'], args )

    ds_sorted = np.round(sorted(list(set(df['d'].values))), 2)
    zs_sorted = np.round(sorted(list(set(df['z'].values))), 2)
    D, Z, E = [], [], []
    model_E = []
    spline_log_E = []
    dre, dke, sre, ske = [], [], [], []
    for d in ds_sorted:
        if d >= 10:
            continue
        for z in zs_sorted:
            # if z > 5: continue
            D.append( d )
            Z.append( z )
            E.append( df[ (df['d'] >= d-0.01) &
                          (df['d'] <= d+0.01) &
                          (df['z'] == z) ]['fa_elec'].values[0] / 4 )
            model_E.append( rosetta_dz_model( d, z ))
            spline_log_E.append( spline_log_df[ (spline_log_df['d'] >= d-0.01) &
                                               (spline_log_df['d'] <= d+0.01) &
                                              (spline_log_df['z'] == z)]['e'].values[0] )

    # gather DSLK data in lists to make plots later
    dslk_D, dslk_Z = [], []
    dre, dke, sre, ske = [], [], [], []
    for d in np.round(sorted(list(set(dslk_df['d'].values))), 2):
        for z in np.round(sorted(list(set(dslk_df['z'].values))), 2):
            dslk_D.append( d )
            dslk_Z.append( z )
            dre.append( dslk_df[ (dslk_df['aa1']=='D')&(dslk_df['aa2']=='R')&(dslk_df['d']==d)&(dslk_df['z']==z) ]['fa_elec'].values[0])
            dke.append( dslk_df[ (dslk_df['aa1']=='D')&(dslk_df['aa2']=='K')&(dslk_df['d']==d)&(dslk_df['z']==z)  ]['fa_elec'].values[0] )
            sre.append( dslk_df[ (dslk_df['aa1']=='D')&(dslk_df['aa2']=='R')&(dslk_df['d']==d)&(dslk_df['z']==z)  ]['fa_elec'].values[0] )
            ske.append( dslk_df[ (dslk_df['aa1']=='D')&(dslk_df['aa2']=='K')&(dslk_df['d']==d)&(dslk_df['z']==z)  ]['fa_elec'].values[0] )


    args['logger'].log('finished preparing lists')
    master_df = pd.DataFrame({'d': D, 'z': Z, 'e_caca': E,
                              'model_e': model_E, 'spline_log_e': spline_log_E})

    if args['plot_type'] == '1ds':
        args['logger'].log('making 1d plots')
        fig = plt.figure()
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.15, hspace=0.45)
        i = 0
        for i, z in enumerate(zs_sorted[::2]):
            plt.subplot(5, 4, 1+i)
            d_df = master_df[ master_df['z'] == z ]
            plt.plot( d_df['d'], d_df['e_caca'].values, label='e_caca', c='r' )
            plt.plot( d_df['d'], d_df['model_e'].values, label='model_e', c='b' )
            plt.plot( d_df['d'], d_df['spline_log_e'].values, label='spline_log_e', c='g' )
            plt.title('z=%.2f' % z)
            plt.xlabel('d')
            plt.ylabel('e')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.show()



    if args['plot_type'] == 'contourf':
        D_, Z_ = np.meshgrid( D, Z )

        levels = np.linspace( np.min(E), np.max(E), 1000 )

        E_ = sc.interpolate.griddata( (D, Z), E, (D_, Z_), method='cubic')
        model_E_ = sc.interpolate.griddata( (D, Z), model_E, (D_, Z_), method='cubic')
        spline_log_E_ = sc.interpolate.griddata( (D, Z), spline_log_E, (D_, Z_), method='cubic' )

        fig = plt.figure()
        plt.subplot(1, 2, 1)
        cs = plt.contourf( D_, Z_, E_, levels=levels, cmap=cm.coolwarm )
        plt.xlabel('d (A)')
        plt.ylabel('z (A)')
        fig.colorbar(cs, format="%.2f")

        plt.subplot(1, 2, 2)
        cs = plt.contourf( D_, Z_, model_E_, levels=levels, cmap=cm.coolwarm  )
        plt.xlabel('d (A)')
        plt.ylabel('z (A)')
        fig.colorbar(cs, format="%.2f")
        plt.show()


def analyse_scan_z0( args, ds, zs ):
    """
    analyse the scores generated for the PDBs
    """
    df = read_all_scores('scores/all_scores.score', args)
    plt.figure()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.15, hspace=0.45)
    i = 0
    for aa1 in ['D', 'S']:
        for aa2 in ['K', 'R']:
            plt.subplot(2, 2, 1+i)
            args['logger'].create_header('working on %s - %s' % (aa1, aa2))
            d_min = args['res_maker'].get_maximal_main_radii( aa1 ) + args['res_maker'].get_maximal_main_radii( aa2 )
            args['logger'].log('the sum of radii for pair is %.2f' % d_min)
            mini_df = df[ (df['aa1'] == aa1) & (df['aa2'] == aa2) & (df['z'] == args['z']) & (df['d'] > d_min) ]
            mini_df.sort_values( 'd', inplace=1 )
            mini_df.to_csv('%s_%s_z%.2f.csv' % (aa1, aa2, args['z']),
                           columns=['aa1', 'aa2', 'd', 'z', 'fa_elec'],
                           sep='\t', index=False)
            # plt.scatter( mini_df['d'], mini_df['fa_elec'])
            plt.plot( mini_df['d'], mini_df['fa_elec'])

            plt.xlim([0, 10])#np.max(ds)])
            plt.ylim([np.min( mini_df['fa_elec'].values ) - 1, 1])
            plt.xlabel('d')
            plt.ylabel('fa_elec')
            plt.title('%s vs. %s at z=%i' % ( aa1, aa2, args['z'] ))

            i += 1
    plt.show()


def read_all_scores( file_name: str, args=dict() ) -> pd.DataFrame:
    """
    if necessary gather scores, and coalesce into dataframe with resiudes, ds and zs
    """
    # if not os.path.exists(args['sc']):
        # args['logger'].logger('bo score file found. gathering scores to %s' % args['sc'])
        # os.system('grep description scores/%s > %s' %
                  # ( args['sc'], [a for a in os.listdir('scores/') if '.sc' in a][0] ))
        # os.system('grep SCORE: scores/*.sc | grep -v description >> %s' % args['sc'])
    df = rf.score_file2df( file_name )
    args['logger'].log('found %i entries in score file' % len( df ))
    desc_spl = df['description'].str.split('_')
    df['aa1'] = desc_spl.str.get( 0 )
    df['aa2'] = desc_spl.str.get( 1 )
    df['d'] = desc_spl.str.get( 2 ).astype( np.float64 )
    df['z'] = desc_spl.str.get( 3 ).astype( np.float64 )
    df = df.round( {'d': 2, 'z': 2} )
    return df


def prepare_scan( args, ds, zs ):
    """
    setup a scan, making pdbs, and jobs
    """
    res_maker = ResMaker()
    create_xml_flags( args['in_membrane'], args['memb_sig'] )
    make_pdbs_d_z_aa_aa( args, ds, zs, res_maker )
    make_jobs_d_z_aa_aa( args, ds, zs )


def make_jobs_d_z_aa_aa( args, ds, zs ):
    """
    create jobs to calculate rosetta scores for pdbs at ds, and zs
    """
    in_membrane = True
    os.mkdir('jobs')
    os.mkdir('scores')
    os.mkdir('pdbs_0001')
    cmd = open('jobs/command', 'w+')
    pwd = os.getcwd()
    for aa1 in ['D', 'S']:
        for aa2 in ['K', 'R']:
            for d in ds:
                for z in zs:
                    name = '%s_%s_%.2f_%.2f' % ( aa1, aa2, d, z  )
                    with open('jobs/job.%s' % name, 'w+') as job:
                        job.write('#!/bin/bash\n. /usr/share/lsf/conf/profile.lsf\ncd %s\n' % pwd)
                        job.write('/home/labs/fleishman/jonathaw/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease  -s pdbs/%s.pdb @%s/%s -out:file:scorefile scores/%s.sc -out:path:pdb pdbs_0001/ -script_vars scfxn=beta_nov16%s\n' % (name, pwd, 'in_memb.flags' if in_membrane else 'no_memb.flags', name, '_elazaridis' if in_membrane else ''))

                        os.system('chmod +x jobs/job.%s' % name)
                        cmd_line = 'bsub -N -u /dev/null -G fleishman-wx-grp-lsf -q fleishman -o /dev/null 2>&1 -e /dev/null 2>&1 %s/jobs/job.%s\n' % (pwd, name)
                        os.system( cmd_line )
                        cmd.write( cmd_line )
    cmd.close()


def make_pdbs_d_z_aa_aa( args, ds, zs, res_maker ):
    """
    prepare pdbs for DS/KR at ds and zs
    """
    os.mkdir('pdbs')
    for aa1 in ['D', 'S']:
        for aa2 in ['K', 'R']:
            for d in ds:
                for z in zs:
                    create_pdb_AA_AA_d_z( aa1, aa2, d, z , res_maker, 'pdbs')


def create_pdb_AA_AA_d_z( aa1: str, aa2: str, d: float, z: float, res_maker, path:str='./' ) -> None:
    """
    create a pdb with residues aa1 and aa2 at the XY plane at Z=z and distance d
    """
    res1 = res_maker.get_residue( aa1 )
    res2 = res_maker.get_residue( aa2 )

    res1.change_chain_name('A')
    res2.change_chain_name('B')

    # rotate around Z to oppose
    axis_z = mp.XYZ(0, 0, 1)
    res2.dot_matrix_me( mpf.rotation_matrix_around_vec( axis_z, np.pi ) )

    # translate to get d distance
    move_d = mp.XYZ(0, res1[ res_maker.main_residue_atoms[res1.res_type][1] ].xyz.y -
                            res2[ res_maker.main_residue_atoms[res2.res_type][1] ].xyz.y + d, 0)
    res2.translate_xyz( move_d )

    # translate all to z
    move_z = mp.XYZ(0, 0, z)
    res1.translate_xyz( move_z  )
    res2.translate_xyz( move_z  )

    # setup in a MyPDB instance, and renumber and write
    pdb = mp.MyPDB()
    for res in [res1, res2]:
        for a in res.values():
            pdb.add_atom( a )
    pdb.renumber()
    mp.write_PDB( '%s/%s_%s_%.2f_%.2f.pdb' % ( path, res1.res_type, res2.res_type, d, z ), pdb )


class ResMaker:
    def __init__( self ):
        self.residues = {}
        self.set_main_residues_atoms()
        self.parse_residues_file()
        self.atom_radii = rosetta_atom_radii

    def set_main_residues_atoms( self ) -> None:
        result = {}
        result['D'] = ['CG', 'OD1', 'OD2']
        result['S'] = ['CB', 'OG', 'HG']
        result['R'] = ['CZ', '2HH1', '2HH2']
        result['K'] = ['NZ', '1HZ', '2HZ']
        self.main_residue_atoms = result

    def get_residue( self, aa: str ) -> mp.Residue:
        return copy.deepcopy( self.residues[ aa ] )

    def parse_residues_file( self ) -> dict:
        pdb = mp.parse_PDB('/home/labs/fleishman/jonathaw/temp_residue_data/RKDS.txt')
        for res in pdb['A'].values():
            mpf.translate_and_rotate_res_to_xy_plane( res, self.main_residue_atoms[res.res_type] )
            self.residues[res.res_type] = res

    def get_main_atom_types( self, aa: str, atom_num: int ) -> str:
        atom = self.main_residue_atoms[ aa ][ atom_num ]
        if atom[0] in ['1', '2']:
            return atom[1]
        else:
            return atom[0]

    def get_main_atom_radius( self, aa: str, atom_num: int ) -> float:
        return self.atom_radii[ self.get_main_atom_types( aa, atom_num ) ]

    def get_maximal_main_radii( self, aa ) -> float:
        return max( [ self.get_main_atom_radius( aa, 1 ), self.get_main_atom_radius( aa, 2 )] )


if __name__ == '__main__':
    main()
