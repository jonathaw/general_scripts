#!/usr/bin/env python3.5
"""
"""
import os
import sys
import copy
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import MyPDB as mp
import MyPDB_funcs as mpf
import RosettaFilter as rf
from draw_dielectric_membrane_mesh import create_xml_flags


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode')
    parser.add_argument('-stage')
    parser.add_argument('-z', default=0, type=float)
    parser.add_argument('-in_membrane', type=bool, default=True)
    parser.add_argument('-memb_sig', type=bool, default=True)
    args = vars(parser.parse_args())

    ds = np.append( np.arange(0.1, 5, 0.1), np.arange(5, 50, 5) )
    zs = np.append( np.arange(0.0, 2.5, 0.1), np.arange(2.5, 25, 2.5) )

    ds = np.arange( 0.1, 50, 0.1 )
    zs = np.arange( 0.0, 25, 1 )
    # zs = [ 0.0 ]

    if args['mode'] == 'setup':
        prepare_scan( args, ds, zs )

    elif args['mode'] == 'analyse':
        analyse_scan( args, ds, zs )

    elif args['mode'] == 'test':
        create_pdb_AA_AA_d_z('D', 'R', 1.4, 0, ResMaker())

    else:
        print('no mode given')


def analyse_scan( args, ds, zs ):
    """
    analyse the scores generated for the PDBs
    """
    df = read_all_scores('scores/all_scores.score')
    plt.figure()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.15, hspace=0.45)
    i = 0
    for aa1 in ['D', 'S']:
        for aa2 in ['K', 'R']:
            plt.subplot(2, 2, 1+i)
            mini_df = df[ (df['aa1'] == aa1) & (df['aa2'] == aa2) & (df['z'] == args['z']) ]
            mini_df.sort_values( 'd', inplace=1 )
            mini_df.to_csv('%s_%s_z%.2f.csv' % (aa1, aa2, args['z']),
                           columns=['aa1', 'aa2', 'd', 'z', 'fa_elec'],
                           sep='\t', index=False)
            # plt.scatter( mini_df['d'], mini_df['fa_elec'])
            plt.plot( mini_df['d'], mini_df['fa_elec'])

            plt.xlim([-1, np.max(ds)])
            plt.ylim([np.min( mini_df['fa_elec'].values ) - 1, 1])
            plt.xlabel('d')
            plt.ylabel('fa_elec')
            plt.title('%s vs. %s at z=%i' % ( aa1, aa2, args['z'] ))

            i += 1
    plt.show()


def read_all_scores(file_name: str) -> pd.DataFrame:
    """
    if necessary gather scores, and coalesce into dataframe with resiudes, ds and zs
    """
    if not os.path.exists('scores/all_scores.sc'):
        os.system('grep description scores/%s > scores/all_scores.score' %
                  [a for a in os.listdir('scores/') if '.sc' in a][0])
        os.system('grep SCORE: scores/*.sc | grep -v description >> scores/all_scores.score')
    df = rf.score_file2df( file_name )
    desc_spl = df['description'].str.split('_')
    df['aa1'] = desc_spl.str.get( 0 )
    df['aa2'] = desc_spl.str.get( 1 )
    df['d'] = desc_spl.str.get( 2 ).astype( float )
    df['z'] = desc_spl.str.get( 3 ).astype( float )
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
                        cmd_line = 'bsub -L /bin/bash -N -u /dev/null -G fleishman-wx-grp-lsf -q fleishman -o /dev/null 2>&1 -e /dev/null 2>&1 %s/jobs/job.%s\n' % (pwd, name)
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


if __name__ == '__main__':
    main()
