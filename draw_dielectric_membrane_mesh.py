#!/usr/bin/env python3.5
import os
import sys
import time
import argparse

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import pandas as pd

from scipy.interpolate import griddata
from Logger import Logger
import MyPDB as mp
import RosettaFilter as rf
import numpy as np


def rosetta_coefficients():
    low_poly_start2_ = 1.8225
    low_poly_end2_ = 3.4225
    hi_poly_start2_ = 20.25
    max_dis2_ = 30.25
    C0_ = 322.064
    C1_ = 322.064
    C2_ = 1.72637
    Min_dis_score_ = 23.3829
    dEfac_ = -322.064
    sigmoidal_die_ = 1
    smooth_fa_elec_ = 1
    no_dis_dep_die_ = 0
    sigmoidal_D_ = 80
    sigmoidal_D0_ = 6
    sigmoidal_S_ = 0.4


def beta_nov16_rosetta_coefficents() -> dict:
    return {'C0_': 322.064, 'C1_': 322.064, 'C2_': 1.51278,
            'min_dis_score_': 20.347,
            'max_dis_': 5.5,
            'min_dis_': 1.6,
            'dEfac_': -322.064,
            'low_poly_start_': 1.35,
            'low_poly_end_': 1.85,
            'sigmoidal_die_': 1,
            'no_dis_dep_die_': 0,
            'hi_poly_start2_': 20.25,
            'sigmoidal_D_': 79.9,
            'sigmoidal_D0_': 6.65,
            'sigmoidal_S_': 0.415,
            'low_poly_': {'c0': -54.9047, 'c1': 132.607, 'c2': -72.5834, 'c3': 11.59},
            'hi_poly_': {'c0': -43.5752, 'c1': 31.3588, 'c2': -7.0817, 'c3': 0.512837}}



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-log')
    parser.add_argument('-mode', default='3d')
    parser.add_argument('-stage', default=0, type=int)
    parser.add_argument('-memb_sig', type=bool, default=True)
    args = vars(parser.parse_args())

    args['logger'] = Logger('logeer_%s_%i_%s.log' % (args['mode'], args['stage'], time.strftime("%d.%0-m"))
                           )
    if args['mode'] == '3d':
        draw_3d_mesh_by_log(args)

    elif args['mode'] == '2d':
        draw_2d_by_log(args)

    elif args['mode'] == 'draw_rosetta_d_model':
        draw_rosetta_d_model(args)

    elif args['mode'] == 'draw_rosetta_d_by_z_model':
        draw_rosetta_d_by_z_model(args)

    elif args['mode'] == 'draw_desired_model':
        draw_desired_model(args)

    else:
        print('no mode given')


def draw_desired_model(args):
    rosetta_coefs = beta_nov16_rosetta_coefficents()
    ds = np.arange(0, 50, 0.5)
    zs = np.arange(0, 25, 0.5)
    ds = np.append( np.arange(0.1, 5, 0.5), np.arange(5, 50, 5) )
    zs = np.append( np.arange(0.1, 2.5, 0.5), np.arange(2.5, 25, 2.5) )
    if args['stage'] == 1:
        get_rosetta_results_for_d_z( np.arange(0, 50, 5), np.arange(0, 25, 2.5), args['stage'] )
        print('when finished submitting jobs, run the same command with -stage 2')
        sys.exit()
    max_chr_arg_n = -2 # -0.8

    df = pd.DataFrame(columns=['d', 'z', 'e', 'f_d', 'g_z'])
    D, Z, E = [], [], []
    fs_of_ds, gs_of_zs = [], []
    for d in ds[1:]:
        for z in zs:
            D.append(d)
            Z.append(z)
            f_d = calc_rosetta_d_func( d, rosetta_coefs )
            g_z = (rosetta_coefs['min_dis_score_'] * rosetta_coefs['low_poly_start_'] *
                   ( 1 - (( z / 10 )**4 / ( 1 + ( z / 10 ) ** 4 ) ) ) / d)
            if d <= rosetta_coefs['low_poly_start_']:
                e = rosetta_coefs['min_dis_score_']
            else:
                e = np.max( [ f_d, g_z  ]  )
            df = df.append({'d': d, 'z': z, 'f_d': f_d, 'g_z': g_z, 'e': e}, ignore_index=True)
            fs_of_ds.append( f_d )
            gs_of_zs.append( g_z )
            E.append( e )
    # plt.plot(df[ df['z'] == 0 ]['d'], df[ df['z'] == 0 ]['e'] * ( max_chr_arg_n ** 2 ), label='z=0')
    # plt.plot(df[ df['z'] == 24.5  ]['d'], df[ df['z'] == 24.5 ]['e'] * ( max_chr_arg_n ** 2 ), label='z=24.5')
    # plt.axhline()
    # plt.legend()
    df.to_csv('model_dze.csv')

    actual_df = get_rosetta_results_for_d_z( ds, zs, args['stage'] )
    actual_df.to_csv('actual_rosetta_dze.csv')

    # print(E)
    # print(len(D), len(Z), len(E))
    E = [e * ( max_chr_arg_n ** 2 ) for e in E]
    D_, Z_ = np.meshgrid( D, Z )
    E_ = griddata((D, Z), E, (D_, Z_), 'cubic')
    # print(D_)
    # print(Z_)
    # print(E_)
    # print(len(D_), len(Z_), len(E_))
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    # surf = ax.plot_surface(D_, Z_, E_, rstride=50, cstride=50, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.scatter(D, Z, E, c='r')
    ax.scatter(actual_df['d'], actual_df['z'], actual_df['fa_elec'], c='k')
    ax.set_xlabel('D')
    ax.set_ylabel('Z')
    ax.set_zlabel('E')
    plt.show()


def get_rosetta_results_for_d_z( ds, zs, stage ) -> pd.DataFrame:
    if stage == 1:
        get_rosetta_dielectric_constants_for_range(dict(), ds, zs, in_membrane=True, memb_sig=True)
    elif stage == 2:
        if not os.path.exists('scores/all_scores.sc'):
            gather_scores()
        ds_, zs_, es_ = get_dielectric_from_rosetta()
        df = pd.DataFrame(columns=['d', 'z', 'fa_elec'])
        for d, z, e in zip( ds_, zs_, es_  ):
            if d == 0: continue
            df = df.append({'d': d, 'z': z, 'fa_elec': e}, ignore_index=True)
        return df
    else:
        print('unknown stage')
        sys.exit()


def calc_rosetta_d_func( d, c ) -> float:
    d2 = d ** 2

    if d2 > c['max_dis_'] ** 2:
        return 0.0

    elif d2 < c['low_poly_start_'] ** 2:
        return c['min_dis_score_']

    elif d2 < c['low_poly_end_'] ** 2:
        return (  (( c['low_poly_']['c3']*d + c['low_poly_']['c2'])*d+c['low_poly_']['c1'])*d + c['low_poly_']['c0'] )

    elif d2 > c['hi_poly_start2_']:
        return (  (( c['hi_poly_']['c3']*d + c['hi_poly_']['c2'] )*d+c['hi_poly_']['c1'])*d + c['hi_poly_']['c0'] )

    elif c['sigmoidal_die_']:
        return ( c['C1_'] / ( d * sigmoid_eps( d, c ) ) - c['C2_'] )

    elif c['no_dis_dep_die_']:
        return ( c['C1_'] / np.sqrt( d2 ) - c['C2_'] )

    else:
        return ( c['C1_'] / d2 - c['C2_'] )


def sigmoid_eps( d: float, c: dict ) -> float:
    return ( c['sigmoidal_D_'] - 0.5 * ( c['sigmoidal_D_'] - c['sigmoidal_D0_'] ) * ( 2 + 2 * d * c['sigmoidal_S_'] + d * d * c['sigmoidal_S_'] * c['sigmoidal_S_'] ) * np.exp( - d * c['sigmoidal_S_'] ) )


def draw_rosetta_d_by_z_model(args):
    ds = np.arange(0, 15, 0.1)
    zs = np.arange(0, 20, 0.5)
    if args['stage'] == 1:
        args['logger'].log('creating pdbs and jobs for \nds %r\nzs %r' % (ds, zs))
        get_rosetta_dielectric_constants_for_range(args, ds, zs, in_membrane=True, memb_sig=args['memb_sig'])
    elif args['stage'] == 2:
        if not os.path.exists('scores/all_scores.sc'):
            gather_scores()
        ds, zs, es = get_dielectric_from_rosetta()
        args['logger'].log('correcting energies by 4')
        es = [a/4.0 for a in es]

        for name, vals in zip(['d', 'z', 'e'], [ds, zs, es]):
            print('for %s mean: %.2f min %.2f max %.2f count %i' % (name, np.mean(vals), np.min(vals), np.max(vals), len(vals)))
        data = [{'d': d, 'z': z, 'e': e} for d, z, e in zip(ds, zs, es)]
        df = pd.DataFrame(data, columns=['d', 'z', 'e'])
        df.to_csv('dze.csv', sep="\t")

        D, Z = np.meshgrid( ds, zs )
        E = griddata((ds, zs,), es, (D, Z), method='cubic')
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(D, Z, E, rstride=50, cstride=50, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        # surf = ax.scatter(D, Z, E, c=E)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        ax.set_xlabel('D')
        ax.set_ylabel('Z')
        ax.set_zlabel('fa_elec')
        plt.show()
    else:
        print('invalid stage, use 1 or 2')


def draw_rosetta_d_model(args):
    ds = np.arange(0, 15, 0.1)
    if args['stage'] == 1:
        es = get_rosetta_dielectric_constants_for_range(args, ds, [0], False, args)
    elif args['stage'] == 2:
        if not os.path.exists('scores/all_scores.sc'):
            gather_scores()
        ds, zs, es = get_dielectric_from_rosetta()
        plt.scatter(ds, es)
        plt.xlabel('distance A')
        plt.ylabel('fa_elec')
        plt.show()
    else:
        print('invalid stage, use 1, than 2')


def get_rosetta_dielectric_constants_for_range(args, ds: list, zs: list, in_membrane: bool, memb_sig: bool=False) -> list:
    # args['logger'].log('creating jobs, pdbs, flags etc with in_membrane %r nad memb_sig %r' % (in_membrane, memb_sig))

    pwd = os.getcwd()
    os.mkdir('pdbs')
    os.mkdir('jobs')
    os.mkdir('scores')
    os.mkdir('pdbs_0001')
    cmd = open('jobs/command', 'w+')
    for z in zs:
        for d in ds:
            name = '%.2f_%.2f' % (d, z)
            create_pdb_from_d_z(d, z)
            with open('jobs/job.%s' % name, 'w+') as job:
                job.write('#!/bin/bash\n. /usr/share/lsf/conf/profile.lsf\ncd %s\n' % pwd)
                job.write('/home/labs/fleishman/jonathaw/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease  -s pdbs/%.2f_%.2f.pdb @%s/%s -out:file:scorefile scores/%s.sc -out:path:pdb pdbs_0001/ -script_vars scfxn=beta_nov16%s\n' % (d, z, pwd, 'in_memb.flags' if in_membrane else 'no_memb.flags', name, '_elazaridis' if in_membrane else ''))

                os.system('chmod +x jobs/job.%s' % name)
                cmd.write('bsub -L /bin/bash -N -u /dev/null -G fleishman-wx-grp-lsf -q fleishman -o /dev/null 2>&1 -e /dev/null 2>&1 %s/jobs/job.%s\n' % (pwd, name))
    cmd.close()
    os.system('sh jobs/command')
    create_xml_flags( in_membrane, memb_sig )
    os.system('sh jobs/command')


def create_xml_flags( in_membrane: bool, memb_sig: bool ):
    name = 'in_memb' if in_membrane else 'no_memb'
    mover = """\n\t\t<AddMembraneMover name="add_memb" membrane_core="10" steepness="4">\n\t\t\t<Span start="1" end="2" orientation="in2out"/>\n\t\t</AddMembraneMover>""" if in_membrane else ''
    protocol = '\n\t\t<Add mover="add_memb"/>\n' if in_membrane else ''
    with open('%s.flags' % name, 'w+') as fo:
        fo.write("-parser:protocol %s_memb.xml\n-overwrite\n-corrections::beta_nov16\n-score::memb_fa_sol\n%s" % ('in' if in_membrane else 'no', '-score::elec_memb_sig_die' if memb_sig else ''))
    with open('%s.xml' % name, 'w+') as fo:
        fo.write("""<ROSETTASCRIPTS>
   <TASKOPERATIONS></TASKOPERATIONS>
   <SCOREFXNS>
     <ScoreFunction name="scfxn" weights="%%%%scfxn%%%%">
      <Reweight scoretype="fa_elec" weight="1"/>
    </ScoreFunction>
  </SCOREFXNS>
  <MOVERS>%s
 </MOVERS>
  <FILTERS>
  </FILTERS>
  <PROTOCOLS>%s
  </PROTOCOLS>
  <OUTPUT scorefxn="scfxn"/>
</ROSETTASCRIPTS>\n""" %  (mover, protocol))


def gather_scores():
    os.system('grep description scores/0.00_0.00.sc > scores/all_scores.score')
    os.system('grep SCORE: scores/*.sc | grep -v description >> scores/all_scores.score')


def get_dielectric_from_rosetta():
    df = rf.score_file2df('scores/all_scores.score')
    ds = [float(a.split('_')[0]) for a in df['description'].values]
    zs = [float(a.split('_')[1]) for a in df['description'].values]
    es = df['fa_elec'].values

    ds_, zs_, es_ = [], [], []
    for d, z, e in zip(ds, zs, es):
        if 1: #d % 1 == 0 or d % 1 == 0.5:
            ds_.append(d)
            zs_.append(z)
            es_.append(e)
    return ds_, zs_, es_


def create_pdb_from_d_z(d: float, z:float) -> None:
    z1, z0 = divmod( z, 1 )
    d1, d0 = divmod( d, 1 )
    d0 *= 10
    z0 *= 10

    z1 = int(z1)
    z0 = int(z0)
    d0 = int(d0)
    d1 = int(d1)

    d_zeros = '00' if len(str(d1)) == 1 else '0'
    z_zeros = '00' if len(str(z1)) == 1 else '0'
    atom1 = "HETATM    1 CA    CA A           0.000   0.000   %i.%i%s  1.00 1.00          CA 1" % (z1, z0, z_zeros)
    atom2 = "HETATM    1 CA    CA B           %i.%i%s   0.000   %i.%i%s  1.00 1.00          CA 1" % (d1, d0, d_zeros, z1, z0, z_zeros)
    with open('pdbs/%.2f_%.2f.pdb' % (d, z), 'w+') as fout:
        fout.write(atom1 + '\n')
        fout.write(atom2 + '\n')


def parse_rosetta_log(args):
    args['logger'].log('reading log file, requires lines to have format\nspline z d dielectric')
    args['logger'].log('MAKE SURE TO HAVE THE spline numbers only ONCE in the log file, rosetta prints them multiple times!!!!')
    z = []
    d = []
    e = []

    for l in open(args['log'], 'r'):
        s = l.split()
        if len(s) < 4: continue
        if s[0] == 'spline' or s[0] == 'spl_':
        # if s[0] == 'spl_':
            # if '.' in s[1] or '.' in s[2] and float(s[2]) > 1.4: continue # skip half A z, to quicken things
            if args['mode'] == '3d':
                if '.' in s[2]: continue
                if '.' in s[1]: continue
            if float(s[1]) > 20: continue
            z.append(float(s[1]))
            d.append(float(s[2]))
            e.append(float(s[3]))

    args['logger'].log('found %i lines' % len(z))
    args['logger'].log('z values had max %.2f, min %.2f, mean %.2f and %i entries' % (np.max(z), np.min(z), np.mean(z), len(z)))
    args['logger'].log('d values had max %.2f, min %.2f, mean %.2f and %i entries' % (np.max(d), np.min(d), np.mean(d), len(d)))
    args['logger'].log('e values had max %.2f, min %.2f, mean %.2f and %i entries' % (np.max(e), np.min(e), np.mean(e), len(e)))

    return z, d, e


if __name__ == '__main__':
    main()
