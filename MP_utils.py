#!/usr/bin/env python3.5
"""
various utility functions for membrane protein related calculations
"""
import numpy as np
import re
import sys
from scipy import interpolate
from CreateEMPBenchmark import POS_Z_DICT_total, Z_total, AAs, get_elazar_scale
from MyPDB import MyPDB


def find_spans_pdb(pdb: MyPDB) -> list:
    memb_vec = []
    for res in pdb.iter_all_res():
        memb_vec.append( res.memb_z is not None )

    spans = []
    chain, start, end, ori = pdb.get_res(1).chain, -1, -1, 'fwd'
    tm_open = False
    for i in range(len(memb_vec)):
        if pdb.get_res(i+1).chain != chain and tm_open:
            end = i
            ori = 'fwd' if pdb.get_res(start).memb_z < 0 and pdb.get_res(end) > 0 else 'rev'
            spans.append([start, end, ori])
            tm_open = False
        if memb_vec[i] and not tm_open:
            start = i + 1
            tm_open = True
        if not memb_vec[i] and tm_open:
            end = i-1
            ori = 'fwd' if pdb.get_res(start).memb_z < 0 and pdb.get_res(end) > 0 else 'rev'
            spans.append([start, end, ori])
            tm_open = False
        chain = pdb.get_res(i+1).chain
    if tm_open:
        spans.append([start, i+1, ori])
    return spans


def draw_splines():
    import matplotlib.pyplot as plt
    fig = plt.figure()
    for i, aa in enumerate(AAs):
        ax = plt.subplot(5, 4, i+1)
        plt.plot(Z_total, interpolate.splev(Z_total, rosetta_splines[aa], der=0), label='spline', color='k')
        plt.vlines(-15, -2, 3, color='grey', linestyles='dashed')
        plt.vlines(15, -2, 3, color='grey', linestyles='dashed')
        plt.title(aa.upper())
    plt.savefig('splines_used.png')
    plt.close()


def calc_rosetta_splines():
    result = {}
    x = Z_total
    for l in open('/home/labs/fleishman/jonathaw/Rosetta/main/database/scoring/score_functions/MembranePotential/elazar_spline_mp_span_ins_fa.txt', 'r'):
        s = l.rstrip().split()
        y = [float(a) for a in s[1:]]
        tck = interpolate.splrep(x, y)
        result[l[0]] = tck
    return result


def calc_seq_dg(seq: str, ori: str, rosetta_splines: dict):
    """
    calculate the dG indsertion as calculated by the splines provided from sequence.
    """
    seq = seq if ori == 'fwd' else seq[::-1]
    membrane_position = np.linspace(-15, 15, endpoint=True, num=len(seq))
    grade = 0
    for i, aa in enumerate(seq):
        g = np.round( interpolate.splev( membrane_position[i], rosetta_splines[aa], der=0 ), 1 )
        grade += g
    return grade


def calc_span_dg(pdb: MyPDB, start: int, end: int, rosetta_splines: dict) -> float:
    """
    calculate a span's dG from structure
    """
    dg = 0.0
    for i in range(start, end+1):
        dg += interpolate.splev( pdb.get_res(i).memb_z, rosetta_splines[ pdb.get_res(i).res_type ], der=0 )
    return dg


def find_topo(ts):
    result = []
    hhh_list = [[a.span()[0], a.span()[1]-1] for a in re.finditer('[Hh]*', ts) if a.span()[1]-a.span()[0] > 1]
    ori = ''
    for h in hhh_list:
        if ts[h[0] - 1] == '1' or ts[h[1] + 1 ] == '2':
            ori = 'fwd'
        elif ts[h[0] - 1] == '2' or ts[h[1] +  1] == '1':
            ori = 'rev'
        else:
            print(ts)
            print('dont know what to do')
            sys.exit()
        result.append([ h[0], h[1], ori ])
        # print(ts[h[0]-1 : h[1]+2], ori)
    return result


def parse_rost_db():
    fin = open('./rostlab_db.txt', 'r')
    cont = fin.read().split('>')
    result = {}
    for par in cont[1:]:
        if 'PDBTM' in par:
            d = par.split('\n')
            name = d[0].split('|')[0].lower()
            result[name] = {'name': name, 'seq': d[1], 'ts': d[2]}
    return result


if __name__ == '__main__':
    main()
