#!/usr/bin/env python3.5
"""
iterate over all entries in the Rost dataset, and claculate the dG of insertion
according to the splines in rosetta (in kcal/mol) that describe dG,
not the calibration of the score funciton to dsTbL.
returns the mean, std and median for the dG. also draws a histogram of the distribution.
"""
import numpy as np
import re
import sys
from scipy import interpolate
from CreateEMPBenchmark import POS_Z_DICT_total, Z_total, AAs, get_elazar_scale

def main():
    global rosetta_splines
    rosetta_splines = calc_rosetta_splines()
    draw_splines()
    rost_db = parse_rost_db()
    all_dgs = []

    if True:
        e = rost_db['p11551']
        t = find_topo(e['ts'])
        for h in t:
            print(h, e['seq'][h[0]: h[1]+1], calc_seq_dg(e['seq'][h[0]: h[1]+1], h[2]), rosetta_splines)
        sys.exit()
        
    for entry in rost_db.values():
        topo = find_topo(entry['ts'])
        for h in topo:
            dg = calc_seq_dg(entry['seq'][h[0]: h[1]+1], h[2], rosetta_splines)
            all_dgs.append( dg )
    plt.hist(all_dgs, bins=50)
    print('the mean is %.2f' % np.mean(all_dgs))
    print('the median is %.2f' % np.median(all_dgs))
    print('the std %.2f' % np.std(all_dgs))
    plt.savefig(natural_TMs_span_dG.png)


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


def calc_seq_dg(seq, ori, rosetta_splines):
    seq = seq if ori == 'fwd' else seq[::-1]
    membrane_position = np.linspace(-15, 15, endpoint=True, num=len(seq))
    grade = 0
    for i, aa in enumerate(seq):
        g = np.round( interpolate.splev( membrane_position[i], rosetta_splines[aa], der=0 ), 1 )
        grade += g
    return grade


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
