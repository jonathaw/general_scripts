#!/usr/bin/env python3.5
"""
go over an MSA, and read al given positions (designed positions). display the charge distribution, Pos, Neg and Tot
"""
__author__ = 'jonathan'
from seq_funcs import read_multi_fastas
from AASeq import AASeq
import matplotlib.pyplot as plt
res2charge = {'K': 'p', 'R': 'p', 'D': 'n', 'E': 'n'}


def extract_charge_configuration(seq: AASeq, positions: list):
    res_in_poses = seq.get_positions(positions)
    charge = [res2charge[a] if a in res2charge.keys() else 'c' for a in res_in_poses]
    return charge


if __name__ == '__main__':
    positions = [32, 33, 35, 37, 39, 63, 66, 68, 70, 73, 75, 77, 79, 81, 82, 83, 85, 87, 115, 116, 118, 119, 121, 123,
                 125, 127]
    fastas = read_multi_fastas('/home/labs/fleishman/jonathaw/data/pssm/cohs/making/1ohz_passed_thresholds.fasta')
    positives, negatives, neutrals, totals = [], [], [], []
    for k, v in fastas.items():
        charge_config = extract_charge_configuration(v, positions)
        positives.append(charge_config.count('p'))
        negatives.append(charge_config.count('n'))
        neutrals.append(charge_config.count('c'))
        totals.append(charge_config.count('p')+charge_config.count('n'))
        if k == '1ohz':
            print('at 1ohz found %i %s' % (charge_config.count('p'), 'positives'))
            print('at 1ohz found %i %s' % (charge_config.count('n'), 'negatives'))
            print('at 1ohz found %i %s' % (charge_config.count('c'), 'neutrals'))
            print('at 1ohz found %i %s' % (charge_config.count('p')+charge_config.count('n'), 'totals'))
    bins = range(max(positives+negatives+neutrals)+1)
    plt.hist(positives, bins=bins, color='b', label='positives')
    plt.hist(negatives, bins=bins, color='r', label='negatives')
    plt.hist(neutrals, bins=bins, color='k', label='neutrals')
    plt.hist(totals, bins=bins, color='g', label='total')
    print('positives', positives)
    print('negatives', negatives)
    print('neutrals', neutrals)
    print('totals', totals)
    plt.legend()
    plt.show()