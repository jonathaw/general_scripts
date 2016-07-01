#!/usr/bin/env python3.5
import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', type=str, default='Kd')
    parser.add_argument('-csv', type=str)
    parser.add_argument('-plate_desc', type=str)
    parser.add_argument('-conc', nargs='+', type=float)
    parser.add_argument('-ligand', type=str)
    args = vars(parser.parse_args())

    if args['mode'] == 'Kd':
        plot_and_calc_Kds(args)

    else:
        print('no mode found')


def plot_and_calc_Kds(args):
    df = pd.read_csv(args['csv'])
    print(df)
    plate_desc, plate_pos = parse_plate_description_file(args['plate_desc'])
    print(plate_desc)
    print(plate_pos)
    result_dict = OrderedDict((strain, []) for strain in plate_desc['column_strains'])
    for strain in plate_desc['column_strains']:
        for conc in plate_desc['conc']:
            well = plate_pos[strain][conc]
            print(strain, conc, well, df[df['Well Name'] == well]['P4 APC-A Mean'].values)
            result_dict[strain].append(df[df['Well Name'] == well]['P4 APC-A Mean'].values[0])
    ax = plt.subplot()
    for k, v in result_dict.items():
        print(k, v)
        plt.plot(plate_desc['conc'], v, label=k)
        # ax.set_xscale('log')
        # ax.set_xlim([])
    # ax.set_yscale('log')
        plt.legend()
        plt.show()



def parse_plate_description_file(file_name):
    """

    """
    result = {'missing_cols': [], 'missing_rows': []}
    for l in open(file_name, 'r'):
        s = l.split()
        if s[0] == 'ligand':
            result[s[0]] = s[1]
        elif s[0] == 'column_strains':
            result['column_strains'] = s[1:]
        elif s[0] == 'conc':
            result['conc'] = [int(a) for a in s[1:]]
        elif s[0] == 'missing_cols':
            result['missing_cols'] = [int(a) for a in s[1:]]
        elif s[0] == 'missing_rows':
            result['missing_rows'] = s[1:]
        else:
            print('unrecognised row!!!', l)
            sys.exit()

    pos_dict = {col: {} for col in result['column_strains']}
    for ni, n in enumerate(range(1, 13)):
        if n in result['missing_cols']:
            continue
        for ai, a in enumerate(list('ABCDEFGH')):
            if a in result['missing_rows']:
                continue
            pos_dict[result['column_strains'][ni]][result['conc'][ai]] = '%s%i' % (a, n)
            # pos_dict['%s%i' % (a, n)] = {'conc': result['conc'][ai],
            #                              'column_strain': result['column_strains'][ni]}

    return result, pos_dict

if __name__ == '__main__':
    main()
