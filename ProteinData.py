#!/usr/bin/env python3.5
"""
script bundle to calculate protein conc. etc.
"""
import argparse
from collections import OrderedDict
import pandas as pd

from AASeq import AASeq


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-in_file', type=str)

    args = vars(parser.parse_args())
    input_dict = parse_input_data(args['in_file'])

    pd.set_option('display.float_format', '{:.2g}'.format)
    df = pd.DataFrame(columns=['name', 'seq', 'dilution_factor', 'absorbance', 'molecular_weight', 'pI',
                               'extinction_coefficient'])

    for k, v in input_dict.items():
        # calculate extinction coefficient
        v['extinction_coefficient'] = v['seq'].calc_extinction_coefficient(reduced=False)

        # calculate Isoelectroc point
        v['pI'] = v['seq'].calc_isoelectric_point()

        # calculate molar concentration
        v['conc'] = v['dilution_factor'] * v['absorbance'] / v['extinction_coefficient']

        # calcualte concentration if dilued by half
        v['glycerol_conc'] = v['conc'] / 2

        # calculate molecular weight
        v['molecular_weight'] = v['seq'].calc_molecular_weight()

        # calculate g/L
        v['g/l'] = v['conc'] / v['molecular_weight']

        print_evernote_format(v)

        v['seq'] = v['seq'].get_seq()
        df = df.append(v, ignore_index=True)

    print(df)


def print_evernote_format(entry: dict) -> None:
    print('%s MW %i, %iAA, £%i, OD %.3f %.2eµM' % (entry['name'], entry['molecular_weight'], len(entry['seq']),
                                                entry['extinction_coefficient'], entry['absorbance'],
                                                entry['glycerol_conc']))


def parse_input_data(in_file: str) -> OrderedDict:
    """
    :param in_file: input table. use the template
    :return: dict of the CSV
    """
    with open(in_file, 'r') as fin:
        cont = fin.read().split('\n')
    result = OrderedDict({})
    for l in cont:
        s = l.split(',')
        if s[0] == 'name':
            continue
        result[s[0]] = {'name': s[0],
                        'seq': AASeq(s[1], name=s[0]),
                        'dilution_factor': float(s[2]),
                        'absorbance': float(s[3])}
    return result


if __name__ == '__main__':
    main()
