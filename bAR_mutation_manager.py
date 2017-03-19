#!/usr/bin/env python3.5
"""
"""
import pandas as pd
import argparse
from RosettaFilter import score_file2df
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode')
    parser.add_argument('-sc')
    parser.add_argument('-sclog')
    args = vars(parser.parse_args())

    if args['mode'] == 'bad_mut_table':
        mut_df = parse_mutations_table()
        print_res_file(mut_df)

    elif args['mode'] == 'parse_bar_plot_based_table':
        muts = parse_bar_plot_based_table()
        print_res_file(muts)

    elif args['mode'] == 'bar_compare':
        bar_compare(args)

    else:
        print('no mode')


def bar_compare(args: dict) -> None:
    exp_dict = parse_bar_plot_based_table()
    sclog_df = pd.read_csv(args['sclog'], sep='\s+', header=None,
                           names=['pos', 'wt_aa', 'mut_aa', 'total_score'])
    sc_df = score_file2df(args['sc'])
    wt_score = sc_df.iloc[[-1]]['total_score'].values[-1]
    print('wt_score', wt_score)
    print('VVVV')
    print(sclog_df)
    print('BBBB')

    sclog_df['delta'] = sclog_df['total_score'] - wt_score
    sclog_df['exp'] = 1000.0
    sclog_df['original_pos'] = 0
    for index, row in sclog_df.iterrows():
        if row['pos']+35 in exp_dict.keys():
            sclog_df.set_value(index, 'exp',
                               exp_dict[row['pos']+35][row['mut_aa']])
            sclog_df.set_value(index, 'original_pos', row['pos']+35)
        elif row['pos'] + 69 in exp_dict.keys():
            sclog_df.set_value(index, 'exp',
                               exp_dict[row['pos']+69][row['mut_aa']])
            sclog_df.set_value(index, 'original_pos', row['pos']+69)
        else:
            print('cant find %i', row['pos'])
    print(sclog_df.to_string())

    exposed_df = sclog_df[sclog_df['original_pos'].isin([160, 151, 89, 68, 338,
                                                         334, 327, 194, 221,
                                                         229, 160])]

    fig, ax = plt.subplots()
    # ax.scatter(sclog_df['exp'], sclog_df['delta'])
    ax.scatter(exposed_df['exp'], exposed_df['delta'])
    for index, row in exposed_df.iterrows():
        ax.annotate(index, (row['exp'], row['delta']))
    plt.axvline(31.7)
    plt.axhline(0)
    plt.xlabel('experimental Tm')
    plt.ylabel('rosetta ∆G')
    print('select muts, resi %s' %
          '+'.join(str(a) for a in set(list(sclog_df['original_pos']))))
    ur = sclog_df[(sclog_df['exp'] > 31.7) & (sclog_df['delta'] > 0)]
    lr = sclog_df[(sclog_df['exp'] > 31.7) & (sclog_df['delta'] < 0)]
    ll = sclog_df[(sclog_df['exp'] < 31.7) & (sclog_df['delta'] < 0)]
    ul = sclog_df[(sclog_df['exp'] < 31.7) & (sclog_df['delta'] > 0)]
    print('select UR, resi %s'
          % '+'.join(str(a) for a in set(list(ur['original_pos']))))
    print('select LR, resi %s'
          % '+'.join(str(a) for a in set(list(lr['original_pos']))))
    print('select LL, resi %s'
          % '+'.join(str(a) for a in set(list(ll['original_pos']))))
    print('select UL, resi %s'
          % '+'.join(str(a) for a in set(list(ul['original_pos']))))

    plt.show()


def parse_filter_scan_sclog(file_name: str) -> dict:
    result = {}
    with open(file_name, 'r') as fin:
        for l in fin:
            s = l.split()
            if len(s) == 4:
                pos = int(s[0])
                if pos not in result.keys():
                    result[pos] = {'wt': s[1]}
                result[pos][s[2]] = float(s[3])
    return result


def parse_bar_plot_based_table() -> dict:
    with open('/home/labs/fleishman/jonathaw/elazaridis/mut_recap_20Feb/'
              'bAR_20Feb/data/experimental_data/bar_plot_table.tsv') as fin:
        result = {}
        for l in fin:
            if 'pos' in l:
                continue
            s = l.split()
            p = int(s[0])
            if p not in result.keys():
                result[p] = {'wt': s[1]}
            # result[ p ][ s[ 2 ] ] = int(s[ 3 ])
            result[p][s[2]] = 31.7 * (1 + (float(s[3]) - 50.0) / 50.0)
        return result


def print_res_file(d: dict) -> None:
    print('\nnataa\nstart')
    for k, v in d.items():
        # print(k, v)
        print('%i\tA\tPIKAA\t%s' % (k, ''.join(a for a in v if a != 'wt')))


def parse_mutations_table() -> dict:
    with open('/home/labs/fleishman/jonathaw/elazaridis/mut_recap_20Feb/'
              'bAR_20Feb/data/experimental_data/processed_table.tsv',
              'r') as fin:
        result = {}
        for l in fin:
            s = l.split()
            if 'pos' in l:
                continue
            if int(s[0]) not in result.keys():
                result[int(s[0])] = {'wt': s[1]}
            result[int(s[0])][s[2]] = float(s[-1])
    return result

if __name__ == '__main__':
    main()
