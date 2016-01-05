#!/usr/bin/env python3.5

import argparse
import pandas as pd

from seq_funcs import read_multi_fastas

aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def main():
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)
    parser = argparse.ArgumentParser()
    parser.add_argument('-gremlin_file', type=str, help='full path to Gremlin output')
    parser.add_argument('-MSA', type=str, help='full path the MSA file')
    parser.add_argument('-probability_threshold', type=float, help='threshold above which probability is considered')
    parser.add_argument('-query_name', type=str, help='query name as it is written in the MSA')
    args = vars(parser.parse_args())

    gremlin = parse_gremlin(args['gremlin_file'], args['probability_threshold'])
    msa = read_multi_fastas(args['MSA'], add_aligned=True)

    with open(args['query_name']+'.gssm', 'w+') as fout:
        for k, v in gremlin.items():
            iden_frq = create_identitiy_frequency_df(k, msa, args['query_name'])
            fout.write('pos_1 %i pos_2 %i probability %f\n' % (k[0], k[1], v))
            fout.write(str(iden_frq) + '\n')


def create_identitiy_frequency_df(positions: tuple, msa_: dict, query_name: str) -> pd.DataFrame:
    """
    :param positions: tuple of two position ints
    :param msa_: {name: AASeq} of MSA
    :param query_name: query name
    :return: dataframe with i positions in rows, j positions in columns, and df_[r][k] = how many times the tuple k
    (at i) and r (at j) occured, devided by the number of time k occured at i.
    """
    df_ = pd.DataFrame(columns=aas, index=aas)
    i_in_aligned = msa_[query_name].non_aligned_position_at_aligned(positions[0])
    j_in_aligned = msa_[query_name].non_aligned_position_at_aligned(positions[1])
    tup_frq = {(aa1, aa2): 0 for aa1 in aas for aa2 in aas}
    counters = {aa: 0 for aa in aas}
    for seq in msa_.values():
        if len(seq.get_aligned()) >= max(i_in_aligned, j_in_aligned):
            iden_1 = seq.get_aligned_positions([i_in_aligned])[0]
            iden_2 = seq.get_aligned_positions([j_in_aligned])[0]
            if iden_1 != '-' and iden_2 != '-':
                tup_frq[(iden_1, iden_2)] += 1
                counters[iden_1] += 1
    for tup, res in tup_frq.items():
        if counters[tup[0]] == 0:
            df_[tup[1]][tup[0]] = 0.0
        else:
            df_[tup[1]][tup[0]] = float(res)/float(counters[tup[0]])
    return df_


def parse_gremlin(gremlin_file: str, probability_threshold: float) -> dict:
    """
    :param gremlin_file: full path to Gremlin output path where columns are: "i j i_id j_id r_sco s_sco prob"
    :param probability_threshold: probability (prob) over which to take
    :return: {(pos_i, pos_j): probability}
    """
    result = {}
    with open(gremlin_file, 'r') as fin:
        for l in fin:
            s = l.split()
            if s[0] == 'i':
                continue
            if float(s[6]) < probability_threshold:
                break
            result[(int(s[0]), int(s[1]))] = float(s[6])
    return result

if __name__ == '__main__':
    main()
