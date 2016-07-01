#!/usr/bin/env python3.5
from __future__ import print_function

from seq_funcs import read_multi_fastas
from AASeq import AASeq
from _binding_data import binding_data

import pandas as pd
from collections import OrderedDict
import matplotlib.pyplot as plt
from matplotlib import colors
from numpy import array, arange
# import sys
# import os
# import subprocess

# import numpy as np
from sklearn.tree import DecisionTreeClassifier, export_graphviz

root_path = '/home/labs/fleishman/jonathaw/decision_tree/'
aa2num = OrderedDict([('-', 0), ('A', 1), ('C', 2), ('D', 3), ('E', 4), ('F', 5), ('G', 6), ('H', 7), ('I', 8),
                       ('K', 9), ('L', 10), ('M', 11), ('N', 12), ('P', 13), ('Q', 14), ('R', 15), ('S', 16), ('T', 17),
                       ('V', 18), ('W', 19), ('Y', 20)])
num2aa = OrderedDict({v: k for k, v in aa2num.items()})
types = ['p', 'n', 'o', 'h']
positions = {'core_coh': [11, 12], 'core_doc': [2, 6, 10]}
positions['rim_coh'] = [i for i in range(1, 21, 1) if i not in positions['core_coh']]
positions['rim_doc'] = [i for i in range(1, 11, 1) if i not in positions['core_doc']]

coh_poses_1ohz = [32, 33, 35, 37, 63, 66, 68, 70, 73, 75, 77, 79, 81, 83, 85, 116, 118, 119, 121, 123, 125, 127]
doc_poses_1ohz = [11, 14, 15, 18, 19, 21, 22, 45, 46, 48, 49, 52, 53]


def main():
    print('parsing data')
    df = parse_binding_data()
    # print(df)

    features = list(df.columns[2:-1])

    y = df['binding']
    X = df[features]

    print('fitting tree')
    dt = DecisionTreeClassifier()#min_samples_split=20)#, random_state=99)
    dt.fit(X, y)

    df['prediction'] = dt.predict(X)
    predict = {coh: {doc: int(df.loc[(df.coh == coh) & (df.doc == doc)].prediction) for doc in df['doc']}
               for coh in set(df['coh'])}

    # with open('decision_tree.dot', 'w') as fout:
    #     export_graphviz(dt, out_file=fout, feature_names=features)

    show_predcition_matrix(predict)
    # command = ['dot', '-Tpng', 'decision_tree.dot', '-o', 'decision_tree.png']
    # subprocess.check_call(command)


def show_predcition_matrix(prediction: dict) -> None:
    obs_pre = {0: {0: 0, 1: 2}, 1: {0: 3, 1: 1}}
    binding_ = binding_data()
    df = pd.DataFrame(data=0, index=list(binding_.values())[0].keys(), columns=binding_.keys())
    for coh, doc_ in binding_.items():
        for doc, obs in doc_.items():
            df[coh][doc] = obs_pre[binding_[coh][doc]][prediction[coh][doc]]

    plt.figure()
    axis = plt.gca()
    cmap = colors.ListedColormap(['white', 'cornflowerblue', 'red', 'darkorange'])
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    heatmap = plt.pcolor(array(df), cmap=cmap, norm=norm, edgecolors='k', linewidth=2)
    for y in range(array(df.shape)[0]):
        for x in range(array(df.shape)[1]):
            if array(df)[y, x] >= 0:
                plt.text(x+0.5, y+0.5, array(df)[y, x], horizontalalignment='center', verticalalignment='center')
    plt.yticks(arange(0.5, len(df.index), 1), df.index)
    plt.xticks(arange(0.5, len(df.columns), 1), df.columns, rotation=70)
    plt.xlabel('Cohesin name', style='oblique')
    plt.ylabel('Dockerin name', style='oblique')
    axis.set_aspect('equal')
    plt.title('Cohesin dockerin cross binding')
    plt.suptitle('0: obs no pred no, 1: obs yes, pred yes\n2: obs no pred yes, 3: obs yes pred no')
    plt.show()


def seq_to_num(seq: str) -> int:
    return int(''.join(str(aa2num[a]) for a in seq))


def parse_binding_data() -> pd.DataFrame:
    # cohs = read_multi_fastas(root_path+'cohs_specific_pos.fasta', suffix_to_remove='/')
    # docs = read_multi_fastas(root_path+'docs_specific_pos.fasta', suffix_to_remove='/')
    cohs, docs = retrive_relevant_poses()
    data = binding_data()

    colums = ['coh', 'doc'] + \
             ['core coh %i %s' % (i, aa) for i in [1, 2] for aa in aa2num.keys()] + \
             ['core doc %i %s' % (i, aa) for i in [1, 2, 3] for aa in aa2num.keys()] + \
             ['rim coh %i %s' % (i, t) for i in range(1, 19, 1) for t in types] + \
             ['rim doc %i %s' % (i, t) for i in range(1, 8, 1) for t in types] + ['binding']
    df = pd.DataFrame(columns=colums)
    i = 1
    for coh, doc_dict in data.items():
        coh_seq = cohs[coh].get_seq
        for doc, res in doc_dict.items():
            doc_seq = docs[doc].get_seq
            df.loc[i] = [coh, doc] + seqs2row(coh_seq, doc_seq) + [1 if res else 0]
            i += 1
    return df


def retrive_relevant_poses() -> (dict, dict):
    """
    :return: seq dicts for cohs and docs, holding only the relevqant positions, determined by 1OHZ
    """
    cohs_old = read_multi_fastas(root_path+'cohesins_from_rachel.fasta_aln', suffix_to_remove='/')
    docs_old = read_multi_fastas(root_path+'dockerins_from_rachel.fasta_aln', suffix_to_remove='/')

    coh_1ohz = cohs_old['1OHZ']
    coh_poses = [coh_1ohz.non_aligned_position_at_aligned(p) for p in coh_poses_1ohz]
    doc_1ohz = docs_old['1OHZ']
    doc_poses = [doc_1ohz.non_aligned_position_at_aligned(p) for p in doc_poses_1ohz]

    cohs_new, docs_new = {}, {}

    for coh, res in cohs_old.items():
        cohs_new[coh] = AASeq(string=''.join(res.get_aligned_positions(coh_poses)), name=coh)
    for doc, res in docs_old.items():
        docs_new[doc] = AASeq(string=''.join(res.get_aligned_positions(doc_poses)), name=doc)
    return cohs_new, docs_new


def aa2rim_type(aa: str) -> str:
    if aa in ['K', 'R', 'H']:
        return 'p'
    elif aa in ['D', 'E']:
        return 'n'
    elif aa in ['N', 'Q', 'T', 'S', 'C']:
        return 'o'
    elif aa in ['L', 'I', 'M', 'V', 'A', 'F', 'W', 'Y', 'P', 'G']:
        return 'h'


def aa2binary_core(aa: str) -> list:
    """
    :param aa: amino acid
    :return: binary for aa
    >>> aa2binary('A')
    [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    >>> aa2binary('W')
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
    """
    res = [0] * len(aa2num.keys())
    res[aa2num[aa]] = 1
    return res


def aa2binary_rim(aa: str) -> list:
    """
    :param aa: an amino acid
    :return: binary for rim type
    >>> aa2binary_rim('A')
    [0, 0, 0, 1]
    >>> aa2binary_rim('E')
    [0, 1, 0, 0]
    """
    t = aa2rim_type(aa)
    res = [0] * len(types)
    for i, ty in enumerate(types):
        if t == ty:
            res[i] = 1
    return res


def seqs2row(coh_seq: str, doc_seq: str) -> list:
    """
    :param coh_seq: cohesin sequence
    :param doc_seq: dockerin sequence
    :return: binary list repersenting the sequences
    >>> seqs2row('SKTRVSSNDYTGLETALGL-', 'SSLSTSNSLA')

    """
    res = []
    for p in positions['core_coh']:
        res += aa2binary_core(coh_seq[p-1])
    for p in positions['core_doc']:
        res += aa2binary_core(doc_seq[p-1])

    for p in positions['rim_coh']:
        res += aa2binary_rim(coh_seq[p-1])
    for p in positions['rim_doc']:
        res += aa2binary_rim(doc_seq[p-1])
    return res


if __name__ == '__main__':
    main()