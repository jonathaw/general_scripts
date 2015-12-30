#!/usr/bin/env python3.5
import pandas as pd
import numpy as np
import sys
from sklearn.tree import DecisionTreeClassifier, export_graphviz
from _binding_data import binding_data

import matplotlib.pyplot as plt
from matplotlib import colors
from numpy import array, arange

from seq_funcs import read_multi_fastas
from AASeq import AASeq

# positions = {'core_coh': [10, 11], 'core_doc': [1, 5, 9]}
# positions['rim_coh'] = [i for i in range(0, 22, 1) if i not in positions['core_coh']]
# positions['rim_doc'] = [i for i in range(0, 13, 1) if i not in positions['core_doc']]

# coh_poses_1ohz = [32, 33, 35, 37, 63, 66, 68, 70, 73, 75, 77, 79, 81, 83, 85, 116, 118, 119, 121, 123, 125, 127]
# doc_poses_1ohz = [11, 14, 15, 18, 19, 21, 22, 45, 46, 48, 49, 52, 53]

aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
aas_dict = {aa: i+1 for i, aa in enumerate(aas)}
coh_core_aas = ['A', 'G', 'I', 'L', 'S', 'T', 'V']
doc_core_aas = ['A', 'C', 'F', 'I', 'L', 'V', 'Y']

types = ['p', 'n', 'o', 'h']
types_dict = {t: i+1 for i, t in enumerate(types)}
rim_types_to_binary = {'p': [1, 0, 0, 0], 'n': [0, 1, 0, 0], 'o': [0, 0, 1, 0], 'h': [0, 0, 0, 1], '-': [0, 0, 0, 0]}
type_to_res = {'p': ['K', 'R', 'H'], 'n': ['D', 'E'], 'o': ['N', 'Q', 'T', 'S', 'C'],
               'h': ['L', 'I', 'M', 'V', 'A', 'F', 'W', 'Y', 'P', 'G'], 'NA': ['-']}

doc_len_reduction = {'2b59': 96, '2ozn': 81}

decision_tree_root = '/home/labs/fleishman/jonathaw/decision_tree/'

ordered_positions = {'coh': ['coh_core_1', 'coh_core_2', 'coh_rim_1', 'coh_rim_2', 'coh_rim_3', 'coh_rim_4',
                             'coh_rim_5', 'coh_rim_6', 'coh_rim_7', 'coh_rim_8', 'coh_rim_9', 'coh_rim_10',
                             'coh_rim_11', 'coh_rim_12', 'coh_rim_13', 'coh_rim_14', 'coh_rim_15', 'coh_rim_16',
                             'coh_rim_17', 'coh_rim_18', 'coh_rim_19', 'coh_rim_20'],
                     'doc': ['doc_core_1', 'doc_core_2', 'doc_core_3', 'doc_rim_1', 'doc_rim_2', 'doc_rim_3',
                             'doc_rim_4', 'doc_rim_5', 'doc_rim_6', 'doc_rim_7', 'doc_rim_8', 'doc_rim_9', 'doc_rim_10']
                     }


def main():
    data_df = parse_binding_data()
    prepared_df, identities_df = prepare_data(data_df)
    # validate_data_frame(data_df, prepared_df)
    decision_tree, features = create_decision_tree(prepared_df)

    # print(identities_df)
    analyse_identity_df(identities_df)
    sys.exit()

    with open('decision_tree.dot', 'w') as fout:
        print('creating decision tree.dot')
        export_graphviz(decision_tree, out_file=fout, feature_names=features)

    compare_observed_to_predicted(decision_tree, data_df, prepared_df[features])


def create_decision_tree(df: pd.DataFrame) -> DecisionTreeClassifier:
    features = list(df.columns[4:-1])

    X = df[features]
    y = df['binders']

    print('fitting tree')
    dt = DecisionTreeClassifier()#min_samples_split=5)#, random_state=99)
    dt.fit(X, y)
    return dt, features


def compare_observed_to_predicted(dt: DecisionTreeClassifier, data_df: pd.DataFrame, X) -> None:
    """
    :param dt: decision tree
    :param data_df: binding data
    :param X: binary data frame slice
    :return: creates the coh->doc->obse_pred dataframe
    """
    obs_pre = {False: {0: 0, 1: 2}, True: {0: 3, 1: 1}}
    data_df['prediction'] = dt.predict(X)

    predict = data_df[[0, 1, -2, -1]]
    df = pd.DataFrame(index=set(data_df.doc_name), columns=set(data_df.coh_name))
    total = 0
    for i in range(1, len(predict.index)+1):
        row = predict.loc[i]
        df[row.coh_name][row.doc_name] = obs_pre[row.binders][row.prediction]
        if row.binders in [True, False]:
            total += 1
    print('found a total of %i entries with known binding' % total)
    show_predcition_matrix(df)


def show_predcition_matrix(prediction: pd.DataFrame) -> None:
    prediction = prediction.sort_index()
    prediction = prediction.reindex_axis(sorted(prediction.columns), axis=1)
    obs_pre = {0: {0: 0, 1: 2}, 1: {0: 3, 1: 1}}
    plt.figure()
    axis = plt.gca()
    cmap = colors.ListedColormap(['white', 'cornflowerblue', 'red', 'darkorange'])
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    heatmap = plt.pcolor(array(prediction), cmap=cmap, norm=norm, edgecolors='k', linewidth=2)
    for y in range(array(prediction.shape)[0]):
        for x in range(array(prediction.shape)[1]):
            if array(prediction)[y, x] == np.nan:
                continue
            if array(prediction)[y, x] >= 0:
                plt.text(x+0.5, y+0.5, array(prediction)[y, x], horizontalalignment='center', verticalalignment='center')
    plt.yticks(arange(0.5, len(prediction.index), 1), prediction.index)
    plt.xticks(arange(0.5, len(prediction.columns), 1), prediction.columns, rotation=70)
    plt.xlabel('Cohesin name', style='oblique')
    plt.ylabel('Dockerin name', style='oblique')
    axis.set_aspect('equal')
    plt.title('Cohesin dockerin cross binding')
    plt.suptitle('0: obs no pred no, 1: obs yes, pred yes\n2: obs no pred yes, 3: obs yes pred no')
    plt.show()


def validate_data_frame(data_df: pd.DataFrame, prepared_df: pd.DataFrame) -> None:
    """
    :param data_df: binding data frame
    :param prepared_df: binary data frame
    :return: prints if there is something wrong...
    """
    rachel_root = '/home/labs/fleishman/jonathaw/decision_tree/'
    cohs = read_multi_fastas(rachel_root+'cohesins_from_rachel_and_vered.fasta_aln', suffix_to_remove='/')
    docs = read_multi_fastas(rachel_root+'dockerins_from_rachel_and_vered.fasta_aln', suffix_to_remove='/')
    # coh_1ohz = cohs['1OHZ']
    # doc_1ohz = docs['1OHZ']
    coh_crys_seqs = [c for c in cohs.values() if c.name in ['1ohz', '2b59', '2ozn', '2vn5', '2y3n', '3ul4',
                                                                    '4fl4', '4fl5', '4dh2', '4uyp', '5new']]
    doc_crys_seqs = [d for d in docs.values() if d.name in ['1ohz', '2b59', '2ozn', '2vn5', '2y3n', '3ul4',
                                                                    '4fl4', '4fl5', '4dh2', '4uyp', '5new']]
    # coh_poses = [coh_1ohz.non_aligned_position_at_aligned(p) for p in coh_poses_1ohz]
    # doc_poses = [doc_1ohz.non_aligned_position_at_aligned(p) for p in doc_poses_1ohz]
    features = list(prepared_df.columns[4:-1])

    interface_positions = parse_interface_positions()
    coh_poses = {coh: {typ: cohs[coh].non_aligned_position_at_aligned(pos) for typ, pos in typos.items()} for coh, typos
                 in interface_positions['coh'].items()}
    doc_poses = {doc: {typ: docs[doc].non_aligned_position_at_aligned(pos) for typ, pos in typos.items()} for doc, typos
                 in interface_positions['doc'].items()}

    for i in range(1, len(data_df.index)):
        # i = len(data_df.index)
        print('i is %i' % i)
        print(data_df.loc[i])
        if data_df.loc[i]['coh_name'] != prepared_df.loc[i]['coh_name'] or \
                        data_df.loc[i]['doc_name'] != prepared_df.loc[i]['doc_name']:
            print('not the same names', data_df.loc[i]['doc_name'], prepared_df.loc[i]['doc_name'])
            sys.exit()
        coh_seq = data_df.loc[i]['coh_seq']
        # coh_q_poses = coh_seq.get_aligned_positions(coh_poses)

        doc_seq = data_df.loc[i]['doc_seq']
        # doc_q_poses = doc_seq.get_aligned_positions(doc_poses)

        prepared_row = row_to_dict(prepared_df.loc[i])

        similar_coh, coh_iden = highest_seq_similarity(coh_crys_seqs, data_df.loc[i]['coh_seq'])
        similar_doc, doc_iden = highest_seq_similarity(doc_crys_seqs, data_df.loc[i]['doc_seq'])
        coh_identities = {typ: data_df.loc[i]['coh_seq'].get_aligned_positions([pos])[0] for typ, pos in
                          coh_poses[similar_coh.name].items()}
        doc_identities = {typ: data_df.loc[i]['doc_seq'].get_aligned_positions([pos])[0] for typ, pos in
                          doc_poses[similar_doc.name].items()}



        # for pos in positions['core_coh']:
        #     if coh_q_poses[pos] != prepared_row['coh_core_%i' % pos]:
        #         print('not the same coh query pos differs from row', pos, coh_q_poses[pos], prepared_row['coh_core_%i' % pos])
        #         sys.exit()

        # for pos in positions['core_doc']:
        #     if doc_q_poses[pos] != prepared_row['doc_core_%i' % pos]:
        #         print('not the same doc query pos differs from row')
        #         sys.exit()

        # for pos in positions['rim_coh']:
        #     if [k for k, v in type_to_res.items() if coh_q_poses[pos] in v][0] != prepared_row['coh_rim_%i' % pos] and \
        #             not ([k for k, v in type_to_res.items() if coh_q_poses[pos] in v][0] == 'NA' and
        #                          prepared_row['coh_rim_%i' % pos] == '-'):
        #         print('breaking', [k for k, v in type_to_res.items() if coh_q_poses[pos] in v][0],
        #               prepared_row['coh_rim_%i' % pos])
        #         sys.exit()

        for fea in features:
            if prepared_df.loc[i][fea] not in [0, 1]:
                print('found problem at row', i, prepared_df.loc[i][fea])

        # break

    print('your df is validated')


def row_to_dict(row) -> dict:
    """
    :param row: row from binary data frame
    :return: dictionary portrayng the row
    """
    result = {'coh_name': row['coh_name'], 'doc_name': row['doc_name'], 'coh_seq': row['coh_seq'],
              'doc_seq': row['doc_seq']}
    for pos in positions['core_coh']:
        result['coh_core_%i' % pos] = binary_to_res([row['coh_core_%i_%s' % (pos, aa)] for aa in aas])
    for pos in positions['core_doc']:
        result['doc_core_%i' % pos] = binary_to_res([row['doc_core_%i_%s' % (pos, aa)] for aa in aas])

    for pos in positions['rim_coh']:
        result['coh_rim_%i' % pos] = [k for k, v in rim_types_to_binary.items() if v == [row['coh_rim_%i_%s' % (pos, t)]
                                                                                         for t in types]][0]
    for pos in positions['rim_doc']:
        result['doc_rim_%i' % pos] = [k for k, v in rim_types_to_binary.items() if v == [row['doc_rim_%i_%s' % (pos, t)]
                                                                                         for t in types]][0]
    return result


def prepare_data(in_df: pd.DataFrame) -> (pd.DataFrame, pd.DataFrame):
    """
    :rtype: (pd.DataFrame, pd.DataFrame)
    """
    rachel_root = '/home/labs/fleishman/jonathaw/decision_tree/'
    cohs_non_aln = read_multi_fastas(rachel_root+'cohesins_from_rachel_and_vered.fasta', suffix_to_remove='/', lower=True)
    docs_non_aln = read_multi_fastas(rachel_root+'dockerins_from_rachel_and_vered.fasta', suffix_to_remove='/', lower=True)

    cohs = read_multi_fastas(rachel_root+'cohesins_from_rachel_and_vered.fasta_aln', suffix_to_remove='/', lower=True)
    docs = read_multi_fastas(rachel_root+'dockerins_from_rachel_and_vered.fasta_aln', suffix_to_remove='/', lower=True)
    interface_positions = parse_interface_positions()
    coh_poses = {coh: {typ: cohs[coh].non_aligned_position_at_aligned(pos) for typ, pos in typos.items()} for coh, typos
                 in interface_positions['coh'].items()}
    doc_poses = {doc: {typ: docs[doc].non_aligned_position_at_aligned(pos) for typ, pos in typos.items()} for doc, typos
                 in interface_positions['doc'].items()}

    validate_aligned_non_aligned_interface_positions(interface_positions['coh'], cohs, cohs_non_aln)
    validate_aligned_non_aligned_interface_positions(interface_positions['doc'], docs, docs_non_aln)

    # coh_1ohz = cohs['1OHZ']
    # doc_1ohz = docs['1OHZ']
    # coh_poses = [coh_1ohz.non_aligned_position_at_aligned(p) for p in coh_poses_1ohz]
    # doc_poses = [doc_1ohz.non_aligned_position_at_aligned(p) for p in doc_poses_1ohz]
    coh_crys_seqs = [c for c in cohs.values() if c.name in ['1ohz', '2b59', '2ozn', '2vn5', '2y3n', '3ul4',
                                                                    '4fl4', '4fl5', '4dh2', '4uyp', '5new']]
    doc_crys_seqs = [d for d in docs.values() if d.name in ['1ohz', '2b59', '2ozn', '2vn5', '2y3n', '3ul4',
                                                                    '4fl4', '4fl5', '4dh2', '4uyp', '5new']]

    # columns = ['coh_name', 'doc_name', 'coh_seq', 'doc_seq'] + \
    #           ['coh_core_%i_%s' % (v, aa) for k, v in ordered_positions.items() if 'core' in k for aa in aas] + \
    #           ['coh_rim_%i_%s' % (v, t) for k, v in ordered_positions.items() if 'rim' in k for t in types] + \
    #           ['doc_core_%i_%s' % (v, aa) for k, v in ordered_positions.items() if 'core' in k for aa in aas] + \
    #           ['doc_rim_%i_%s' % (v, t) for k, v in ordered_positions.items() if 'rim' in k for t in types] + \
    #           ['binders']

    columns = ['coh_name', 'doc_name', 'coh_seq', 'doc_seq'] + \
              ['%s_%s' % (typ, aa) for typ in ordered_positions['coh'] if 'core' in typ for aa in aas] + \
              ['%s_%s' % (typ, aa) for typ in ordered_positions['coh'] if 'rim' in typ for aa in types] + \
              ['%s_%s' % (typ, aa) for typ in ordered_positions['doc'] if 'core' in typ for aa in aas] + \
              ['%s_%s' % (typ, aa) for typ in ordered_positions['doc'] if 'rim' in typ for aa in types] + ['binders']
    out_df = pd.DataFrame(index=range(1, len(in_df.index)), columns=columns)

    id_columns = ['coh_name', 'doc_name', 'coh_seq', 'doc_seq'] + ordered_positions['coh'] + ordered_positions['doc']
    identities_df = pd.DataFrame(index=range(1, len(in_df.index)), columns=id_columns)

    for i in range(1, len(in_df.index)+1):
        # find which crystal coh+doc are most similar
        # get aligned positions accrotding to interface_positions
        similar_coh, coh_iden = highest_seq_similarity(coh_crys_seqs, in_df.loc[i]['coh_seq'])
        similar_doc, doc_iden = highest_seq_similarity(doc_crys_seqs, in_df.loc[i]['doc_seq'])
        coh_identities = {typ: in_df.loc[i]['coh_seq'].get_aligned_positions([pos])[0] for typ, pos in coh_poses[similar_coh.name].items()}
        doc_identities = {typ: in_df.loc[i]['doc_seq'].get_aligned_positions([pos])[0] for typ, pos in doc_poses[similar_doc.name].items()}
        # coh_ = in_df.loc[i]['coh_seq'].get_aligned_positions(coh_poses[similar_coh])
        # doc_ = in_df.loc[i]['doc_seq'].get_aligned_positions(doc_poses[similar_doc])

        coh_core = [core_res_to_identity(coh_identities[v]) for v in ordered_positions['coh'] if 'core' in v]
        coh_rim = [rim_res_to_type_binary(coh_identities[v]) for v in ordered_positions['coh'] if 'rim' in v]
        doc_core = [core_res_to_identity(doc_identities[v]) for v in ordered_positions['doc'] if 'core' in v]
        doc_rim = [rim_res_to_type_binary(doc_identities[v]) for v in ordered_positions['doc'] if 'rim' in v]
        # coh_core, coh_rim = seq_to_binary(coh_, 'coh')
        # doc_core, doc_rim = seq_to_binary(doc_, 'doc')

        # print(coh_core)
        # print(coh_rim)
        # print(doc_core)
        # print(doc_rim)

        coh_core_list, coh_rim_list = [], []
        [coh_core_list.append(a) for b in coh_core for a in b]
        [coh_rim_list.append(a) for b in coh_rim for a in b]

        doc_core_list, doc_rim_list = [], []
        [doc_core_list.append(a) for b in doc_core for a in b]
        [doc_rim_list.append(a) for b in doc_rim for a in b]

        out_df.loc[i] = [in_df.loc[i]['coh_name'], in_df.loc[i]['doc_name'], 0, 0] + coh_core_list + coh_rim_list + \
                        doc_core_list + doc_rim_list + [1 if in_df.loc[i]['binders'] else 0]
        identities_df.loc[i] = [in_df.loc[i]['coh_name'], in_df.loc[i]['doc_name'], 0, 0] + \
                               [coh_identities[v] for v in ordered_positions['coh'] if 'core' in v] + \
                               [coh_identities[v] for v in ordered_positions['coh'] if 'rim' in v] + \
                               [doc_identities[v] for v in ordered_positions['doc'] if 'core' in v] + \
                               [doc_identities[v] for v in ordered_positions['doc'] if 'rim' in v]
    return out_df, identities_df


def validate_aligned_non_aligned_interface_positions(interface_positions: dict, aligned_seqs: dict, seqs: dict) -> True:
    for name, typos in interface_positions.items():
        for typ, pos in typos.items():
            aa_non_aln_res = seqs[name][pos]
            aa_aln_pos = aligned_seqs[name].non_aligned_position_at_aligned(pos)
            aa_aln_res = aligned_seqs[name].get_aligned_positions([aa_aln_pos])[0]
            assert aa_non_aln_res == aa_aln_res, 'could not validate aligned to non aligned interface residues'


def highest_seq_similarity(crys_seqs: list, query: AASeq) -> (AASeq, float):
    """
    :param crys_seqs: list of AASeq instances of crystalised seqs
    :param query: a query AASeq
    :return: the most sequence-similar sequence
    """
    best_seq, best_iden = AASeq(), 0.0
    for seq in crys_seqs:
        iden_ = query.aligned_identity(seq)
        if iden_ > best_iden:
            best_iden = iden_
            best_seq = seq
    return best_seq, best_iden


def seq_to_binary(seq: AASeq, coh_doc: str) -> (list, list):
    """
    :param seq: query seq as list of identities in the 1ohz alignments
    :param coh_doc: either coh or doc
    :return: to list of binary, one for coh, the other for doc
    """
    if coh_doc == 'coh':
        core = [core_res_to_identity(seq[i]) for i in positions['core_coh']]
        rim = [rim_res_to_type_binary(seq[i]) for i in positions['rim_coh']]
    elif coh_doc == 'doc':
        core = [core_res_to_identity(seq[i]) for i in positions['core_doc']]
        rim = [rim_res_to_type_binary(seq[i]) for i in positions['rim_doc']]
    return core, rim


def rim_res_to_type_binary(res: str) -> list:
    """
    :param res: residues as str
    :return: binary for type
    >>> rim_res_to_type_binary('K')
    [1, 0, 0, 0]
    >>> rim_res_to_type_binary('L')
    [0, 0, 0, 1]
    """
    if res == '-':
        return [0, 0, 0, 0]
    return rim_types_to_binary[[k for k, v in type_to_res.items() if res in v][0]]


def core_res_to_identity(res: str) -> list:
    """
    :param res: a residue
    :return: [0, 1] list portraying the residue
    >>> core_res_to_identity('A')
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    >>> core_res_to_identity('Y')
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
    """
    if res == '-':
        return [0] * 20
    return [1 if aa == res else 0 for aa in aas]


def binary_to_res(binary: list) -> str:
    """
    :param binary: binary list [0, 0...1]
    :return: res as string
    >>> binary_to_res([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    'A'
    >>> binary_to_res([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    'Y'
    """
    if binary == [0] * 20:
        return '-'
    return [[k for k, v in aas_dict.items() if i+1 == v][0] for i, n in enumerate(binary) if n == 1][0]


def parse_binding_data() -> pd.DataFrame:
    """
    :return: data frame 'coh_name', 'doc_name', 'coh_seq', 'doc_seq', 'binders' for Rachel's data
    """
    from _binding_data import binding_data
    rachel_root = '/home/labs/fleishman/jonathaw/decision_tree/'
    cohs = read_multi_fastas(rachel_root+'cohesins_from_rachel_and_vered.fasta_aln', suffix_to_remove='/', lower=True)
    docs = read_multi_fastas(rachel_root+'dockerins_from_rachel_and_vered.fasta_aln', suffix_to_remove='/', lower=True)
    rachel_bind = binding_data()
    vered_bind = parse_vered_binding()
    result = pd.DataFrame(columns=['coh_name', 'doc_name', 'coh_seq', 'doc_seq', 'binders'])
    i = 1
    for coh, docs_dict in rachel_bind.items():
        for doc, res in docs_dict.items():
            result.loc[i] = [coh, doc, cohs[coh], docs[doc], rachel_bind[coh][doc]]
            i += 1
    for coh, docs_dict in vered_bind.items():
        for doc, res in docs_dict.items():
            result.loc[i] = [coh, doc, cohs[coh], docs[doc], vered_bind[coh][doc] == 1]
            i += 1
    print('there are %i rows in the data' % (i-1))
    for name in ['1ohz', '2b59', '2ozn', '2vn5', '2y3n', '3ul4', '4fl4', '4fl5', '4dh2', '4uyp', '5new']:
        result.loc[i] = [name, name, cohs[name], docs[name], True]
        i += 1

    return result


def parse_vered_binding() -> dict:
    """
    :return: {coh: {doc: 0/1}} dict of binding results from Vered
    """
    cohs = ['ScaA1', 'ScaB2', 'ScaB4', 'ScaB6', 'ScaB9', 'ScaC', 'ScaE', 'ScaF', 'ScaG', 'ScaH', 'ScaI', 'ScaJ1',
            'ScaJ2', 'ScaO']
    docs = ['1132', '1222', '3070', '3925', '4079', '4293', '3113', '4069', '341', '614', '794', '3116', '3129', '3115',
            '3114', '1965', '1541']
    result = {coh.lower(): {doc.lower(): 0 for doc in docs} for coh in cohs}
    for l in open('/home/labs/fleishman/jonathaw/decision_tree/experimental_results.txt', 'r'):
        s = l.split()
        for d in s[1:]:
            result[s[0].lower()][d] = 1
    return result


def parse_interface_positions() -> dict:
    """
    :return: parses the information in the interface position analysis I made that tells which position in every crystal
    dimer corresponds to which coh/doc core/rim position. {'coh/doc': {name: {coh/doc_core/rim_#: i}}
    """
    result, fields = {'coh': {}, 'doc': {}}, {}
    coh_doc = 'coh'
    with open(decision_tree_root+'interface_positions.txt', 'r') as fin:
        for l in fin:
            s = l.split()
            if len(s) > 1:
                if s[0] in ['coh', 'doc']:
                    fields = {a: i for i, a in enumerate(s)}
                    if s[0] == 'doc':
                        coh_doc = 'doc'
                else:
                    if coh_doc != 'doc' or s[0] not in doc_len_reduction.keys():
                        result[coh_doc][s[0].lower()] = {k: int(s[v]) for k, v in fields.items() if v != 0}
                    else:
                        result[coh_doc][s[0].lower()] = {k: int(s[v])-doc_len_reduction[s[0]] for k, v in fields.items()
                                                         if v != 0}
                    # print('parser', coh_doc, s[0], result[coh_doc][s[0].lower()])
    return result


def analyse_identity_df(idf: pd.DataFrame):
    """
    :param idf: identitites dataframe. every cannonical positions, every residue pair, and the identity
    :return:
    """
    crys_names = ['1ohz', '2b59', '2ozn', '2vn5', '2y3n', '3ul4', '4fl4', '4fl5', '4dh2', '4uyp', '5new']
    coh_columns = [a for a in idf.columns if 'coh_core' in a] + [a for a in idf.columns if 'coh_rim' in a]
    doc_columns = [a for a in idf.columns if 'doc_core' in a] + [a for a in idf.columns if 'doc_rim' in a]
    known_cohs, known_docs = {}, {}
    for i in range(1, len(idf.index)+1):
        if idf.loc[i]['coh_name'] not in known_cohs.keys() and idf.loc[i]['coh_name'].lower() in crys_names:
            known_cohs[idf.loc[i]['coh_name']] = ''.join(idf.ix[i, coh_columns].values.tolist())
        if idf.loc[i]['doc_name'] not in known_docs.keys() and idf.loc[i]['doc_name'].lower() in crys_names:
            known_docs[idf.loc[i]['doc_name']] = ''.join(idf.ix[i, doc_columns].values.tolist())

    print('cohesins:')
    for k, v in known_cohs.items():
        print('>%s\n%s' % (k, v))
    print('dockerins:')
    for k, v in known_docs.items():
        print('>%s\n%s' % (k, v))


if __name__ == '__main__':
    main()