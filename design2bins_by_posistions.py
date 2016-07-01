#!/usr/bin/env python3.5
import sys
from RosettaFilter import RunFilters, Filter, score2dict
from DoCohResultProcessor import generate_run_filters, all_who_pass_run_filters
from seq_funcs import read_multi_fastas
from AASeq import AASeq
# from DoCohResultProcessor import best_n_structures
# from rosetta_score_files import score2dict
import operator
import pickle
import os
from collections import Counter, OrderedDict
import sys
import re
import networkx as nx

type_dict = {'D': 'n', 'E': 'n', 'K': 'p', 'R': 'p'}
doc_symm_poses = OrderedDict({11: 45, 12: 46, 15: 49, 18: 52, 19: 53, 22: 56, 45: 11, 46: 12, 49: 15, 52: 18, 53: 19,
                              56: 22})
positions_dict = {'1anu': [32, 36, 62, 65, 69, 82, 115, 126],
                  '1aoh': [33, 37, 63, 66, 70, 83, 119, 130],
                  '1ohz': [33, 35, 37, 39, 63, 66, 68, 70, 73, 75, 77, 79, 81, 83, 121, 125, 127],
                  '2ccl': [33, 37, 63, 66, 70, 83, 116, 127]}
doc_symm_dict = {'1ohz': OrderedDict({11: 45, 12: 46, 15: 49, 18: 52, 19: 53, 22: 56, 45: 11, 46: 12, 49: 15, 52: 18,
                                      53: 19, 56: 22}),
                 '2vn5': OrderedDict({14: 42, 15: 43, 17: 45, 18: 46, 21: 49, 22: 50, 24: 52, 25: 53, 42: 14, 43: 15,
                                      45: 17, 46: 18, 49: 21, 50: 22, 53: 25}),
                 # '2y3n': {OrderedDict({})},
                 '3ul4': {},
                 '4dh2': {},
                 '4fl4': {},
                 '4fl5': {},
                 '4uyp': {},
                 '5new': {}}
# this dict is for normalizing between the docs for comparisons. {doc_pos: 1ohz_pos}.
doc21ohz = {'1ohz': OrderedDict({11: 11, 12: 12, 15: 15, 18: 18, 19: 19, 22: 22, 45: 45, 46: 46, 49: 49, 52: 52, 53: 53,
                                 56: 56}),
            '2vn5': OrderedDict({14: 11, 15: 12, 18: 15, 21: 18, 22: 19, 25: 22, 42: 45, 43: 46, 46: 49, 49: 52, 50: 53,
                                 53: 56}),
            '2y3n': OrderedDict({17: 11, 18: 12, 21: 15, 24: 18, 25: 19, 27: 22, 44: 45, 45: 46, 48: 49, 51: 52, 52: 53,
                                 54: 56}),
            '3ul4': OrderedDict({14: 11, 15: 12, 18: 15, 21: 18, 22: 19, 25: 22, 47: 45, 48: 46, 51: 49, 54: 52, 55: 53,
                                 58: 56}),
            '4dh2': OrderedDict({14: 11, 15: 12, 18: 15, 21: 18, 22: 19, 25: 22, 50: 45, 51: 46, 54: 49, 57: 52, 58: 53,
                                 61: 56}),
            '4fl4': OrderedDict({13: 11, 14: 12, 17: 15, 20: 18, 21: 19, 24: 22, 49: 45, 50: 46, 53: 49, 56: 52, 57: 53,
                                 60: 56}),
            '4fl5': OrderedDict({5 : 11, 6 : 12, 9 : 15, 12: 18, 13: 19, 15: 22, 38: 45, 39: 46, 42: 49, 45: 52, 46: 53,
                                 48: 56}),
            '4uyp': OrderedDict({15: 11, 16: 12, 19: 15, 22: 18, 23: 19, 26: 22, 51: 45, 52: 46, 55: 49, 58: 52, 59: 53,
                                 62: 56}),
            '5new': OrderedDict({13: 11, 14: 12, 17: 15, 20: 18, 21: 19, 24: 22, 44: 45, 45: 46, 48: 49, 51: 52, 52: 53,
                                 55: 56})
            }
ohz2doc = {k1: OrderedDict({v2: k2 for k2, v2 in v1.items()}) for k1, v1 in doc21ohz.items()}
switch_symm = OrderedDict({0: 6, 1: 7, 2: 8, 3: 9, 4: 10, 5: 11, 6: 0, 7: 1, 8: 2, 9: 3, 10: 4, 11: 5})

# clusters made by examining the resemblance of the dockerins around the hot-spot binding area
doc_bb_clusters = {1: ['1ohz', '4fl4', '4uyp', '2vn5', '5new'], 2: ['2y3n'], 3: ['4dh2', '3ul4'], 4: ['4fl5']}


class Result():
    def __init__(self, name: str, coh_AASeq: AASeq, doc_AASeq: AASeq, purples: int, j=False, originals=None, doc_wt=False):
        self.name = name
        self.coh_AASeq = coh_AASeq
        self.doc_AASeq = doc_AASeq
        self.purples = purples

        if not j:
            name_split = name.split('_')
            a_ind = name_split.index('A')
            self.coh_wt = name_split[a_ind-1]
            self.doc_wt = name_split[a_ind+1]
        elif doc_wt:
            self.coh_wt = '1ohz'
            self.doc_wt = doc_wt
        else:
            self.coh_wt = '1ohz'
            self.doc_wt = originals[name[:-3]+'.pdb.gz'].split('_A_')[1].split('_')[0]

        pos_str = coh_AASeq.get_positions(positions_dict['1ohz'])
        self.coh_switch = ''.join([type_dict[a] if a in type_dict.keys() else 'c' for a in pos_str])
        doc_str = doc_AASeq.get_positions(list(doc21ohz[self.doc_wt].keys()))
        self.doc_switch = ''.join([type_dict[a] if a in type_dict.keys() else 'c' for a in doc_str])

    def __repr__(self):
        return '%s #purples: %i coh_switch: %s, doc_switch: %s doc_wt: %s' % \
               (self.name, self.purples, self.coh_switch, self.doc_switch, self.doc_wt)

    def get_doc_rel_pos(self, pos):
        """
        :param pos: a position in the doc, from doc_symm_poses
        :return: the identity of my doc's relative position
        """
        assert pos in doc_symm_poses.keys(), 'pos %i not in keys' % pos
        return self.doc_AASeq.get_positions([ohz2doc[self.doc_wt][pos]])[0]


def make_switches(args, scores, run_filters, seq_dict):
    """
    A script that takes a fasta and score files and assembles bins. a bin is a stack of
    sequences that have similar AAs at the same positions (set by type_dict and
    positions_dict respectively. each sequence is read, and it's score is examined whether
    it passes my thresholds (purple). if it is it is assigned a bin, where only Negative
    and Positive make a difference. the set of all generated bins (basically a list of
    strings of n/p/c) is subsetted to get the longest subsets of bins that differe from
    all other bins in their subset in at least 1 position N <> P.
    INPUT: 1st cmd argument fasta file
           2nd cmd argument score file
    :return:
    """
    ### this positions dict is for the 1st, 8 parts switches design
    # positions_dict = {'1anu': [36, 38, 114, 115, 117, 120, 124, 126], '1ohz': [37, 39, 115, 116, 118, 121, 125, 127],
    #                   '2ccl': [37, 39, 115, 116, 118, 121, 125, 127]}
    ### this positions dict if for the second, 10 parts switches design from 1-2.3.2015
    positions_dict = {'1anu': [32, 36, 62, 65, 69, 82, 115, 126],
                      '1aoh': [33, 37, 63, 66, 70, 83, 119, 130],
                      '1ohz': [33, 35, 37, 39, 63, 66, 68, 70, 73, 75, 77, 79, 81, 83, 121, 125, 127],
                      '2ccl': [33, 37, 63, 66, 70, 83, 116, 127]}
    type_dict = {'D': 'n', 'E': 'n', 'K': 'p', 'R': 'p'}

    bins = {}
    num_structs_bin = {}
    name_list = list(scores.keys())
    for name, seq in seq_dict.items():
        if name[:-6] not in name_list:
            continue
        if run_filters.test_all(scores[name[:-6]]):
            coh_name = seq_dict[name].get_positions(positions_dict['1ohz'])

            switches = ''.join(type_dict[a] if a in type_dict.keys() else 'c' for a in coh_name)
            ### adding a condition where total #charges is <= 7, and distributes 5/2 or 4/3:
            # counter = Counter(switches)
            # print('chounter', counter)
            # if counter['n'] < 2 or counter['p'] < 2:
            #     continue
            ###
            # if switches.count('n') < 4 or switches.count('p') < 4:
            #     continue
            if switches not in bins.keys():
                bins[switches] = []
                num_structs_bin.update({switches: 0})
            bins[switches].append(name)
            num_structs_bin[switches] += 1
    return bins, num_structs_bin


def swithces_from_diagonal(args, run_filters, coh_seq_dict, doc_seq_dict):
    """
    :param args: run arguments
    :param run_filters: run filters
    :param coh_seq_dict: {name: AASeq()} of cohesins
    :param doc_seq_dict: {name: AASeq()} of dockerins
    :return: {switch_name: {design_name: #purples}
    """
    # positions_dict = {'1anu': [32, 36, 62, 65, 69, 82, 115, 126],
    #                   '1aoh': [33, 37, 63, 66, 70, 83, 119, 130],
    #                   '1ohz': [33, 35, 37, 39, 63, 66, 68, 70, 73, 75, 77, 79, 81, 83, 121, 125, 127],
    #                   '2ccl': [33, 37, 63, 66, 70, 83, 116, 127]}
    # type_dict = {'D': 'n', 'E': 'n', 'K': 'p', 'R': 'p'}
    results = {}
    sc_files = [a for a in os.listdir(args['score_dir']) if a[-6:] == '.score']
    bins = {}
    for sc_file in sc_files:
        score_dict = score2dict(sc_file)
        passed, failed = all_who_pass_run_filters(args, score_dict, run_filters)
        # results[sc_file] = len(list(passed.keys()))
        if len(list(passed.keys())) <= args['purples_threshold']:
            continue
        ### these kick out the date from the score names, and the makes it into the proper names:
        if '_11.10' in sc_file:
            coh_seq = coh_seq_dict[''.join(str(sc_file[4:-6]+'_0001.pdb.A').split('_11.10'))]
            doc_seq = doc_seq_dict[''.join(str(sc_file[4:-6]+'_0001.pdb.B').split('_11.10'))]
        elif '_12.10' in sc_file:
            coh_seq = coh_seq_dict[''.join(str(sc_file[4:-6]+'_0001.pdb.A').split('_12.10'))]
            doc_seq = doc_seq_dict[''.join(str(sc_file[4:-6]+'_0001.pdb.B').split('_12.10'))]

        result = Result(sc_file[4:-6], coh_seq, doc_seq, len(list(passed.keys())))

        # pos_str = coh_seq.get_positions(positions_dict['1ohz'])
        # switch_str = ''.join([type_dict[a] if a in type_dict.keys() else 'c' for a in pos_str])
        #
        # doc_str = doc_seq.get_positions(list(doc_symm_poses.keys()))
        # doc_switch = ''.join([type_dict[a] if a in type_dict.keys() else 'c' for a in doc_str])

        # if switch_str not in bins.keys():
        #     bins[switch_str] = {}
        # if doc_switch not in bins[switch_str].keys():
        #     bins[switch_str][doc_switch] = {}

        # bins[switch_str][doc_switch][''.join(str(sc_file[4:-6]+'_0001.pdb.A').split('_11.10'))] \
        #     = {'doc_seq': doc_seq, 'purples': len(list(passed.keys()))}
        results[sc_file[4:-6]] = result
    # return bins
    return results


def diagonal_bins_to_cliques(args, run_filters, bins):
    G = nx.Graph()
    [G.add_node(a) for a in bins.values()]

    for n1 in G.nodes_iter():
        for n2 in G.nodes_iter():
            if n1.name == n2.name:
                continue

            docs_diff = docs_differ_symmetry(n1.doc_AASeq, n2.doc_AASeq, '1ohz')
            if switches_differ(args, n1.coh_switch, n2.coh_switch) and docs_diff >= args['doc_diff_by']:
                G.add_edge(n1, n2)

            else:
                print('\n')
                print(n1.coh_switch)
                print(n2.coh_switch)
                print(n1.doc_switch)
                print(n2.doc_switch, docs_diff)


    cliques = [a for a in nx.find_cliques(G)]
    max_len = max([len(a) for a in cliques])
    max_cliques = [a for a in cliques if len(a) == max_len]
    print('there are %i cliques with %i structures in each for diff_by=%i doc_diff_by=%i' % (len(max_cliques), max_len, args['diff_by'], args['doc_diff_by']))
    return max_cliques


def minidiagonal_cliques(args):
    original_names = parse_name_translation('/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/recliques_4Nov/clique_6_pdbs/mini_diagonal_11Nov/minidiagonal_pdbs/translate_names.txt')
    coh_seqs = read_multi_fastas(args['coh_seqs_file'], suffix_to_remove='.A')
    doc_seqs = read_multi_fastas(args['doc_seqs_file'], suffix_to_remove='.B')
    if not os.path.isfile('all_results.obj') or not os.path.isfile('all_bins.obj'):
        print('creating bins and results')
        all_results, bins = {}, {}
        for design in coh_seqs.keys():
            r = Result(design, coh_seqs[design], doc_seqs[design], 0, j=True, originals=original_names)
            if not 2 < r.coh_switch.count('n') + r.coh_switch.count('p') < 8 and not 2 <= r.coh_switch.count('p') <= 3:
                continue
            all_results[design] = r
            d_sw = r.coh_switch+'-'+r.doc_switch+'-'+r.doc_wt
            if d_sw not in bins.keys():
                bins[d_sw] = []
            bins[d_sw].append(r)
        with open('all_results.obj', 'wb') as w_obj:
            pickle.dump(all_results, w_obj)
        with open('all_bins.obj', 'wb') as w_obj:
            pickle.dump(bins, w_obj)
    else:
        print('reading results')
        with open('all_results.obj', 'rb') as r_obj:
            all_results = pickle.load(r_obj)
        with open('all_bins.obj', 'rb') as r_obj:
            bins = pickle.load(r_obj)

    print('found %i bins' % len(bins))

    if not os.path.isfile('graph_%i_%i.obj' % (args['diff_by'], args['doc_diff_by'])):
        print('creating graph')
        G = nx.Graph()
        [G.add_node(a) for a in bins.keys()]
        print('found %i nodes' % G.number_of_nodes())
        for n1 in G.nodes_iter():
            for n2 in G.nodes_iter():
                if n1 != n2:
                    coh_sw_1, coh_sw_2 = n1.split('-')[0], n2.split('-')[0]
                    doc_sw_1, doc_sw_2 = n1.split('-')[1], n2.split('-')[1]
                    doc_wt_1, doc_wt_2 = n1.split('-')[2], n2.split('-')[2]
                    doc_diff = 1 if are_docs_from_diff_clusters(doc_wt_1, doc_wt_2) else 0
                    symm_switch = switch_symm_changer(doc_sw_2)
                    if switches_differ({'diff_by': args['diff_by']}, coh_sw_1, coh_sw_2) >= args['diff_by'] and \
                            switches_differ({'diff_by': args['doc_diff_by']}, doc_sw_1, doc_sw_2) + doc_diff >= args['doc_diff_by'] and \
                            switches_differ({'diff_by': args['doc_diff_by']}, doc_sw_1, symm_switch) + doc_diff >= args['doc_diff_by']:
                        G.add_edge(n1, n2)
        with open('graph_%i_%i.obj' % (args['diff_by'], args['doc_diff_by']), 'wb') as w_obj:
            pickle.dump(G, w_obj)
    else:
        print('reading graph')
        with open('graph_%i_%i.obj' % (args['diff_by'], args['doc_diff_by']), 'rb') as r_obj:
            G = pickle.load(r_obj)

    if not os.path.isfile('max_cliques_%i_%i.obj' % (args['diff_by'], args['doc_diff_by'])):
        cliques = [a for a in nx.find_cliques(G)]
        max_len = max([len(a) for a in cliques])
        max_cliques = [a for a in cliques if len(a) == max_len]
        print('there are %i cliques with %i structures in each for diff_by=%i doc_diff_by=%i' %
              (len(max_cliques), max_len, args['diff_by'], args['doc_diff_by']))
        with open('max_cliques_%i_%i.obj' % (args['diff_by'], args['doc_diff_by']), 'wb') as w_obj:
            pickle.dump(max_cliques, w_obj)
    else:
        print('reading cliques')
        with open('max_cliques_%i_%i.obj' % (args['diff_by'], args['doc_diff_by']), 'rb') as r_obj:
            max_cliques = pickle.load(r_obj)

    occurences = {a.name: 0 for clq in max_cliques for k in clq for a in bins[k]}
    for clq in max_cliques:
        print('in clq', clq)
        for k in clq:
            # print(bins[k][0].name)
            print('\n'.join(set([a.name for a in bins[k]])))
            for a in bins[k]:
                occurences[a.name] += 1
    for w in sorted(occurences, key=occurences.get, reverse=False):
        print(w, occurences[w])


def primary_design_to_cliques(args, run_filters, coh_seqs, doc_seqs, scores):
    """
    :param args: run arguments
    :param run_filters: run filters
    :param coh_seqs: {name: AASeq}
    :param doc_seqs: {name: AASeq}
    :param scores: {name: {filter: score}}
    :return:
    """
    ### first segment: if there's no results yet, will make them. goes over all scores,
    if not os.path.isfile('all_results.obj') or not os.path.isfile('all_bins.obj'):
        print('creating results')
        all_results = {}
        bins = {}
        passed = 0
        for score_name, score in scores.items():
            # if not run_filters.test_all(score):
            #     continue
            passed += 1
            r = Result(score_name, coh_seqs[score_name], doc_seqs[score_name], 0)
            # gate decided upon by using MSA_charge_distribution.py
            if not 2 < r.coh_switch.count('n') + r.coh_switch.count('p') < 8 and not 2 <= r.coh_switch.count('p') <= 3:
                continue
            all_results[score_name] = r
            d_sw = r.coh_switch+'-'+r.doc_switch+'-'+r.doc_wt
            if d_sw not in bins.keys():
                bins[d_sw] = []
            bins[d_sw].append(r)
        with open('all_results.obj', 'wb') as w_obj:
            pickle.dump(all_results, w_obj)
        with open('all_bins.obj', 'wb') as w_obj:
            pickle.dump(bins, w_obj)
        print('found %i designs that passed the thresholds' % passed)
    else:
        print('reading results')
        with open('all_results.obj', 'rb') as r_obj:
            all_results = pickle.load(r_obj)
        with open('all_bins.obj', 'rb') as r_obj:
            bins = pickle.load(r_obj)

    if not os.path.isfile('graph_%i_%i.obj' % (args['diff_by'], args['doc_diff_by'])):
        print('creating graph')
        G = nx.Graph()
        [G.add_node(a) for a in bins.keys()]
        print('found %i nodes' % G.number_of_nodes())
        for n1 in G.nodes_iter():
            for n2 in G.nodes_iter():
                if n1 != n2:
                    coh_sw_1, coh_sw_2 = n1.split('-')[0], n2.split('-')[0]
                    doc_sw_1, doc_sw_2 = n1.split('-')[1], n2.split('-')[1]
                    doc_wt_1, doc_wt_2 = n1.split('-')[2], n2.split('-')[2]
                    doc_diff = 1 if are_docs_from_diff_clusters(doc_wt_1, doc_wt_2) else 0
                    symm_switch = switch_symm_changer(doc_sw_2)
                    if switches_differ({'diff_by': args['diff_by']}, coh_sw_1, coh_sw_2) >= args['diff_by'] and \
                            switches_differ({'diff_by': args['doc_diff_by']}, doc_sw_1, doc_sw_2) + doc_diff >= args['doc_diff_by'] and \
                            switches_differ({'diff_by': args['doc_diff_by']}, doc_sw_1, symm_switch) + doc_diff >= args['doc_diff_by']:
                        G.add_edge(n1, n2)
        with open('graph_%i_%i.obj' % (args['diff_by'], args['doc_diff_by']), 'wb') as w_obj:
            pickle.dump(G, w_obj)
    else:
        print('reading graph')
        with open('graph_%i_%i.obj' % (args['diff_by'], args['doc_diff_by']), 'rb') as r_obj:
            G = pickle.load(r_obj)

    if not os.path.isfile('max_cliques.obj'):
        cliques = [a for a in nx.find_cliques(G)]
        max_len = max([len(a) for a in cliques])
        max_cliques = [a for a in cliques if len(a) == max_len]
        print('there are %i cliques with %i structures in each for diff_by=%i doc_diff_by=%i' %
              (len(max_cliques), max_len, args['diff_by'], args['doc_diff_by']))
        with open('cliques.obj', 'wb') as w_obj:
            pickle.dump(max_cliques, w_obj)
    else:
        print('reading cliques')
        with open('cliques.obj', 'rb') as r_obj:
            max_cliques = pickle.load(r_obj)

    for clq in max_cliques:
        print('in clq', clq)
        for k in clq:
            # print(bins[k][0].name)
            print('\n'.join([a.name for a in bins[k]]))


def post_pred_cliques(args):

    run_filters = generate_run_filters(args={'ddg': 25.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6, 'buried_2': 3})

    if not os.path.isfile('./all_data.obj'):
        sc_files = [a for a in os.listdir('./') if '.score' in a]
        cohs_seqs = read_multi_fastas('/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/recliques_4Nov/all_cohs.fasta')
        docs_seqs = read_multi_fastas('/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/recliques_4Nov/all_docs.fasta')
        results = []
        for sc_file in sc_files:
            seq_name = '_'.join(sc_file.split('_')[1:8])
            coh_name = seq_name+'.pdb.gz.A'
            doc_name = seq_name+'.pdb.gz.B'
            sc_dict = score2dict(sc_file)
            ynum = re.search('y[0-9]{3}', sc_file).group(0)
            passed, failed = all_who_pass_run_filters(args, sc_dict, run_filters)
            if len(passed) >= args['purples_threshold']:
                r = Result(seq_name, cohs_seqs[coh_name], docs_seqs[doc_name], len(passed))
                results.append(r)
        with open('./all_data.obj', 'wb') as fout:
            pickle.dump(results, fout)
    else:
        with open('./all_data.obj', 'rb') as fin:
            results = pickle.load(fin)

    if not os.path.isfile('./graph.obj'):
        result_dict = {i+1: r for i, r in enumerate(results)}
        G = nx.Graph()
        [G.add_node(a) for a in result_dict.keys()]
        for n1 in G.nodes_iter():
            for n2 in G.nodes_iter():
                if n1 != n2:
                    coh_sw_1, coh_sw_2 = result_dict[n1].coh_switch, result_dict[n2].coh_switch
                    doc_sw_1, doc_sw_2 = result_dict[n1].doc_switch, result_dict[n2].doc_switch
                    doc_wt_1, doc_wt_2 = result_dict[n1].doc_wt, result_dict[n2].doc_wt
                    doc_diff = 1 if are_docs_from_diff_clusters(doc_wt_1, doc_wt_2) else 0
                    symm_switch = switch_symm_changer(doc_sw_2)
                    if switches_differ({'diff_by': args['diff_by']}, coh_sw_1, coh_sw_2) >= args['diff_by'] and \
                            switches_differ({'diff_by': args['doc_diff_by']}, doc_sw_1, doc_sw_2) + doc_diff >= args['doc_diff_by'] and \
                            switches_differ({'diff_by': args['doc_diff_by']}, doc_sw_1, symm_switch) + doc_diff >= args['doc_diff_by']:
                        G.add_edge(n1, n2)
                        print('adding edge\n', result_dict[n1], '\n', result_dict[n2])
                    else:
                        print('NOT\n', result_dict[n1], '\n', result_dict[n2])
        cliques = [a for a in nx.find_cliques(G)]
        max_len = max([len(a) for a in cliques])
        max_cliques = [a for a in cliques if len(a) == max_len]
        for clq in max_cliques:
            print(clq, '\n', '\n'.join([str(result_dict[a]) for a in clq]))


def are_docs_from_diff_clusters(doc_1, doc_2):
    """
    :param doc_1: doc 1 WT
    :param doc_2: doc 2 WT
    :return: True if they come from differnet clusters as described in doc_bb_clusters
    >>> are_docs_from_diff_clusters('1ohz', '1ohz')
    False
    >>> are_docs_from_diff_clusters('1ohz', '4fl4')
    False
    >>> are_docs_from_diff_clusters('1ohz', '4dh2')
    True
    """
    clst_1 = [k for k, v in doc_bb_clusters.items() if doc_1 in v][0]
    clst_2 = [k for k, v in doc_bb_clusters.items() if doc_2 in v][0]
    return clst_1 != clst_2


def switch_symm_changer(switch_str):
    """
    :param switch_str: a dockerin's switch string (c,p,n)
    :return: it's symmetric conterpart
    >>> switch_symm_changer('nnnnnnpppppp')
    'ppppppnnnnnn'
    >>> switch_symm_changer('nnppnpppnnpn')
    'ppnnpnnnppnp'
    """
    sw = ['-'] * 12
    for i, s in enumerate(switch_str):
        sw[switch_symm[i]] = s
    return ''.join(sw)


def different_docs_differ_symm(res1: Result, res2: Result):
    counter = 1
    # print('\n')
    # print(res1.doc_AASeq.get_seq())
    # print(res2.doc_AASeq.get_seq())
    for p1, p2 in doc_symm_poses.items():
        counter += 1 if poses_differ(res1.get_doc_rel_pos(p1), res2.get_doc_rel_pos(p1)) and \
                        poses_differ(res1.get_doc_rel_pos(p1), res2.get_doc_rel_pos(p2)) else 0
        # print(p1, p2, res1.get_doc_rel_pos(p1), res2.get_doc_rel_pos(p1), res2.get_doc_rel_pos(p2))
        # print(poses_differ(res1.get_doc_rel_pos(p1), res2.get_doc_rel_pos(p1)) and \
        #                 poses_differ(res1.get_doc_rel_pos(p1), res2.get_doc_rel_pos(p2)))
    # print('\n')
    return counter


def docs_differ_symmetry(seq1, seq2, doc_wt):
    """
    :param args: run arguments
    :param seq1: doc seq 1
    :param seq2: doc seq 2
    :return: True if they symmetrically differ by args['doc_diff_by']
    >>> a = AASeq(string='KKKDDDAAA')
    >>> b = AASeq(string='KKKDDDAAA')
    >>> docs_differ_symmetry({'doc_diff_by': 2}, a, b)
    False
    >>> b = AASeq(string='DDKDDDAAA')
    >>> docs_differ_symmetry({'doc_diff_by': 2}, a, b)
    True
    """
    counter = 0
    for p1, p2 in doc_symm_dict[doc_wt].items():
        counter += 1 if poses_differ(seq1[p1], seq2[p1]) and poses_differ(seq1[p1], seq2[p2]) else 0
        # print(p1, p2, seq1[p1], seq2[p1], seq2[p2], poses_differ(seq1[p1], seq2[p1]) and poses_differ(seq1[p1], seq2[p2]))
    return counter


def poses_differ(p1, p2):
    """
    :param p1: aa1
    :param p2: aa2
    :return: True is differ in charge
    >>> poses_differ('K', 'D')
    True
    >>> poses_differ('K', 'R')
    False
    >>> poses_differ('K', 'c')
    False
    """
    return (p1 in ['K', 'R'] and p2 not in ['K', 'R']) or (p1 in ['D', 'E'] and p2 not in ['D', 'E'])


def best_struct_in_bin(args: dict, score_dict: dict, switches: dict) -> dict:
    """
    :param args: run arguments
    :param switches: swithches dict {switch_name: [strucutres]}
    :return: {swithch: best n structs in switch bin}
    """
    switches_best = {}
    for k, v in switches.items():
        switch_scores = {name: score for name, score in score_dict.items() if name+'.pdb.A' in v}
        switches_best[k] = best_n_structures({'filter': 'ddg', 'n': args['n']}, switch_scores)
    return switches_best


def best_cliques(args: dict, switch_names: list):
    """
    :param args: run arguments
    :param switches: [switch names]
    :return: the biggest clique with args['diff_by'] charge differences
    """
    import networkx as nx
    G = nx.Graph()
    [G.add_node(sw) for sw in switch_names]
    for sw1 in G.nodes_iter():
        for sw2 in G.nodes_iter():
            if sw1 != sw2:
                if switches_differ(args, sw1, sw2):
                    G.add_edge(sw1, sw2)
    cliques = [a for a in nx.find_cliques(G)]
    max_len = max([len(a) for a in cliques])
    max_cliques = [a for a in cliques if len(a) == max_len]
    print('there are %i cliques with %i structures in each for diff_by=%i' % (len(max_cliques), max_len, args['diff_by']))
    return max_cliques


def switches_differ(args: dict, sw1: str, sw2: str):
    """
    :param args: run arguments
    :param sw1: switch name 1
    :param sw2: swithc name 2
    :return: True iff switches differ by at least args['diff_by']
    >>> switches_differ({'diff_by': 5}, 'nppnppn', 'pppppnp')
    False
    >>> switches_differ({'diff_by': 5}, 'ppppnnnn', 'nnnnpppp')
    True
    """
    diff = 0
    for i in range(len(sw1)):
        # if i == 3:
        #     continue
        if (sw1[i] == 'p' and sw2[i] == 'n') or (sw1[i] == 'n' and sw2[i] == 'p'):
            diff += 1
    return diff


def choose_clique_by_seq_identity(args, switches_dict, cliques, coh_seqs, doc_seqs):
    """
    :param args: run arguments
    :param cliques: [[name1, name2,], [name3, name3]] list of lists of names of different cliques
    :param coh_seqs: {name: AASeq()}
    :param doc_seqs: {name: AASeq()}
    :return: best 5 cliques, chosen by seq identity difference
    """
    from numpy import mean
    identity_clq = []
    for i, clq in enumerate(cliques):
        coh_identities, doc_identities = [], []
        ### because there are multiple structures in each bin, I examine only the first one for determining seq identity
        clq_names = [switches_dict[name][0] for name in clq]
        for n1 in clq_names:
            for n2 in clq_names:
                if n1 != clq_names:
                    coh_identities.append(coh_seqs[n1].non_aligned_identity(coh_seqs[n2]))

                    dn1, dn2 = n1[:-1], n2[:-1]
                    dn1 += 'B'
                    dn2 += 'B'
                    doc_identities.append(doc_seqs[dn1].non_aligned_identity(doc_seqs[dn2]))
        identity_clq.append({'clique': clq, 'coh_iden_avg': mean(coh_identities), 'doc_iden_avg': mean(doc_identities)})
    return identity_clq
    #
    # identity_clq.sort(key=lambda k: (k['doc_iden_seq'], k['coh_iden_seq']))
    # for i in identity_clq[:5]:
    #     print(i)


def choose_orthogonal_cliques(args, best_cliques):
    """
    :param args: run arguments
    :param best_cliques: the best cliques be average doc seq identity
    :return: cliques
    """


def binnem(args, bins):
    bin_subsets = bin_subsetter(list(bins.keys()))
    bin_subsets_struct = {i: [num_structs_bin[j] for j in bini] for i, bini in enumerate(bin_subsets)}
    least_ones = 1000000
    chosen_bsss = 'no bin subsets'
    for bsss_key, bsss_val in bin_subsets_struct.items():
        if len([i for i in bsss_val if i == 1]) < least_ones:
            chosen_bsss = bsss_key
            least_ones = len([i for i in bsss_val if i == 1])

    best_bins_structs = []
    for best_bin in bin_subsets[chosen_bsss]:
        score_list = [scores[i] for i in bins[best_bin]]
        score_list.sort(key=operator.itemgetter('ddg'))
        best_bins_structs.append({best_bin: [i['description'] for i in score_list]})
    for biner in best_bins_structs:
        print(biner.keys()[0].upper())
        for struct in biner.values()[0]:
            print(struct)
        print('\n')


def bin_subsetter(bins_names):
    from random import shuffle
    from operator import itemgetter
    from itertools import groupby
    bins_sets = []
    already_picked = []
    for i in range(100000):
        bins_rearrange = bins_names[:]
        shuffle(bins_rearrange)
        iter_set = [bins_rearrange.pop()]
        # while iter_set[0] in already_picked:
        #     bins_rearrange = bins_names[:]
        #     shuffle(bins_names)
        #     iter_set = [bins_rearrange.pop()]
        for bin_i in bins_rearrange:
            if seq_differ_set(bin_i, iter_set):
                iter_set.append(bin_i)
        bins_sets.append(sorted(iter_set))
        already_picked.extend(iter_set)
        # if all([bn in already_picked for bn in bins_names]):
        #     print 'all were already picked', len(already_picked)
        #     break
    # sort and throw duplicate bins
    sorts = sorted(bins_sets, key=len, reverse=True)
    sorts = list(map(itemgetter(0), groupby(sorts)))
    return [i for i in sorts if len(i) >= len(sorts[0])]


def seq_differ_set(seq, seq_set):
    return all([seq_differ_seq(seq, i) for i in seq_set])


def seq_differ_seq(seq1, seq2):
    return any([seq1[i] != seq2[i] for i in range(len(seq1)) if seq1[i] != 'c' and seq2[i] != 'c'])


def differ_by_n_p_from_all(query, seq_set):
    for seq in seq_set:
        for i in range(len(query)):
            if not(query[i] == 'n' and seq[i] == 'p' or query[i] == 'p' and seq[i] == 'n'):
                return False
    return True


def differ_by_n_p_2_poses(query, seq_set):
    for seq in seq_set:
        diff = 0
        diff_pos = []
        for i in range(len(query)):
            if query[i] == 'n' and seq[i] == 'p' or query[i] == 'p' and seq[i] == 'n':
                diff_pos.append(i)
                diff += 1
        if diff >= 2:
            print('%s != %s' % (query, seq))
            return True
    return False


def parse_name_translation(file_name) -> dict:
    """
    :param file_name: full path to name translation file
    :return: {name: original_name}
    """
    result = {}
    for l in open(file_name, 'r'):
        s = l.split()
        if len(s) == 2:
            result[s[1]] = s[0]
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-score_file')
    parser.add_argument('-coh_seqs_file')
    parser.add_argument('-doc_seqs_file')
    parser.add_argument('-mode')
    parser.add_argument('-n', type=int, default=1)
    parser.add_argument('-diff_by', type=int, default=2)
    parser.add_argument('-doc_diff_by', type=int, default=1)
    parser.add_argument('-score_dir', type=str, default='./')
    parser.add_argument('-purples_threshold', type=int, default=50)
    args = vars(parser.parse_args())

    if args['mode'] != 'bins_diagonal' and args['mode'] != 'post_pred_cliques' and args['mode'] != 'minidiagonal_cliques':
        scores = score2dict(args['score_file'])
        run_filters = ''#generate_run_filters()
        coh_seq_dict = read_multi_fastas(args['coh_seqs_file'], suffix_to_remove='.pdb.gz')
        doc_seq_dict = read_multi_fastas(args['doc_seqs_file'], suffix_to_remove='.pdb.gz')

    if args['mode'] == 'switches_n_cliques':
        switches, num_bins = make_switches(args, scores, run_filters, coh_seq_dict)
        max_cliques = best_cliques(args, list(switches.keys()))

        with open('switches.obj', 'wb') as sw_file:
            pickle.dump(switches, sw_file)
        with open('max_cliques.obj', 'wb') as clq_file:
            pickle.dump(max_cliques, clq_file)

    elif args['mode'] == 'choose_by_identity':
        with open('switches.obj', 'rb') as sw_in:
            switches = pickle.load(sw_in)
        with open('max_cliques.obj', 'rb') as clq_file:
            max_cliques = pickle.load(clq_file)
        iden_cliq_list = choose_clique_by_seq_identity(args, switches, max_cliques, coh_seq_dict, doc_seq_dict)
        print('found %i in iden clique' % len(iden_cliq_list))
        with open('iden_clq.obj', 'wb') as iden_file:
            pickle.dump(iden_cliq_list, iden_file)

    elif args['mode'] == 'best_by_identity':
        with open('iden_clq.obj', 'rb') as iden_file:
            iden_cliques = pickle.load(iden_file)
        iden_cliques.sort(key=lambda k: (k['doc_iden_avg'], k['coh_iden_avg']))
        with open('switches.obj', 'rb') as sw_in:
            switches = pickle.load(sw_in)
        all_of_them = []
        for clique in iden_cliques:
            for swi_str in clique['clique']:
                for struct in switches[swi_str]:
                    all_of_them.append(struct[:-2])
        all_of_them = list(set(all_of_them))
        for i in all_of_them:
            print(i)

    elif args['mode'] == 'bins_diagonal':
        run_filters = generate_run_filters()
        coh_seq_dict = read_multi_fastas(args['coh_seqs_file'])
        doc_seq_dict = read_multi_fastas(args['doc_seqs_file'])

        if not os.path.isfile('bins.obj'):
            bins = swithces_from_diagonal(args, run_filters, coh_seq_dict, doc_seq_dict)
            with open('bins.obj', 'wb') as b_file:
                pickle.dump(bins, b_file)
        else:
            with open('bins.obj', 'rb') as b_file:
                bins = pickle.load(b_file)

        # for k, v in bins.items():
        #     print(k, v)

        if not os.path.isfile('best_clq.obj'):
            best_clique = diagonal_bins_to_cliques(args, run_filters, bins)
            with open('best_clq.obj', 'wb') as c_file:
                pickle.dump(best_clique, c_file)
        else:
            with open('best_clq.obj', 'rb') as c_file:
                best_clique = pickle.load(c_file)

    elif args['mode'] == 'primary_design_to_cliques':
        coh_seq_dict = read_multi_fastas(args['coh_seqs_file'], suffix_to_remove='.pdb.gz')
        doc_seq_dict = read_multi_fastas(args['doc_seqs_file'], suffix_to_remove='.pdb.gz')
        # run_filters = generate_run_filters(args={'ddg': 25.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6})
        run_filters = {}
        scores = {}
        primary_design_to_cliques(args, run_filters, coh_seq_dict, doc_seq_dict, scores)

    elif args['mode'] == 'test':
        with open('all_bins.obj', 'rb') as r_obj:
            bins = pickle.load(r_obj)
        for k, v in bins.items():
            print(k, len(v))

    elif args['mode'] == 'post_pred_cliques':
        post_pred_cliques(args)

    elif args['mode'] == 'minidiagonal_cliques':
        minidiagonal_cliques(args)

    else:
        print('no mode found')
    # first_clique = max_cliques[0]
    # best_switches_clique = {name: switches[name] for name in first_clique}
    # for k, v in best_switches_clique.items():
    #     print(v[0][:-2])
    # switches_best = best_struct_in_bin(args, scores, switches)