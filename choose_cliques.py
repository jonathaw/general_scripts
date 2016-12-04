#!/usr/bin/env python3.5
"""
a script to analyse design2bins clique output, and choose cliques such that the minimal number of desings covers the
maximal number of cliques
"""
import os
import pickle
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx

from design2bins_by_posistions import Result, parse_name_translation, switches_differ, are_docs_from_diff_clusters, switch_symm_changer
from DoCohResultProcessor import generate_run_filters, all_who_pass_run_filters, what_coh, what_doc, show_prediction_heat_map
from seq_funcs import read_multi_fastas
from RosettaFilter import score2dict
__author__ = 'jonathan'


def using_graph():
    """
    reads files that start with analyse. the number analyse_clique_#.txt will be the number of structs in the clique.
    chooses designs that create as many cliques as possible
    :return:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-weight_cutoff', type=int)
    parser.add_argument('-maximal', default=False)
    args = vars(parser.parse_args())
    weight_cutoff = args['weight_cutoff']
    if not os.path.isfile('graph_%i.obj' % weight_cutoff):
        # analyse_files = [a for a in os.listdir('./') if 'analyse' in a and '.obj' not in a]
        analyse_files = ['analyse_clique_6.txt']
        analyses, num = {}, 1
        for f in analyse_files:
            parsed = parse_anlyse(f)
            for l in parsed:
                analyses[num] = l
                num += 1
        G = nx.Graph()
        [G.add_node(a) for a in analyses.keys()]

        for n1 in G.nodes_iter():
            for n2 in G.nodes_iter():
                if n1 != n2:
                    wt = len(set(analyses[n1]) & set(analyses[n2]))
                    if wt > weight_cutoff:
                        G.add_edge(n1, n2, weight=wt)
        print('finished building graph with %i nodes and %i edges' % (G.number_of_nodes(), G.number_of_edges()))
        with open('graph_%i.obj' % weight_cutoff, 'wb') as out:
            pickle.dump(G, out)
        with open('analyses_%i.obj' % weight_cutoff, 'wb') as out:
            pickle.dump(analyses, out)
    else:
        print('reading graph')
        with open('graph_%i.obj' % weight_cutoff, 'rb') as fin:
            G = pickle.load(fin)
        with open('analyses_%i.obj' % weight_cutoff, 'rb') as fin:
            analyses = pickle.load(fin)
        print('finished reading graph with %i nodes and %i edges' % (G.number_of_nodes(), G.number_of_edges()))

    if args['maximal']:
        if not os.path.isfile('clq_size_maximal_%i.obj' % weight_cutoff):
            clq_size, clq_grade = [], []
            for clq in nx.find_cliques(G):
                clq_size.append(len(clq))
                clq_grade.append(len(set([a for b in clq for a in analyses[b]])))
            with open('clq_size_maximal_%i.obj' % weight_cutoff, 'wb') as fout:
                pickle.dump(clq_size, fout)
            with open('clq_grade_maximal_%i.obj' % weight_cutoff, 'wb') as fout:
                pickle.dump(clq_grade, fout)
        else:
            with open('clq_size_maximal_%i.obj' % weight_cutoff, 'rb') as fin:
                clq_size = pickle.load(fin)
            with open('clq_grade_maximal_%i.obj' % weight_cutoff, 'rb') as fin:
                clq_grade = pickle.load(fin)

    if not os.path.isfile('clq_size_%i.obj' % weight_cutoff):
        clq_size, clq_grade = [], []
        for clq in nx.enumerate_all_cliques(G):
            clq_size.append(len(clq))
            clq_grade.append(len(set([a for b in clq for a in analyses[b]])))
        with open('clq_size_%i.obj' % weight_cutoff, 'wb') as fout:
            pickle.dump(clq_size, fout)
        with open('clq_grade_%i.obj' % weight_cutoff, 'wb') as fout:
            pickle.dump(clq_grade, fout)
    else:
        with open('clq_size_%i.obj' % weight_cutoff, 'rb') as fin:
            clq_size = pickle.load(fin)
        with open('clq_grade_%i.obj' % weight_cutoff, 'rb') as fin:
            clq_grade = pickle.load(fin)
    plt.scatter(clq_grade, clq_size)
    plt.show()


def another_option():
    analyse_files = ['analyse_clique_6.txt']
    analyses, num = {}, 1
    for f in analyse_files:
        parsed = parse_anlyse(f)
        for l in parsed:
            analyses[num] = l
            num += 1


def parse_anlyse(analyse_file: str) -> list:
    with open(analyse_file, 'r') as fin:
        cont = fin.read().split('in clq')
    cliques = []
    for l in cont:
        cliques.append(l.split('\n')[1:-1])
    return cliques


def choose_by_cliques():
    with open('/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/recliques_4Nov/cliques.obj', 'rb') as fin:
        from_pickle = pickle.load(fin)
    with open('/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/recliques_4Nov/all_bins.obj', 'rb') as r_obj:
            bins = pickle.load(r_obj)
    dirs = [a.split('_y')[0] for a in
            os.listdir('/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/recliques_4Nov/clique_6_pdbs/diagonal_for_cliques_8Nov')
            if 'dz_1ohz' in a]
    the_set = set()
    for clq in from_pickle:
        clq_names = [a.name for k in clq for a in bins[k]]
        overlap = [a for a in clq_names if a in dirs]
        if len(overlap) < 4:
            the_set = the_set ^ set(clq_names)
    print(the_set)
    print(len(the_set))
    for k in the_set:
        print(k)


def parse_cliques_lists(file_name: str, remove: str=None) -> list:
    """
    :param remove:
    :param file_name: file name for design2bins run (clq in [] names)
    :return: list of lists of cliques
    """
    with open(file_name, 'r') as fin:
        cont = fin.read()
    results = []
    for para in cont.split('in clq'):
        results.append([])
        for l in para.split('\n')[1:]:
            if l != '':
                if remove is not None:
                    results[-1].append(l.replace(remove, ''))
                else:
                    results[-1].append(l)
    return results


def all_i_know():
    results_path = '/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/reclique_18Nov/cliques_prediction/results/'
    sc_files = [a for a in os.listdir(results_path) if '.score' in a]
    result = []
    for sc in sc_files:
        coh_name = what_coh(sc, args={'naming': 'coh_on_doc'})
        doc_name = what_doc(sc, args={'naming': 'coh_on_doc'})
        result.append([coh_name, doc_name])
    return result


def creat_coh_doc_purples():
    results_path = '/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/reclique_18Nov/cliques_prediction/results/'
    sc_files = [a for a in os.listdir(results_path) if '.score' in a]
    run_filters = generate_run_filters(args={'ddg': 24.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6, 'buried_2': 3,
                                             'hbonds': -10.})

    if not os.path.isfile(results_path+'analysed.obj'):
        coh_doc_purples = {}
        for sc_file in sc_files:
            sc_dict = score2dict(results_path+sc_file)
            passed, failed = all_who_pass_run_filters({}, sc_dict, run_filters)
            coh_name = what_coh(sc_file, args={'naming': 'coh_on_doc'})
            doc_name = what_doc(sc_file, args={'naming': 'coh_on_doc'})
            if coh_name not in coh_doc_purples.keys():
                coh_doc_purples[coh_name] = {}
            coh_doc_purples[coh_name][doc_name] = len(passed)
        pickle.dump(coh_doc_purples, open(results_path+'analysed.obj', 'wb'))

    else:
        print("reading coh_doc_purples")
        coh_doc_purples = pickle.load(open(results_path+'analysed.obj', 'rb'))
    return coh_doc_purples


def predict_who():
    """
    :return: goes over all results and print for every clique in the design2bins output what is up with it.
     at the bottom also prints the ones that are finished.
    """
    coh_doc_purples = creat_coh_doc_purples()
    clique_list = parse_cliques_lists('/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/reclique_18Nov/stabilisation/all_stabilised/cliques_2_1.txt', remove='_st')
    finished_diagonal, finished_all = [], []
    for clq in clique_list:
        diagonal = []
        missing = []
        off_diagonal = []

        # test the diagonal
        for mem in clq:
            if mem in coh_doc_purples.keys():
                if mem in coh_doc_purples[mem].keys():
                    diagonal.append(coh_doc_purples[mem][mem] >= 10)
            else:
                missing.append(mem)

        # test the off-diagonal
        for mem1 in clq:
            for mem2 in clq:
                if mem1 != mem2:
                    if mem1 in coh_doc_purples.keys():
                        if mem2 in coh_doc_purples[mem1].keys():
                            off_diagonal.append(coh_doc_purples[mem1][mem2] >= 10)

                    if mem2 in coh_doc_purples.keys():
                        if mem1 in coh_doc_purples[mem2].keys():
                            off_diagonal.append(coh_doc_purples[mem2][mem1] >= 10)

        if all(diagonal) and diagonal:
            # print('clq %s missing these %s' % (', '.join([a for a in clq]), ','.join(missing)), diagonal, off_diagonal)
            # print()
            # if len(missing) == 1:
            #     print("MISSING ONE !!!!!!!!!!!!!!!!!!")
            # if len(diagonal) >= 4 and not any(off_diagonal):
            #     finished_diagonal.append(clq)
            # if len(diagonal) >= 4 and len(off_diagonal) >= len(diagonal)**2-len(diagonal) and all(off_diagonal):
            #     print('for clq %s found this in diagonal %s and off-diagonal %s' % (', '.join([a for a in clq]), diagonal, off_diagonal))
            #     finished_all.append(clq)
            if not any(off_diagonal) and 0 < len(off_diagonal) < 12:
                print('all off diagonal False', clq, len(off_diagonal), off_diagonal)
    print('found these finished diagonal')
    print('\n'.join([str(a) for a in finished_diagonal]))
    print('found these to be finished diagonal and all')
    print('\n'.join([str(a) for a in finished_all]))


def find_predicted_cliques():
    """
    :return:[[member1, member2....], ...] all cliques that are completely predicted and are orthogonal
    """
    coh_doc_purples = creat_coh_doc_purples()
    G = nx.Graph()
    all_cohs = list(coh_doc_purples.keys())
    all_docs = list(set([doc for coh in all_cohs for doc in coh_doc_purples[coh].keys()]))
    # [G.add_node((coh, doc)) for coh in all_cohs for doc in all_docs if coh_doc_purples[coh][doc] >= 10]
    for coh in all_cohs:
        for doc in all_docs:
            if doc in coh_doc_purples[coh].keys():
                if coh_doc_purples[coh][doc] >= 10:
                    G.add_node((coh, doc))

    for coh1, doc1 in G.nodes_iter():
        for coh2, doc2 in G.nodes_iter():
            if (coh1, doc1) != (coh2, doc2):
                if doc1 in coh_doc_purples[coh2].keys() and doc2 in coh_doc_purples[coh1].keys():
                    if coh_doc_purples[coh1][doc2] < 10 and coh_doc_purples[coh2][doc1] < 10:
                        G.add_edge((coh1, doc1), (coh2, doc2))

    cliques = list(nx.find_cliques(G))
    print('found the following cliques:')
    for clq in cliques:
        print(clq, len(clq))

    print('the grapg had %i nodes, and %i edges' % (G.number_of_nodes(), G.number_of_edges()))
    return cliques


def analyse_cliques(cliques):
    """
    :param cliques: list of cliques
    :return: prints an anlysis
    """
    coh_seqs = read_multi_fastas('/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/reclique_18Nov/stabilisation/all_stabilised/all_j_st_cohs.fasta', '_st.A')
    doc_seqs = read_multi_fastas('/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/reclique_18Nov/stabilisation/all_stabilised/all_j_st_docs.fasta', '_st.B')
    cliques_by_charges = parse_cliques_lists('/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/reclique_18Nov/stabilisation/all_stabilised/cliques_2_1.txt', remove='_st')
    original_names = parse_name_translation('/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/recliques_4Nov/clique_6_pdbs/mini_diagonal_11Nov/minidiagonal_pdbs/translate_names.txt')
    clqs_by_len = {k: [] for k in range(1, 8)}
    for clq in cliques:
        clqs_by_len[len(clq)].append(clq)
    designs_in_all_clqs = []
    for length in range(10, 6, -1):
        if length not in clqs_by_len.keys():
            continue
        for clq in clqs_by_len[length]:
            coh_diffs, doc_diffs, doc_diffs_symm, doc_bb_diffs = [], [], [], []
            print('\n\n\nfor clq', clq)
            for mem1 in clq:
                designs_in_all_clqs.append(mem1)
                wt_doc = original_names[mem1[1]+'.pdb.gz'][10:14]
                res1 = Result(mem1, coh_seqs[mem1[0]], doc_seqs[mem1[1]], 1, j=True, doc_wt=wt_doc)
                res1_doc_symm = switch_symm_changer(res1.doc_switch)
                for mem2 in clq:
                    if mem1 != mem2:
                        wt_doc = original_names[mem2[1]+'.pdb.gz'][10:14]
                        res2 = Result(mem2, coh_seqs[mem2[0]], doc_seqs[mem2[1]], 1, j=True, doc_wt=wt_doc)
                        doc_bb_diffs.append(are_docs_from_diff_clusters(res1.doc_wt, res2.doc_wt))
                        coh_diffs.append(switches_differ({}, res1.coh_switch, res2.coh_switch))
                        doc_diffs.append(switches_differ({}, res1.doc_switch, res2.doc_switch))
                        doc_diffs_symm.append(switches_differ({}, res1_doc_symm, res2.doc_switch))
                        print('results', res1)
                        print('results', res2)
                        print('docs diff', switches_differ({}, res1.doc_switch, res2.doc_switch))
                        print('doc symm diff', switches_differ({}, res1_doc_symm, res2.doc_switch))
                        print('doc BB dif', are_docs_from_diff_clusters(res1.doc_wt, res2.doc_wt))
                        print('coh diff', switches_differ({}, res1.coh_switch, res2.coh_switch))
                        N
            print('for clq %r found the following results:' % clq)
            print('doc_bb_diffs', doc_bb_diffs)
            print('doc_diffs', doc_diffs)
            print('doc_diffs_symm', doc_diffs_symm)
            print('coh_diffs', coh_diffs)
            print('total', sum([1 for a in doc_bb_diffs if a] + doc_diffs + doc_diffs_symm + coh_diffs))
    all_cohs = list(set([a[0] for a in designs_in_all_clqs]))
    all_docs = list(set([a[1] for a in designs_in_all_clqs]))
    print('these are all the cohs i need: %s, total %i' % (', '.join(all_cohs), len(all_cohs)))
    print('these are all the docs i need: %s, total %i' % (', '.join(all_docs), len(all_docs)))
    print('LONGEST CLIQUES FOUND ARE %i' % max([len(clq) for clq in cliques]))
    coh_doc_purples = creat_coh_doc_purples()
    for clq in clqs_by_len[max(list(clqs_by_len.keys()))]:
        print('clq', clq)
        cohs = [a[0] for a in clq]
        docs = [a[1] for a in clq]
        df = pd.DataFrame(index=docs, columns=cohs)
        for coh in cohs:
            for doc in docs:
                df[coh][doc] = coh_doc_purples[coh][doc]
        show_prediction_heat_map(df)

def draw_heatmap_for_chosen_designs(cliques):
    chosen_desins = ['j5517', 'j4518', 'j3622', 'j3983', 'j1647', 'j4286', 'j4653', 'j5711', 'j5093', 'j3626', 'j4398', 'j829', 'j5106', 'j1526']
    coh_doc_purples = creat_coh_doc_purples()
    df = pd.DataFrame(index=chosen_desins, columns=chosen_desins)
    for coh in chosen_desins:
        for doc in chosen_desins:
            df[coh][doc] = coh_doc_purples[coh][doc]
    print(df)
    show_prediction_heat_map(df)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', default='choose_by_cliques')
    args = vars(parser.parse_args())

    if args['mode'] == 'choose_by_cliques':
        #obsulete
        print('obsulete')
        choose_by_cliques()

    elif args['mode'] == 'predict_who':
        predict_who()

    elif args['mode'] == 'find_predicted_cliques':
        find_predicted_cliques()

    elif args['mode'] == 'analyse_cliques':
        cliques = find_predicted_cliques()
        analyse_cliques(cliques)

    elif args['mode'] == 'draw_chosen':
        cliques = find_predicted_cliques()
        draw_heatmap_for_chosen_designs(cliques)

    else:
        print('no mode found')
