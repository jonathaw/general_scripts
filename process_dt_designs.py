#!/usr/bin/env python3.5
"""
a script to work with decision tree for clique design
"""
import argparse
import os
import subprocess
import pandas as pd
import pickle
import networkx as nx


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode')
    args = vars(parser.parse_args())

    if args['mode'] == 'compare_minidiagonal_decision_tree':
        """
        goes over the decision tree and minimal diagonal predicitions and makes a file with the intersection (predicted
        to bind by both)
        """
        mini_diagonal = set(parse_mini_diagonal_results())
        dt_result = set(parse_decision_tree_results())
        intersection = mini_diagonal.intersection(dt_result)
        with open('/home/labs/fleishman/jonathaw/decision_tree/design_data/mini_diagonal_with_decision_tree_prediciton.txt', 'w+') as fout:
            fout.write('\n'.join(intersection))

    elif args['mode'] == 'predict_all_vs_all':
        """
        make all Vs. all prediciton for designs that pass decision tree and minidiagonal
        """
        job_root = '/home/labs/fleishman/jonathaw/decision_tree/design_data/all_vs_all_decision_tree_6Jan/'
        passed_diagonal = open('/home/labs/fleishman/jonathaw/decision_tree/design_data/mini_diagonal_with_decision_tree_prediciton.txt',
                               'r').read().split('\n')
        with open(job_root+'command', 'w+') as cmd:
            for dzn in passed_diagonal:
                with open('%sjob_%s.j' % (job_root, dzn), 'w+') as fout:
                    fout.write('#!/bin/bash\n')
                    fout.write('. /usr/share/lsf/conf/profile.lsf\n')
                    fout.write('cd %s\n' % os.getcwd())
                    fout.write('/apps/RH6U4/python/3.5.0/bin/python3.5 /home/labs/fleishman/jonathaw/scripts/general_scripts/binding_decision_tree2.py -mode predict_by_coh -coh_name %s' % dzn)
                subprocess.call(['chmod', '+x', '%sjob_%s.j' % (job_root, dzn)])
                cmd.write(str('bsub -N -u /dev/null -G fleishman-wx-grp-lsf -q fleishman -o out.%s -e err.%s %s \n' %
                              (dzn, dzn, '%sjob_%s.j' % (job_root, dzn))))
        subprocess.call(['chmod', '+x', '%scommand' % job_root])

    elif args['mode'] == 'all_vs_all':
        job_root = '/home/labs/fleishman/jonathaw/decision_tree/design_data/all_vs_all_decision_tree_6Jan/'
        if not os.path.isfile(job_root+'all_vs_all_df.obj'):
            df = get_dt_all_vs_all_results()
            with open(job_root+'all_vs_all_df.obj', 'wb') as obj:
                pickle.dump(df, obj)
        else:
            with open(job_root+'all_vs_all_df.obj', 'rb') as obj:
                df = pickle.load(obj)

        if not os.path.isfile(job_root+'all_cliques.obj'):
            all_cliques, max_cliques = find_cliques(df)
            with open(job_root+'all_cliques.obj', 'wb') as obj:
                pickle.dump(all_cliques, obj)
            with open(job_root+'max_cliques.obj', 'wb') as obj:
                pickle.dump(max_cliques, obj)
        else:
            with open(job_root+'all_cliques.obj', 'rb') as obj:
                all_cliques = pickle.load(obj)
            with open(job_root+'max_cliques.obj', 'rb') as obj:
                max_cliques = pickle.load(obj)


def find_cliques(df_):
    G = nx.Graph()
    print(df_.index)
    # add node for every coh-doc predicted to bind
    for coh in df_.index:
        for doc in df_.index:
            if df_[doc][coh] == 1:
                G.add_node((coh, doc))

    # add edge if no cross-binding
    for n1 in G.nodes_iter():
        coh1, doc1 = n1[0], n1[1]
        for n2 in G.nodes_iter():
            coh2, doc2 = n2[0], n2[1]
            if n1 != n2:
                if df_[doc1][coh2] != 1 and df_[doc2][coh1] != 1:
                    G.add_edge(n1, n2)
    cliques = nx.find_cliques(G)
    max_clique_len = max([len(c) for c in cliques])
    max_cliques = [c for c in cliques if len(c) == max_clique_len]
    print('AAAAA', max_cliques)
    return cliques, max_cliques


def get_dt_all_vs_all_results() -> pd.DataFrame:
    """
    :return: cohs by docs data frame with 1s and 0s
    """
    all_vs_all_root = '/home/labs/fleishman/jonathaw/decision_tree/design_data/all_vs_all_decision_tree_6Jan/'
    all_vs_all_files = [a for a in os.listdir(all_vs_all_root) if '.txt' in a]
    df_ = pd.DataFrame()
    for all_vs_all_file in all_vs_all_files:
        df_ = df_.append(parse_dt_all_vs_all_result(all_vs_all_root, all_vs_all_file))
    print(df_)
    return df_


def parse_dt_all_vs_all_result(root_path: str, file_name: str) -> dict:
    """
    :param root_path: dir path
    :param file_name: file name
    :return: {coh: {doc: 1/0}}
    """
    res = {}
    for l in open(root_path+file_name, 'r'):
        s = l.split()
        if len(s) == 4:
            if '...' in s:
                continue
            res[s[2]] = int(s[3])
    return pd.DataFrame(index=[s[1]], columns=res.keys(), data=res)


def parse_mini_diagonal_results():
    with open('/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/recliques_4Nov/clique_6_pdbs/mini_diagonal_11Nov/minidiagonal.txt', 'r') as fin:
        cont = fin.read().split('\n')
    result = []
    for l in cont:
        s = l.split()
        if len(s) == 2:
            name = s[0].split('_dz')[0]
            num = int(s[1])
            if num > 5:
                result.append(name)
    return result


def parse_decision_tree_results():
    with open('/home/labs/fleishman/jonathaw/decision_tree/design_data/diagonal_prediciton.txt', 'r') as fin:
        return [a.split()[1] for a in fin.read().split('\n') if 'coh_name' not in a]


if __name__ == '__main__':
    main()
