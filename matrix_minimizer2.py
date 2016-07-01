def main():
    from pandas import DataFrame, Series
    from rosetta_score_files import how_many_purples_in_file
    import os
    import re
    # from matplotlib import pyplot as plt
    import networkx as nx
    score_file_list = [x for x in os.listdir('./')
                       if re.match('.*\.score', x)]
    coh_name_list = sorted(list(set([a.split('_')[1] for a in score_file_list])))
    doc_name_list = sorted(list(set([a.split('_')[3] for a in score_file_list])))
    df = DataFrame({coh_name: Series([-1 * len(doc_name_list)], index=doc_name_list) for coh_name in coh_name_list})
    df_true_score = DataFrame({coh_name: Series([-1 * len(doc_name_list)], index=doc_name_list) for coh_name in coh_name_list})
    for score_file in score_file_list:
        coh_name = score_file.split('_')[1]
        doc_name = score_file.split('_')[3]
        purple_num = int(how_many_purples_in_file('./'+score_file))
        df[coh_name][doc_name] = 1 if purple_num >= 10 else 0
        df_true_score[coh_name][doc_name] = purple_num

    G = nx.Graph()
    # labels = {}
    for coh in coh_name_list:
        for doc in doc_name_list:
            if df[coh][doc] == 1:
                G.add_node((coh, doc))
                # labels[(coh, doc)] = '%s<>%s' % (coh, doc)
    for c1, d1 in G.nodes_iter():
        for c2, d2 in G.nodes_iter():
            if df[c1][d2] == 0 and df[c2][d1] == 0:
                G.add_edge((c1, d1), (c2, d2))
    # pos = nx.spring_layout(G)
    # for node in labels:
    #     plt.annotate(labels[node], xy=pos[node])
    cliques = [a for a in nx.find_cliques(G)]
    max_len = max([len(a) for a in cliques])
    max_cliques = [a for a in cliques if len(a) == max_len]
    print len(max_cliques)
    clique_coh_list, clique_doc_list = coh_doc_set_span_maximal_cliques(max_cliques)
    print 'cohs that span entire clique list', clique_coh_list
    print 'docs that span entire clique list', clique_doc_list
    # best_ranker, best_rank = best_clique_by_overlapp(max_cliques, clique_coh_list, clique_doc_list)
    best_ranker, best_rank = best_clique_by_purples(max_cliques, df_true_score)
    print 'best ranker\n', best_ranker, best_rank
    ### find least similar clique:
    min_similarity = min(clique_similarity(best_ranker, a) for a in max_cliques)
    min_similars = []
    for clique in max_cliques:
        similarity = clique_similarity(best_ranker, clique)
        if similarity == min_similarity:
            min_similars.append(clique)
    best_min_similar_ranker, best_min_similar_rank = best_clique_by_purples(min_similars, df_true_score)
    print best_min_similar_ranker, best_min_similar_rank
    print min_similarity
    ### show true-score heat map for the best ranks clique:
    show_clique_heatmap(best_ranker, df_true_score, coh_name_list, doc_name_list)
    ### show true-score heat map for the least similar clique:
    show_clique_heatmap(best_min_similar_ranker, df_true_score, coh_name_list, doc_name_list)


def show_clique_heatmap(clique, df_tg, coh_name_list, doc_name_list):
    clique_cohs = [a[0] for a in clique]
    clique_docs = [a[1] for a in clique]
    df_tg = df_tg.drop([a for a in coh_name_list if a not in clique_cohs], axis=1)
    df_tg = df_tg.drop([a for a in doc_name_list if a not in clique_docs], axis=0)
    show_prediction_heat_map(df_tg)


def clique_similarity(c1, c2):
    '''
    :param c1:a clique
    :param c2:a clique
    :return:number of overlapping clique members
    '''
    cohs_c1 = [a[0] for a in c1]
    docs_c1 = [a[1] for a in c1]
    cohs_c2 = [a[0] for a in c2]
    docs_c2 = [a[1] for a in c2]
    score = 0
    for coh in cohs_c1:
        if coh in cohs_c2:
            score += 1
    for doc in docs_c1:
        if doc in docs_c2:
            score += 1
    return score


def best_clique_by_purples(max_cliques, df_tg):
    '''
    :param max_cliques:list of maximal cliques
    :param df_tg: data frame of true purple scores
    :return:the clique with the highest total purple count. < 15 is punished with -30
    '''
    best_rank = 0
    best_ranker = ''
    for clique in max_cliques:
        rank = 0
        for a in clique:
            coh = a[0]
            doc = a[1]
            if df_tg[coh][doc] < 15:
                rank -= 30
            rank += df_tg[coh][doc]
        if rank > best_rank:
            best_rank = rank
            best_ranker = clique
    return best_ranker, best_rank


def best_clique_by_overlapp(max_cliques, clique_coh_list, clique_doc_list):
    '''
    :param max_cliques: list of maximal cliques
    :param clique_coh_list: coh names
    :param clique_doc_list: doc names
    :return:clique with highest rank. rank being sum of its cph/doc members ranks
    '''
    coh_rank, doc_rank = rank_coh_or_doc(max_cliques, clique_coh_list, clique_doc_list)
    best_rank = 0
    best_ranker = ''
    for clique in max_cliques:
        rank = sum(coh_rank[a[0]] for a in clique) + sum(doc_rank[a[1]] for a in clique)
        if rank > best_rank:
            best_rank = rank
            best_ranker = clique
    return best_ranker, best_rank
    # print_names_for_cp(best_ranker)


def coh_doc_set_span_maximal_cliques(max_cliques):
    '''
    :param max_cliques:list of maximal cliques
    :return:dicts of ranks for cohs and docs, where rank is in how many cliqeus in max_cliques they are members
    '''
    clique_coh_list = set(b[0] for a in max_cliques for b in a)
    clique_doc_list = set(b[1] for a in max_cliques for b in a)
    return clique_coh_list, clique_doc_list


def rank_coh_or_doc(max_cliques, coh_list, doc_list):
    '''
    :param max_cliques:list of maximal cliques
    :param coh_list: list of cohs
    :param doc_list: list of docs
    :return:dicts for cohs and docs name:
    '''
    coh_res = {a: 0 for a in coh_list}
    doc_res = {a: 0 for a in doc_list}
    for clique in max_cliques:
        for pair in clique:
            coh_res[pair[0]] += 1
            doc_res[pair[1]] += 1
    return coh_res, doc_res


def print_names_for_cp(g):
    msg = ''
    for a in g:
        for b in g:
            msg += '*' + a[0] + '*' + b[1] + '*' + ' '
    print msg


def show_prediction_heat_map(df):
    import matplotlib.pyplot as plt
    from matplotlib import colors
    from numpy import array, arange
    cmap = colors.ListedColormap(['white', 'red', 'green', 'purple'])
    bounds = [-100, 0, 10, 15, 100]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    heatmap = plt.pcolor(array(df), cmap=cmap, norm=norm, edgecolors='k')
    for y in range(array(df.shape)[0]):
        for x in range(array(df.shape)[1]):
            if array(df)[y, x] >= 0:
                plt.text(x+0.5, y+0.5, array(df)[y, x], horizontalalignment='center', verticalalignment='center')
    plt.yticks(arange(0.5, len(df.index), 1), df.index)
    plt.xticks(arange(0.5, len(df.columns), 1), df.columns, rotation='vertical')
    plt.xlabel('cohesin name')
    plt.ylabel('dockerin name')
    plt.title('cohesin dockerin cross binding')
    plt.show()


if __name__ == '__main__':
    main()