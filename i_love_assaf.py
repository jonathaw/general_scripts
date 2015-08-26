"""
A script to find a group of nucleotide sequences that translate to the same AA sequence (AA_SEQ), and have
HAMMING_DIST distance between them
"""


def main():
    import networkx as nx
    global GENETIC_CODE, AA_SEQ, AA_LENGTH, NUC_LENGTH, HAMMING_DIST
    GENETIC_CODE = {'L': ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                    'S': ['TCT', 'TCC', 'TCA', 'TCG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG']}
    AA_SEQ = 'RLRLRLRLRL'
    # AA_SEQ = 'RRRR'
    HAMMING_DIST = 8
    AA_LENGTH = len(AA_SEQ)
    NUC_LENGTH = 3*AA_LENGTH
    current_opts = [nucs for nucs in GENETIC_CODE[AA_SEQ[0]]]
    for aa in range(AA_LENGTH-1):
        new_opts = []
        for opt in current_opts:
            for nucs in GENETIC_CODE[AA_SEQ[len(current_opts[0])/3-1]]:
                new_opts.append(opt+nucs)
        current_opts = new_opts[:]

    print 'there are %i options' % len(current_opts)

    G = nx.Graph()
    [G.add_node(a) for a in current_opts]
    [G.add_edge(a, b) for a in current_opts for b in current_opts if different_enough(a, b)]

    print 'found %i edges' % len(G.edges())

    cliques = [a for a in nx.find_cliques(G)]

    print 'now printing cliques larger than 4'
    for clique in cliques:
        if len(clique) >= 10000:
            print clique
    # print cliques
    # max_len = max([len(a) for a in cliques])
    # max_cliques = [a for a in cliques if len(a) == max_len]
    # print len(max_cliques)
    # print max_cliques


def different_enough(x, y):
    diff = 0
    for i, n in enumerate(x):
        diff += 1 if n != y[i] else 0
    return diff >= HAMMING_DIST


if __name__ =='__main__':
    main()
