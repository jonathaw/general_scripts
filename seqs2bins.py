def main():
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
    import sys
    from rosetta_score_files import score2dict
    import operator
    from collections import Counter
    ### this positions dict is for the 1st, 8 parts switches design
    # positions_dict = {'1anu': [36, 38, 114, 115, 117, 120, 124, 126], '1ohz': [37, 39, 115, 116, 118, 121, 125, 127],
    #                   '2ccl': [37, 39, 115, 116, 118, 121, 125, 127]}
    ### this positions dict if for the second, 10 parts switches design from 1-2.3.2015
    positions_dict = {'1anu': [32, 36, 62, 65, 69, 82, 115, 126],
                      '1aoh': [33, 37, 63, 66, 70, 83, 119, 130],
                      '1ohz': [33, 37, 63, 66, 70, 83, 116, 127],
                      '2ccl': [33, 37, 63, 66, 70, 83, 116, 127]}
    type_dict = {'D': 'n', 'E': 'n', 'K': 'p', 'R': 'p'}
    scores = score2dict(sys.argv[2])
    bins = {}
    num_structs_bin = {}
    f = open(sys.argv[1], 'r')
    cont = f.read().split('>')
    seq_dict = {i.split('\n')[0]: i.split('\n')[1] for i in cont if len(i) > 0}
    for name, seq in seq_dict.items():
        if scores[name]['purple']:
            coh_name = name.split('_')[0]
            switches = ''.join([type_dict[seq[i-1]] if seq[i-1] in type_dict.keys()
                                else 'c' for i in positions_dict[coh_name]])
            ### adding a condition where total #charges is <= 7, and distributes 5/2 or 4/3:
            counter = Counter(switches)
            if counter['n']+counter['p'] != 7 or counter['n'] < 2 or counter['p'] < 2:
                continue
            ###
            if switches not in bins.keys():
                bins[switches] = []
                num_structs_bin.update({switches: 0})
            bins[switches].append(name)
            num_structs_bin[switches] += 1

    bin_subsets = bin_subsetter(bins.keys())
    bin_subsets_struct = {i: [num_structs_bin[j] for j in bini] for i,bini in enumerate(bin_subsets)}
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
        print biner.keys()[0].upper()
        for struct in biner.values()[0]:
            print struct
        print '\n'


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
            print '%s != %s' % (query, seq)
            return True
    return False

if __name__ == '__main__':
    main()