#!/usr/bin/env python3.5
"""
script for analysing membrane protein designs
"""
import os
import argparse
import math
import matplotlib.pyplot as plt
from collections import OrderedDict

from AASeq import read_seqs
import RosettaFilter as Rf


wt_aa_freq = OrderedDict([('L', 0.16), ('I', 0.1), ('A', 0.09), ('F', 0.08),
                          ('G', 0.08), ('V', 0.08), ('S', 0.07), ('T', 0.07),
                          ('M', 0.04), ('P', 0.04), ('Y', 0.03), ('W', 0.03),
                          ('N', 0.03), ('H', 0.02), ('C', 0.02), ('Q', 0.01),
                          ('E', 0.01), ('R', 0.01), ('K', 0.01), ('D', 0.01)])
aas = list('ACDEFGHIKLMNPQRSTVWY')

# color_map = {'A': 'yellowgreen', 'C': 'gold', 'D': 'lightskyblue',
             # 'E': 'lightcoral', 'F': 'green', 'G': 'red', 'H': 'tomato',
             # 'I': 'orange', 'K': 'lightblue', 'L': 'white', 'M': 'yellow',
             # 'N': 'lightsalmon', 'P': 'violet', 'Q': 'lime',
             # 'R': 'royalblue', 'S': 'tan', 'T': 'beige', 'V': 'grey',
             # 'W': 'plum', 'Y': 'orchid'}

color_map = {'A': 'palegoldenrod', 'C': 'gold', 'D': 'indianred',
             'E': 'firebrick', 'F': 'burlywood', 'G': 'red', 'H': 'deepskyblue',
             'I': 'blanchedalmond', 'K': 'steelblue', 'L': 'white',
             'M': 'azure', 'N': 'darkgreen', 'P': 'maroon', 'Q': 'olivedrab',
             'R': 'navy', 'S': 'forestgreen', 'T': 'lightgreen', 'V': 'snow',
             'W': 'beige', 'Y': 'floralwhite'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', default='compare_designs')
    parser.add_argument('-sc')
    parser.add_argument('-fasta')
    parser.add_argument('-native')
    parser.add_argument('-show', default='show')
    parser.add_argument('-dir')
    args = vars(parser.parse_args())

    if args['mode'] == 'compare_designs':
        compare_designs(args)

    elif args['mode'] == 'meta_pie':
        meta_pie(args)

    else:
        print('no mode')


def calc_freq_dist(f1: OrderedDict, f2: OrderedDict, aas=None) -> float:
    if aas is None:
        to_test = list('ACDEFGHIKLMNPQRSTVWY')
    else:
        to_test = aas if isinstance(aas, list) else list(aas)
    dists = 0
    for aa in to_test:
        if aa in f1.keys() and aa in f2.keys():
            dists += (f1[aa] - f2[aa]) ** 2
    return math.sqrt(dists) / float(len(to_test))


def meta_pie(args: dict):
    mut_aa_freqs = {aa: 0 for aa in aas}
    score_files = [a for a in os.listdir(args['dir']) if '.score' in a]
    for fasta_file in [a for a in os.listdir(args['dir']) if '.fasta' in a]:
        temp_seqs = read_seqs('%s/%s' % (args['dir'], fasta_file), remove_suffix='.pdb')
        temp_sc = Rf.score_file2df('%s/%s.score' %
                                   (args['dir'], fasta_file.split('.')[0]))
        names = list(Rf.get_best_num_by_term(temp_sc, 10, 'a_tms_aa_comp')['description'])
        for n, aaseq in temp_seqs.items():
            if n in names:
                for aa in aas:
                    mut_aa_freqs[aa] += (aaseq.aa_frequency(aa) * len(aaseq))
    mut_aa_freqs_srt = OrderedDict(sorted(mut_aa_freqs.items(),
                                          key=lambda t: t[1]))

    ori_aa_freqs = {aa: 0 for aa in aas}
    for aaseq in read_seqs('/home/labs/fleishman/jonathaw/elazaridis/design/' +
                           'polyA_13Nov/chosen_from_all_27Feb/pdbs/' +
                           'all_dzns.fasta').values():
        for aa in aas:
            ori_aa_freqs[aa] += (aaseq.aa_frequency(aa) * len(aaseq))
    ori_aa_freqs_srt = OrderedDict(sorted(ori_aa_freqs.items(),
                                          key=lambda t: t[1]))

    plt.figure()
    plt.subplot(1, 3, 1)
    plt.title('natural TMs')
    plt.pie(list(wt_aa_freq.values()), labels=list(wt_aa_freq.keys()),
            autopct='%1.1f%%',
            colors=[color_map[a] for a in list(wt_aa_freq.keys())])
    plt.axis('equal')

    plt.subplot(1, 3, 2)
    plt.title('original designs')
    plt.pie(list(ori_aa_freqs_srt.values()),
            labels=list(ori_aa_freqs_srt.keys()),
            autopct='%1.1f%%',
            colors=[color_map[a] for a in list(ori_aa_freqs_srt.keys())])
    plt.axis('equal')

    plt.subplot(1, 3, 3)
    plt.title('mutated designs')
    plt.pie(list(mut_aa_freqs_srt.values()),
            labels=list(mut_aa_freqs_srt.keys()),
            autopct='%1.1f%%',
            colors=[color_map[a] for a in list(mut_aa_freqs_srt.keys())])
    plt.axis('equal')

    plt.show()


def compare_designs(args):
    all_seqs = read_seqs(args['fasta'], '.pdb')
    sc_df = Rf.score_file2df(args['sc'])
    best_by_aa_freq_df = Rf.get_best_num_by_term(sc_df, 10, 'a_tms_aa_comp')
    names_to_use = list(best_by_aa_freq_df['description'])
    names_to_use += [args['sc'].split('.')[0].split('all_')[1]]
    all_seqs = {k: v for k, v in all_seqs.items() if k in names_to_use}

    seqs_aa_freqs = {n: {} for n in all_seqs.keys() if n in names_to_use}
    freq_dists = {}
    for n, s in all_seqs.items():
        seqs_aa_freqs[n] = s.all_aas_frequencies(clean_zeros=True)
        # freq_dists[n] = calc_freq_dist(seqs_aa_freqs[n], wt_aa_freq,
                                        # 'AGFILPSTVWY')
        freq_dists[n] = calc_freq_dist(seqs_aa_freqs[n], wt_aa_freq)

    freq_dists = OrderedDict(sorted(freq_dists.items(), key=lambda t: t[1]))
    fig = plt.figure(figsize=(18, 18))
    nrows = math.floor(len(all_seqs.keys()) / 4)
    nrows += 1 if len(all_seqs.keys()) % 4 > 0 else 0
    i = 0
    # for n, s in all_seqs.items():
    first = True
    for n, dist_freq in freq_dists.items():
        if first:
            print(n)
            first = False
        s = all_seqs[n]
        plt.subplot(nrows, 4, 1+i)
        fracs = [a for a in seqs_aa_freqs[n].values()]
        labels = [a for a in seqs_aa_freqs[n].keys()]
        colors = [color_map[a] for a in labels]
        plt.pie(fracs, labels=labels, autopct='%1.1f%%', colors=colors)
        name = n.split('.pdb')[0]
        if name.count('poly') == 2:
            p_name = name.split('_poly')[0]
        else:
            p_name = name
        if name in list(sc_df['description']):
            title = ('%s\nscore: %.0f ∆∆G: %.0f\nfreq_dist: %.2f\naa_comp %.2f' %
                     (p_name, sc_df[sc_df['description'] == name]['score'],
                      sc_df[sc_df['description'] == name]['a_ddg'],
                      dist_freq*100,
                      best_by_aa_freq_df[best_by_aa_freq_df['description'] == name]['a_tms_aa_comp']))
        else:
            title = 'freq_dist %.2f' % (dist_freq*100)
        plt.title(title)

        i += 1

    if args['show'] == 'show':
        plt.show()
    else:
        plt.savefig('%s_aa_freq.png' % args['fasta'])


if __name__ == '__main__':
    main()
