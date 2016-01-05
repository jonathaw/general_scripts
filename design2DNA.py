#!/usr/bin/env python3.5
import argparse
import re
from design2bins_by_posistions import parse_name_translation
from seq_funcs import read_multi_fastas

genetic_code = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "s", "TCA": "S", "TCG": "S", "TAT":
                "Y", "TAC": "Y", "TAA": "STOP", "TAG": "STOP", "TGT": "C", "TGC": "C", "TGA": "STOP", "TGG": "W",
                "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "ATT":
                "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAT": "N",
                "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GTT": "V",
                "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "GAT": "D", "GAC":
                "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G", }


def add_flanks(args):
    doc_flanks = {'1ohz': {'start': 'ESSSVLL', 'end': 'RVIDKFPVAENP'},
                  '2vn5': {'start': 'V', 'end': 'SKLPSN'},
                  '3ul4': {'start': 'V', 'end': ''},
                  '4dh2': {'start': 'WNK', 'end': 'NSAPTF'},
                  '5new': {'start': '', 'end': 'Y'},
                  }
    name_path = '/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/reclique_18Nov/stabilisation/'
    translation = parse_name_translation(name_path + 'translate_names.txt')

    if args['type'] == 'coh':
        return 'D' + args['seq'] + 'NAT'
    elif args['type'] == 'doc':
        doc_name = translation[args['name'] + '.pdb.gz'][10:14]
        return doc_flanks[doc_name]['start'] + args['seq'] + doc_flanks[doc_name]['end']


def DNA2AA(dna: str) -> str:
    """
    :param dna: a dna string 
    :return: AA seq
    >>> DNA2AA('TTTACT')
    'FT'
    """
    assert len(dna) % 3 == 0
    res = ''
    for j in range(0, len(dna) - 3 + 1, 3):
        res += genetic_code[dna[j:j + 3]]
    return res


def validate(args):
    # original_seqs = read_multi_fastas(args['original_seqs_file'], suffix_to_remove='_')
    DNA_seqs = read_multi_fastas(args['DNA_seqs_file'], suffix_to_remove='.')

    for k, v in DNA_seqs.items():
        # assert original_seqs[k].get_seq() in DNA2AA(v.get_seq())
        if not gen9_standards(v.get_seq):
            print('seq name %s does not comply with Gen9 standards' % k)


def check_homopolymer(dna: str, nuc: str, threshold: int) -> bool:
    reg = re.compile('[%s%s]*' % (nuc.lower(), nuc.upper()))
    for AA in reg.findall(dna):
        if len(AA) >= threshold:
            print('too many %ss' % nuc)
            return False
    return True


def check_internal_repeats(dna):
    forty_mers = []
    for i in range(0, len(dna) - 40 + 1):
        if dna[i:i + 41] in forty_mers:
            print('repeating seq too long %s' % dna[i:i + 41])
            return False
    return True


def calc_gc_content(sli: str) -> float:
    return float(sli.count('G') + sli.count('C')) / float(len(sli))


def check_100_gc(dna: str) -> bool:
    for i in range(len(dna) - 100 + 1):
        if not 0.3 <= calc_gc_content(dna[i:i + 101]) <= 0.7:
            print('slice GC content out of range %s' % dna[i:i + 101])
            return False
    return True


def reverse_complement(seq: str) -> str:
    """
    :param seq: dna sequence
    :return: the reverse complement
    >>> reverse_complement('AAATTTGGGCCC')
    'GGGCCCAAATTT'
    """
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join([comp[a] for a in seq[::-1]])


def gen9_standards(dna):
    homomers = {'A': 8, 'C': 8, 'G': 5, 'T': 8}
    for nuc, thre in homomers.items():
        res = check_homopolymer(dna, nuc, thre)
        if not res:
            print('failed homopolymers %s' % nuc)
            return False

    if 'GGTCTC' in dna or reverse_complement('GGTCTC') in dna:
        print('unallowed restriction sites GGTCTC')
        return False

    if 'CACCTGC' in dna or reverse_complement('CACCTGC') in dna:
        print('unallowed restriction sites CACCTGC')
        return False

    if not check_internal_repeats(dna):
        print('failed internal repeats')
        return False

    if not 0.4 <= calc_gc_content(dna) <= 0.65:
        print('overall GC content not in range, %f' % calc_gc_content(dna))
        return False

    if not check_100_gc(dna):
        print('failed 100mers GC contents')
        return False

    return True


def add_primers(seq_: str, t_: str) -> str:
    primer_flanks = {'coh': {'start': 'CGTCAGATGATCCGAATGCAGGATCC', 'end': 'TAACTCGAGCACCACCACCACCAC'},
                     'doc': {'start': 'TGGGCTATTATCGACCACAAAGTGGTACCA', 'end': 'TAAGGATCCGGCTGCTAACAAAGCCCG'}}
    return primer_flanks[t_]['start'] + seq_ + primer_flanks[t_]['end']


def add_primers_to_all(args):
    DNA_seqs = read_multi_fastas(args['DNA_seqs_file'])
    for k, v in DNA_seqs.items():
        print('>%s' % k)
        print(add_primers(v.get_seq, args['type']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode')
    parser.add_argument('-seq')
    parser.add_argument('-type')
    parser.add_argument('-name')
    parser.add_argument('-original_seqs_file')
    parser.add_argument('-DNA_seqs_file')
    args = vars(parser.parse_args())

    if args['mode'] == 'seq2flanks':
        print(add_flanks(args))

    elif args['mode'] == 'add_primers':
        add_primers_to_all(args)

    elif args['mode'] == 'validate':
        validate(args)

    elif args['mode'] == 'test':
        print(reverse_complement('AAATTTGGGCCC'))

    else:
        print('no mode chosen')
