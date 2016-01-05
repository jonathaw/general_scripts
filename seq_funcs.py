#!/usr/bin/env python3.5
"""
a collection of sequence functions
"""
import pandas as pd
from AASeq import AASeq


def read_multi_fastas(fastas_file: str, suffix_to_remove: str=None, lower=False, add_aligned=False) -> dict:
    """
    :param fastas_file: file address
    :return: {name: AASeq}
    """
    with open(fastas_file, 'r') as f:
        cont = f.read().split('>')
    result = {}
    for entry in cont:
        split_entry = entry.split('\n')
        if len(split_entry) < 2:
            continue
        name = '_'.join(split_entry[0].rstrip().split())
        if name == '':
            continue
        if suffix_to_remove is not None:
            name = name.split(suffix_to_remove)[0]
        seq = ''.join(a.rstrip() for a in split_entry[1:])
        if '-' in seq or add_aligned:
            aln = seq
            seq = aln.replace('-', '')
            if lower:
                result[name.lower()] = AASeq(string=seq, name=name.lower(), aligned=aln)
            else:
                result[name] = AASeq(string=seq, name=name, aligned=aln)
        else:
            if lower:
                result[name.lower()] = AASeq(string=seq, name=name.lower())
            else:
                result[name] = AASeq(string=seq, name=name)
    return result


def parse_msa_to_df(file_name: str) -> pd.DataFrame:
    """
    :param file_name: file name
    :return: data frame with seq name as index, and column for every aligned position
    """
    seq_dict = read_multi_fastas(file_name)
    first_seq = list(seq_dict.values())[0]
    df_ = pd.DataFrame(columns=range(1, len(first_seq.get_aligned())+1), index=list(seq_dict.keys()))
    for seq in seq_dict.values():
        df_.loc[seq.name] = list(seq.get_aligned())
    return df_


def pair_wise_aln_from_seqs(seq1, seq2, matrix=None, gap_open=-10, gap_extend=-0.5):
    """
    :param seq1: a fasta sequence
    :param seq2: a fasta sequence
    :param matrix: a substitution matrix
    :param gap_open: gap_open penalty
    :param gap_extend: gap_extend penalty
    :return: aln for seq1, aln for seq2, score, beginning and end for the best alignment
    """
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
    if matrix is None:
        matrix = matlist.blosum62
    alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)[0]
    return alns


def aln_identity(aln1: AASeq, aln2: AASeq) -> float:
    """
    :param aln1: alignment sequence (with gaps)
    :param aln2: alignment sequence (with gaps)
    :return: the identity calculated by: (# identities)/(aln1 length, no gaps)
    >>> a = AASeq(aligned='-ABC-D')
    >>> b = AASeq(aligned='-ABCED')
    >>> aln_identity(a, b)
    1.0
    >>> b = AASeq(string='-BBCED')
    >>> aln_identity(a, b)
    0.75
    """
    res = 0.
    for i, aa in enumerate(aln1):
        res += 1. if aln1.get_aligned()[i] == aln2.get_aligned()[i] != '-' else 0.
    length = float(len(len(aln1)))
    return res/length


def write_multi_seqs_to_file(seqs: dict, out_file: str, query: AASeq=None, no_dups: bool=True):
    """
    :param seqs: {name: AASeq}
    :param out_file: file to write
    :param query: query AASeq
    :param no_dups: whether to refrain from duplicates or not
    :return: write fasta file
    """
    written_seqs = []
    with open(out_file, 'w+') as fout:
        if query is not None:
            fout.write('%s\n' % query.write())
            written_seqs.append(query)
        for s in seqs.values():
            if query is not None:
                if query == s:
                    continue
            if no_dups:
                if s in written_seqs:
                    continue
            fout.write('%s\n' % s.write())
            written_seqs.append(s)


def write_seq_from_seq(name, seq, seq_out):
    with open(seq_out, 'wr+') as fout:
        fout.write('>%s\n%s' % (name, seq))


def run_muscle(fastas_file, muscle_out):
    from os import system
    system('muscle -in %s -out %s' % (fastas_file, muscle_out))


def run_psi_blast_for_pssm(seq_file, msa_file, pssm_file):
    import os
    os.system('psiblast -subject %s -in_msa %s -out_ascii_pssm %s' % (seq_file, msa_file, pssm_file))


def arrange_for_pssm(msa_file, out_file, query):
    from shutil import move
    fastas = read_multi_fastas(msa_file)
    write_multi_seqs_to_file(fastas, msa_file+'_tmp', query)
    move(msa_file+'_tmp', msa_file)
    write_seq_from_seq(query, fastas[query]['seq'], query+'.fasta')
    print('running pssm')
    run_psi_blast_for_pssm(query+'.fasta', msa_file, out_file)


if __name__ == '__main__':
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument('-name', type=str)
    parser.add_argument('-path', type=str, default=os.getcwd()+'/')
    parser.add_argument('-mode', type=str)
    parser.add_argument('-seq1', type=str)
    parser.add_argument('-seq2', type=str)
    parser.add_argument('-aln_matrix', type=str)
    parser.add_argument('-gap_open', type=int)
    parser.add_argument('-gap_extend', type=int)
    parser.add_argument('-fastas_file', type=str)
    parser.add_argument('-muscle_out', type=str)
    parser.add_argument('-seq_out', type=str)
    parser.add_argument('-pssm_file', type=str)
    parser.add_argument('-msa_file', type=str)
    parser.add_argument('-seq_file', type=str)
    args = vars(parser.parse_args())

    if args['mode'] == 'test':
        parse_msa_to_df('/home/labs/fleishman/jonathaw/decision_tree/cohesins_from_rachel_and_vered.fasta_aln')