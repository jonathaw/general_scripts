#!/usr/bin/env python3.5
import math
import argparse
from collections import OrderedDict

from List1 import List1
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq


class AASeq():
    def __init__(self, string: str=None, name: str=None, aligned: str=None):
        self.seq = '0'+string if string is not None else None
        self.aligned = aligned if aligned is None else '0'+aligned
        self.name = name

    def __str__(self) -> str:
        msg = '>%s: %s' % (self.name, self.seq[1:])
        if self.aligned is not None:
            msg += ' aligned: %s' % self.aligned[1:]
        return msg

    def __repr__(self) -> str:
        return self.__str__()

    def __getitem__(self, item: int) -> None:
        """
        >>> seq = AASeq(name='john', string='ACDEFGHIKLMNPQRSTVWY')
        >>> seq[1]
        'A'
        >>> seq[-1]
        'Y'
        >>> seq[3]
        'D'
        >>> seq[5]
        'F'
        """
        return self.seq[item]

    def __setitem__(self, key: int, value: str) -> None:
        """
        >>> seq = AASeq(name='john', string='AAA')
        >>> seq[1] = 'B'
        >>> seq.seq
        '0BAA'
        """
        self.seq = self.seq[:key] + value + self.seq[key+1:]

    def __eq__(self, other) -> bool:
        return self.seq == other.seq and self.name == other.name

    def __len__(self) -> int:
        """
        >>> seq = AASeq(name='john', string='AAA')
        >>> len(seq)
        3
        """
        return len(self.get_seq())

    def __iter__(self):
        return self.seq[1:].__iter__()

    def __contains__(self, item: str) -> bool:
        return self.seq[1:].__contains__(item)

    def add_suffix(self, suffix: str) -> None:
        self.seq += suffix

    def add_prefix(self, prefix: str) -> None:
        self.seq = '0' + prefix + self.get_seq()

    def index(self, aa: str) -> list:
        """
        :param aa: an AA
        :return: list of indices of aa in self
        >>> s = AASeq(string='ABCDABCDAAA')
        >>> s.index('A')
        [1, 5, 9, 10, 11]
        """
        return [i+1 for i, a in enumerate(self) if a == aa]

    def write(self) -> str:
        return ">%s\n%s" % (self.name, self.get_seq())

    def split(self, d: str) -> list:
        return self.get_seq().split(d)

    def enumerate(self) -> (int, str):
        for i, aa in enumerate(self.seq[1:]):
            yield i+1, aa

    def enumerate_aligned(self) -> (int, str):
        for i, aa in enumerate(self.aligned[1:]):
            yield i+1, aa

    def add_aa(self, other: str):
        self.seq += other

    def set_name(self, name: str) -> None:
        self.name = name

    def set_seq(self, seq: str) -> None:
        self.seq = '0'+seq

    def get_seq(self) -> str:
        """
        :return: the sequence as a string
        >>> s = AASeq(string='ABCDEFG')
        >>> s.get_seq()
        'ABCDEFG'
        """
        return self.seq[1:]

    def get_aligned(self) -> str:
        if self.aligned is None:
            return self.get_seq()
        return self.aligned[1:]

    def write_file(self, file_name: str=None) -> None:
        if file_name is None:
            import os
            file_name = '%s/%s.fasta' % (os.getcwd(), self.name)
        with open(file_name, 'w+') as fout:
            fout.write('>%s\n%s\n' % (self.name, self.get_seq()))

    def align(self, other, matrix=None, gap_open=-10,
              gap_extend=-0.5) -> (float, int, int):
        """
        :param other: another AASeq
        :param matrix: a matrix
        :return: score, begin, end. updates self and other with .aligned
        """
        from Bio import pairwise2
        from Bio.SubsMat import MatrixInfo as matlist
        if matrix is None:
            matrix = matlist.blosum62
        alns = pairwise2.align.globalds(self.get_seq(), other.get_seq(), matrix,
                                        gap_open, gap_extend)[0]
        aln1, aln2, score, begin, end = alns
        self.aligned = '0' + aln1
        other.aligned = '0' + aln2
        return score, begin+1, end

    def non_aligned_position_at_aligned(self, pos: int) -> int:
        """
        :param pos: a position in the non-aligned sequence (int)
        :return: the respective position in the aligned version
        >>> s = AASeq(string='ABCDEFGHIJKLMNOP', aligned='ABC-DEFG-HIJK-LMNO-P')
        >>> s.non_aligned_position_at_aligned(1)
        1
        >>> s.non_aligned_position_at_aligned(4)
        5
        >>> s.non_aligned_position_at_aligned(16)
        20
        """
        assert self.aligned is not None, "there is no aligned sequence for AASeq instance %s" % str(self)
        no_gaps = 0
        for i, aa in enumerate(self.aligned[1:]):
            no_gaps += 1 if aa != '-' else 0
            if no_gaps == pos:
                return i+1

    def aligned_position_at_non_aligned(self, pos: int) -> int:
        """
        :param pos: a position at the aligned sequence
        :return: the respective position at the non_aligned version
        >>> s = AASeq(string='ABCDEFGHIJKLMNOP', aligned='ABC-DEFG-HIJK-LMNO-P')
        >>> s.aligned_position_at_non_aligned(1)
        1
        >>> s.aligned_position_at_non_aligned(5)
        4
        >>> s.aligned_position_at_non_aligned(10)
        8
        >>> s.aligned_position_at_non_aligned(20)
        16
        """
        return len(self.aligned[:pos].replace('-', ''))

    def aligned_identity(self, other) -> float:
        """
        :param other: another AASeq
        :return: identity percentage,
        calculated by number of identical positions in alignment,
        devided by length of self.seq
        >>> a = AASeq(string='ABCDEFGTTT', aligned='ABCDEFGTTT')
        >>> b = AASeq(string='ABCFGRRR',   aligned='ABC--FGRRR')
        >>> a.aligned_identity(b)
        0.5
        """
        identical = 0
        for i, aa in self.enumerate_aligned():
            if aa == '-':
                continue
            identical += 1 if aa == other.aligned[i] else 0
        return float(identical) / float(len(self))

    def get_positions(self, pos_list: list) -> list:
        """
        :param pos_list: a list of positions in 1-num index
        :return: a list of residues ath the requested positions
        >>> a = AASeq(string='ABCDEFG')
        >>> a.get_positions([1, 4, 6])
        ['A', 'D', 'F']
        """
        return [self[int(pos)] for pos in pos_list]

    def get_aligned_positions(self, pos_list: list) -> list:
        """
        :param pos_list: list of integers positions
        :return: list of identities in the aligned version
        >>> a = AASeq(string='ACDEF', aligned='ACD-EF')
        >>> a.get_aligned_positions([2, 4, 6])
        ['C', '-', 'F']
        """
        return [self.aligned[int(pos)] for pos in pos_list]

    def non_aligned_identity(self, other) -> float:
        """
        :param other: another AASeq
        :return: identity percentage, calculated by number of identical
        positions in non_alignment,
        devided by length of self.seq
        >>> a = AASeq(string='ABCDEFGTTT')
        >>> b = AASeq(string='ABCDTTGTTT')
        >>> a.non_aligned_identity(b)
        0.8
        """
        identical = 0
        for i, aa in self.enumerate():
            identical += 1 if aa == other[i] else 0
        return float(identical) / float(len(self))

    def calc_molecular_weight(self) -> float:
        """
        :return: protein seq molecular weight, float
        """
        return molecular_weight(self.get_seq(), seq_type='protein')

    def calc_extinction_coefficient(self, reduced=True) -> float:
        """
        data and equation from web.expasy.org/protparam/protparam-doc.html
        :param reduced: whether to calculate cysteins as SS bonds or not. if
        True then cys contribution is 0
        :return: computed extinction coefficient
        """
        # AA extinction coefficients:
        ext_tyr = 1490.
        ext_trp = 5500.
        ext_cystine = 125.

        # number of residues:
        num_tyr = self.seq.count('Y')
        num_trp = self.seq.count('W')
        num_cys = math.floor(self.seq.count('C') / 2) if not reduced else 0

        return ext_tyr*num_tyr + ext_trp*num_trp + ext_cystine*num_cys

    def calc_isoelectric_point(self) -> float:
        """
        using biopython.org/DIST/docs/api/Bio.SeqUtils.ProtParam-pysrc.html
        :return: calculates the sequence's isoelectric point
        """
        protein_analysis = ProteinAnalysis(self.get_seq())
        return protein_analysis.isoelectric_point()

    def aa_frequency(self, aa: str, start: int=1, end: int=-1) -> float:
        """
        >>> a= AASeq(string='ABCDEFGHIJ')
        >>> a.aa_frequency('A')
        0.1
        >>> a.aa_frequency('A', 1)
        0.1
        >>> a.aa_frequency('A', 2, 5)
        0.0
        >>> a.aa_frequency('J', 9)
        0.5
        """
        if end != -1:
            subseq = self.get_seq()[start-1:end+1]
        else:
            subseq = self.get_seq()[start-1:]
        aa_num = float(subseq.count(aa))
        return aa_num / float(len(subseq))

    def all_aas_frequencies(self, start: int=1, end: int=-1,
                            clean_zeros: bool=False) -> OrderedDict:
        d = {}
        for aa in list('ACDEFGHIKLMNPQRSTVWY'):
            v = self.aa_frequency(aa, start, end)
            if not clean_zeros:
                d[aa] = v
            elif clean_zeros and v != 0:
                d[aa] = v
        return OrderedDict(sorted(d.items(), key=lambda t: t[1]))


def compare_2_seqs(seq_1: AASeq, seq_2: AASeq, start=0) -> None:
    """
    :param seq_1: an AASeq instance
    :param seq_2: an AASeq instance
    :return: prints a comparison
    """
    print('comparing %s to %s' % (seq_1.name, seq_2.name))
    changes = []
    for i, aa in seq_1.enumerate():
        if aa != seq_2[i]:
            print('%s%i%s' % (aa, i, seq_2[i]))
            changes.append(i)
    print('found %i changes, thats is %.2f%% differences' %
          (len(changes), 100*len(changes)/len(seq_1)))
    print('select changes, resi %s' % '+'.join([str(change+start)
                                                for change in changes]))


def read_seq(file_name: str) -> AASeq:
    with open(file_name, 'r') as fin:
        cont = fin.read().split('\n')
        name = cont[0][1:]
        seq = ''.join(cont[1:])
    return AASeq(string=seq, name=name)


def read_seqs(file_name: str, remove_suffix=None) -> OrderedDict:
    result = OrderedDict()
    with open(file_name, 'r') as fin:
        cont = fin.read().split('>')
        for p in cont[1:]:
            s = p.split('\n')
            if remove_suffix is None:
                name = s[0]
            else:
                name = s[0].split(remove_suffix)[0]
            result.update({name: AASeq(string=s[1], name=s[0])})
    return result


if  __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-seq1')
    parser.add_argument('-seq2')
    parser.add_argument('-mode', default='compare')
    args = vars(parser.parse_args())

    if args['mode'] == 'compare':
        s1 = AASeq(string=args['seq1'])
        s2 = AASeq(string=args['seq2'])
        compare_2_seqs(s1, s2)

    # a = AASeq(string='ABCDEFG')
    # b = AASeq(string='ABCFG')
    # a.align(b)
    # print(a)
    # print(b)
    # print(a.aligned_identity(b))
    # c = AASeq(string='ABCDEG')
    # a.align(c)
    # print(a)
    # print(a.aligned_identity(c))
