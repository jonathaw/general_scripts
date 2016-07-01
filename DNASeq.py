#!/usr/bin/env python3.5
import re
from AASeq import AASeq
genetic_code = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "TAT":
                "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", 
                "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", 
                "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "ATT": 
                "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAT": "N", 
                "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GTT": "V", 
                "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "GAT": "D", "GAC": 
                "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G", }
complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


class DNASeq:
    def __init__(self, seq: str=None, name: str=None):
        self.seq = '0' + seq
        self.name = name

    def __repr__(self) -> str:
        return '>%s\n%s' % (self.name, self.get_seq)

    def __getitem__(self, item) -> str:
        """
        :param item: position
        :return: nucleotide at position
        >>> dna_seq = DNASeq('ATGC')
        >>> dna_seq[1]
        'A'
        """
        return self.seq[item]

    def __setitem__(self, key, value) -> None:
        """
        :param key: position
        :param value: new for position
        :return: None
        >>> dna_seq = DNASeq('ATGC')
        >>> dna_seq[2] = 'G'
        >>> dna_seq == DNASeq('AGGC')
        True
        """
        self.seq = self.seq[:key] + value + self.seq[key+1:]

    def __eq__(self, other) -> bool:
        return self.name == other.name and self.seq == other.seq

    @property
    def __len__(self) -> int:
        """
        :return: length of sequence
         >>> len(DNASeq('ATGC'))
         4
        """
        return len(self.get_seq)

    @property
    def enumerate(self) -> (int, str):
        for i, nuc in enumerate(self.seq):
            if nuc != '0':
                yield (i, nuc)

    def get_seq(self) -> str:
        """
        :return: sequence as string
         >>> DNASeq('ATGC').get_seq
         'ATGC'
        """
        return self.seq[1:]

    @property
    def get_GC(self) -> float:
        """
        :return: the GC content in float
         >>> DNASeq('ATGCATGCAA').get_GC
         0.4
        """
        return float(self.get_seq.count('G') + self.get_seq.count('C')) / float(len(self))

    def has_seq(self, seq: str) -> bool:
        """
        :param seq: a sub sequence
        :return: if seq is in DNASeq
        >>> DNASeq('ATGCATGC').has_seq('CAT')
        True
        """
        return seq in self.get_seq

    def find_seq(self, seq: str) -> list:
        """
        :param seq: sequence
        :return: a list of where sequence occurs
        >>> DNASeq('ATGCATGCATGC').find_seq('AT')
        [1, 5, 9]
        """
        return [m.start() for m in re.finditer(seq, self.seq)]

    @property
    def reverse_complement(self) -> str:
        """
        :return: reverse complement of the sequence
         >>> DNASeq('ATGC').reverse_complement
         'GCAT'
        """
        return ''.join([complement[a] for a in self.get_seq[::-1]])

    def reverese(self) -> str:
        """
        :return: reverse with 0 in the beginning
         >>> dnaseq = DNASeq('ATGC')
         >>> dnaseq.reverese()
         '0CGTA'
        """
        return '0' + self[::-1][:-1]

    def split(self, d: str) -> list:
        return self.get_seq().split(d)


def translate(seq: str, name=None) -> AASeq:
    """
    :param seq: a nucleotide seq
    :return: amino acid seq
    >>> translate('TTTCATAAG').get_seq()
    'FHK'
    """
    return AASeq(string=''.join([genetic_code[seq[i:i+3]] for i in range(0, len(seq)-3+1, 3)]), name=name)


def find_protein_coding_regions(dna_seq: DNASeq) -> dict:
    result = {}
    names = {1: 'F1', 2: 'F2', 3: 'F3', -1: 'R1', -2: 'R2', -3: 'R3'}
    rev_seq = dna_seq.reverese()
    for start in [1, 2, 3]:
        result[names[start]] = translate(dna_seq[abs(start):], name=names[start])
        result[names[-start]] = translate(rev_seq[abs(start):], name=names[-start])
    return result

if __name__ == '__main__':
    dnaseq = DNASeq('TGGGCTATTATCGACCACAAAGTGGTACCAGTGTACGGCGACCTGGACGGCGACGGAGAGGTGAACAGCTTCGACAGGCTGAGGATGGCCCAGGCCGTGCTGAACGGCGACACCGAGAGGTTCCAGGCCGCCGACCTGAACTGCGACGGCGTGATCGACGAGGACGACCTGAAGTACCACAGCGAGTACCTGCTGGGCAAGAGGAAGACCCTGCCCGTGGAGTACTAAGGATCCGGCTGCTAACAAAGCCCG')
    print(find_protein_coding_regions(dnaseq))
