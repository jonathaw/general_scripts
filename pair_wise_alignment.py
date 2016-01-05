#!/usr/bin/env python3.5
__author__ = 'jonathan'
from AASeq import read_seq
import sys
seq1 = read_seq(sys.argv[1])
seq2 = read_seq(sys.argv[2])
print('score %f, start %i, end %i' % seq1.align(seq2))
print('seq1 aligned:', seq1.get_aligned())
print('seq2 aligned:', seq2.get_aligned())
print('identity %f' % seq1.aligned_identity(seq2))