#!/usr/bin/env python3.5
__author__ = 'jonathan'
from AASeq import compare_2_seqs
import seq_funcs as sf
import MyPDB


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', type=str, default='pdbs', help='different running modes.\npdbs=from two pdbs')
    parser.add_argument('-1', type=str, help='first entry')
    parser.add_argument('-2', type=str, help='second entry')
    args = vars(parser.parse_args())
    if args['mode'] == 'pdbs':
        pdb_1 = MyPDB.parse_PDB(args['1'], args['1'])
        pdb_2 = MyPDB.parse_PDB(args['2'], args['2'])
        seq_1 = MyPDB.extract_seq(pdb_1)
        seq_2 = MyPDB.extract_seq(pdb_2)
        start = 0
        sorted_keys = sorted(seq_1.keys())
        for k in sorted_keys:
            compare_2_seqs(seq_1[k], seq_2[k], start=start)
            start += len(seq_1[k])
            print(start, len(seq_1[k]))
