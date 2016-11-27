#!/usr/bin/env python3.5
"""
calculate a span's ∆G of insertion. if a sequence is entered, it's ∆G is calculated by sequence only
if a PDB is entered (has '.pdb' in the file name) than the ∆G is calculated using actual Z positions
"""
import numpy as np
import sys
from scipy import interpolate
# from CreateEMPBenchmark import POS_Z_DICT_total, Z_total, AAs, get_elazar_scale
from retrive_natural_TMs_scores import calc_rosetta_splines, calc_seq_dg
from MyPDB import parse_PDB, memb_residues
from MP_utils import find_spans_pdb, calc_span_dg

def main():
    global rosetta_splines
    rosetta_splines = calc_rosetta_splines()
    arg1 = sys.argv[1]
    if '.pdb' not in arg1:
        seq = arg1
        ori = sys.argv[2]
        dg = calc_seq_dg(seq, ori, rosetta_splines)
        print('%s %s %.2f' % (seq, ori, dg))

    else:
        pdb = parse_PDB(arg1)
        spans = find_spans_pdb(pdb)
        for span in spans:
            dg = calc_span_dg(pdb, span[0], span[1], rosetta_splines)
            print('span %i->%i ∆G=%.2f' % (span[0], span[1], dg))

if __name__ == '__main__':
    main()
