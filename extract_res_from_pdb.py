#!/usr/bin/env python
"""
A script to extract specific residues from a pdb file
"""
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-pdb', type=str)
parser.add_argument('-residues', nargs='+', type=int, default=[151])
parser.add_argument('-atoms', nargs='+', type=str, default=['N', 'CA', 'C', 'O'])
args = vars(parser.parse_args())
pdb_in = args['pdb']
residues = args['residues']
atoms = args['atoms']
pdb_path = '/'.join(pdb_in.split('/')[:-1])+'/'
pdb_name = pdb_in.split('/')[-1].split('.')[0]
in_file = open(pdb_in, 'r')
out_file = open(pdb_path+pdb_name+'_res.pdb', 'wr+')
for line in in_file:
    split = line.split()
    if len(split) < 5: continue
    if int(split[5]) in residues and split[2] in atoms:
        out_file.write(line)
in_file.close()
out_file.close()