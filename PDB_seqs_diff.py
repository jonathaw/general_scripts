#!/usr/bin/env python3.5
from MyPDB import parse_PDB, extract_seq
import sys

in_1 = sys.argv[1]
in_2 = sys.argv[2]

pdb_1 = parse_PDB(in_1, in_1.split('.')[0])
pdb_2 = parse_PDB(in_2, in_2.split('.')[0])

seqs_1 = extract_seq(pdb_1)
seqs_2 = extract_seq(pdb_2)
assert seqs_1.keys() == seqs_2.keys(), 'PDBs have different chain names'

seq_1, seq_2 = '', ''
for chain_name in sorted(seqs_1.keys()):
	seq_1 += seqs_1[chain_name].get_seq
	seq_2 += seqs_2[chain_name].get_seq
assert len(seq_1) == len(seq_2), 'sequences not the same length'
diffs = []
for i in range(len(seq_1)):
	if seq_1[i] != seq_2[i]:
		print('%s%i%s' % (seq_1[i], i+1, seq_2[i]))
		diffs.append(i+1)
print('found %i changes, over %i chains' % (len(diffs), len(seqs_1.keys())))
print('select diffs, resi %s' % '+'.join([str(a) for a in diffs]))