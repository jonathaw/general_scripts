#!/usr/bin/env python3.5
from MyPDB import parse_PDB
from sys import argv
name = argv[1]
pdb = parse_PDB(file_in=name, name=name)
seq = pdb.seqs
if len(argv) > 2:
    chains = argv[2:]
    for c in chains:
        print('>%s.%s' % (name, c))
        print(seq[c].get_seq)
else:
    for c in sorted(seq.keys()):
        if 'non_res' not in c:
            print('>%s.%s' % (name, c))
            print(seq[c].get_seq)
