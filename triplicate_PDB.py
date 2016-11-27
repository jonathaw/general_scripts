#!/usr/bin/env python3.5
"""
a script to add a membrane cage to a pdb
"""
from MyPDB import parse_PDB, XYZ, write_PDB
import sys
import copy

def main():
    pdb = parse_PDB(sys.argv[1])
    for name, xyz in zip(['B', 'C', 'X'], [XYZ(7.5, 0, 0), XYZ(0, 7.5, 0), XYZ(7.5, 7.5, 7.5)]):
        new_chain = copy.deepcopy(pdb.chains['A'])
        new_chain.change_chain_name(name)
        new_chain.translate_xyz(xyz)
        pdb.add_chain(new_chain)
        pdb.renumber()
    write_PDB('test.pdb', pdb)

if __name__ == '__main__':
    main()
