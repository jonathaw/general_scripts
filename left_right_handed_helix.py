#!/usr/bin/env python3.5
"""
a script to calc if a helix in a PDB is left or right handed
"""
import MyPDB
import sys


def main():
    pdb_name = sys.argv[1]
    chain = sys.argv[2]
    start, end = sys.argv[3], sys.argv[4]

    pdb = MyPDB.parse_PDB(pdb_name)
    
    for rid, r in pdb[chain]:
        if rid > 1:
            print(rid, r.phi(pdb[chain][rid-1]))
        if rid < len(pdb[chain]):
            print(rid, r.psi(pdb[chain][rid+1]))

    MyPDB.draw_ramachadran(pdb)



if __name__ == '__main__':
    main()
