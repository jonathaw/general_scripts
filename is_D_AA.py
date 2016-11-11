#!/usr/bin/env python3.5
"""
a script to check if an amino acid in a pdb is D or L amino acid
"""
import MyPDB
import sys
import numpy as np


def main():
    pdb_name = sys.argv[1]
    residue_num = int(sys.argv[2])
    pdb = MyPDB.parse_PDB(pdb_name)

    res = pdb['A'][residue_num]
    CO = np.array([res['C'].xyz.x, res['C'].xyz.y, res['C'].xyz.z])
    CA = np.array([res['CA'].xyz.x, res['CA'].xyz.y, res['CA'].xyz.z])
    CB = np.array([res['CB'].xyz.x, res['CB'].xyz.y, res['CB'].xyz.z])
    N = np.array([res['N'].xyz.x, res['N'].xyz.y, res['N'].xyz.z])
    # HA = np.array([res['HA'].xyz.x, res['HA'].xyz.y, res['HA'].xyz.z])

    v1 = N - CO
    v2 = CA - CO

    cp = np.cross(v1, v2)
    # HA_infront = cp.dot(HA-CA) > 0
    CB_infront = cp.dot(CB-CA) > 0
    print('residue %r, is %s ' % (res, res.D_or_L()))

if __name__ == '__main__':
    main()
