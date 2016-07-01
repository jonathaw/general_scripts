#!/usr/bin/env python3.5
"""
script to parse a pdb file and output the X (membrane Z) coordinate and ResSolv energy of every residue
"""
import sys
import numpy as np
from MyPDB import parse_PDB, parse_energy_table


def main():
    pdb_file = sys.argv[1]
    pdb = parse_PDB(pdb_file)
    energy_table = parse_energy_table(pdb_file)
    polyval = MakeHydrophobicityGrade()

    for res in pdb.iter_all_res():
        cnt_6 = pdb.count_atoms_near_res(res, 6)
        cnt_12 = pdb.count_atoms_near_res(res, 12)
        print(res, burial_sigmoid(cnt_6, cnt_12))


    # for res in pdb.iter_all_res():
    #     dg = round(np.polyval(polyval[res.res_type], res['CA'].xyz.z), 1)
    #     print('%s %i z: %.2f res_solv: %.2f dG: %.2f' % (res.res_type, res.res_num, res['CA'].xyz.z,
    #                                                      energy_table[energy_table['res_type_num'] == '%s_%i' % (res.res_type_3, res.res_num)].res_solv.values[0], dg))


def burial_sigmoid(cnt_6_, cnt_12_):
    sig_6 = 1. / ( 1. + np.exp(cnt_6_ - 30.) * 0.05)
    sig_12 = 1. / (1. + np.exp(cnt_12_ - 350.) * 0.05)
    return sig_6 * sig_12

def MakeHydrophobicityGrade():
    """
    :return: returns a dictionary of the polynom values for each residue
    """
    hydrophobicity_grade = '/home/labs/fleishman/jonathaw/membrane_prediciton/polyval_val_by_ile_symm_3Apr.txt'
    hydrophobicity_grade = open(hydrophobicity_grade, 'r')
    hydrophobicity_polyval = {}
    for line in hydrophobicity_grade:
        split = line.split()
        hydrophobicity_polyval[split[0]] = [float(n) for n in split[1:6]]
    hydrophobicity_grade.close()
    return hydrophobicity_polyval

if __name__ == '__main__':
    main()
