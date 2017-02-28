#!/usr/bin/env python3.5
"""
a script to add a membrane cage to a pdb
"""
import MyPDB
import sys


def main():
    pdb_name = sys.argv[1]
    output_file = pdb_name.split('.pdb')[0] + '_MBR.pdb'
    pdb = MyPDB.parse_PDB(pdb_name)

    memb_chain = determine_membrane_chain(pdb)
    highet_atom_number = determine_highest_serial_num(pdb)
    create_and_add_MBR_residue(pdb, memb_chain, highet_atom_number)

    MyPDB.write_PDB(output_file, pdb)

    # append_CONECT_to_PDB(output_file, highet_atom_number)


def determine_membrane_chain(pdb):
    memb_chain = 'Z'
    if 'M' not in pdb.chains.keys():
        memb_chain = 'M'
    else:
        for n in list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')[::-1]:
            if n not in pdb.chains.keys():
                memb_chain = n
                break
    return memb_chain


def determine_highest_serial_num(pdb):
    return max([a.serial_num for c in pdb.chains.values() for r in c.values() for a in r.values()])


def create_and_add_MBR_residue(pdb, memb_chain, highest_atom_number):
    inn = MyPDB.Residue('INN', 1, memb_chain)
    atom_num = 1
    i=1
    for x in range(-30, 30):
        for y in range(-30, 30):
            inn.add_atom(MyPDB.Atom(serial_num=highest_atom_number+atom_num, name='%i' % i , res_type_3='INN', chain=memb_chain, res_seq_num=1, x=x, y=y, z=-15, element='N', charge='', occupancy=1, temp=1, header='HETATM', alternate='', achar='', si=''))
            atom_num += 1
            i += 1

    out = MyPDB.Residue('OUT', 2, memb_chain)
    i = 1
    for x in range(-30, 30):
        for y in range(-30, 30):
            out.add_atom(MyPDB.Atom(serial_num=highest_atom_number+atom_num, name='%i' % i, res_type_3='MID', chain=memb_chain, res_seq_num=2, x=x, y=y, z=15, element='O', charge='', occupancy=1, temp=1, header='HETATM', alternate='', achar='', si=''))
            atom_num += 1
            i += 1

    mid = MyPDB.Residue('MID', 3, memb_chain)
    i = 1
    for x in range(-30, 30):
        for y in range(-30, 30):
            mid.add_atom(MyPDB.Atom(serial_num=highest_atom_number+atom_num, name='%i' % i, res_type_3='MID', chain=memb_chain, res_seq_num=3, x=x, y=y, z=0, element='C', charge='', occupancy=1, temp=1, header='HETATM', alternate='', achar='', si=''))
            atom_num += 1
            i += 1

    chn = MyPDB.Chain(memb_chain, {1: inn, 2: out, 3: mid})
    chn.add_residue(inn)
    chn.add_residue(out)
    chn.add_residue(mid)

    pdb.add_chain(chn)


def append_CONECT_to_PDB(pdb_file, highet_atom_number):
    with open(pdb_file, 'a') as fout:
        for i, j in zip([1, 2, 3, 4], [2, 3, 4, 1]):
            fout.write('CONECT %i %i\n' % (highet_atom_number+i, highet_atom_number+j))
        for i, j in zip([5, 6, 7, 8], [6, 7, 8, 5]):
            fout.write('CONECT %i %i\n' % (highet_atom_number + i, highet_atom_number + j))
        fout.write('CONECT %i %i\n' % (highet_atom_number+1, highet_atom_number+5))
        fout.write('CONECT %i %i\n' % (highet_atom_number + 2, highet_atom_number + 6))
        fout.write('CONECT %i %i\n' % (highet_atom_number + 3, highet_atom_number + 7))
        fout.write('CONECT %i %i\n' % (highet_atom_number + 4, highet_atom_number + 8))


if __name__ == '__main__':
    main()
