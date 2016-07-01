#!/usr/bin/env python3.5
import argparse
import sys
from MyPDB import parse_PDB, MyPDB, memb_residues


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', default='select_memb')
    parser.add_argument('-in_file', default=sys.argv[1])
    parser.add_argument('-name', default=None)
    args = vars(parser.parse_args())

    if args['mode'] == 'select_memb_old':
        # load pdb
        pdb = parse_PDB(args['in_file'], args['name'])
        res_in_memb = pymol_select_memb_old(pdb)
        print('select membrane_residues, resi %s' % '+'.join([str(a) for a in res_in_memb]))

    elif args['mode'] == 'print_num_memb_old':
        pdb = parse_PDB(args['in_file'], args['name'])
        res_in_memb = pymol_select_memb_old(pdb)
        print('num_in_memb %i' % len(res_in_memb))

    elif args['mode'] == 'select_memb':
        pdb = parse_PDB(args['in_file'], args['name'])
        memb_ress = memb_residues(pdb)
        print('select membrane_residues, resi %s' % '+'.join([str(r.res_num) for r in memb_ress]))

    else:
        print('no mode selected')




def pymol_select_memb_old(pdb: MyPDB) -> set():
    """
    print a pymol selection line for all residues that are in the membrane
    !!! this assumes that cntr is at 0 0 0 and norm at 15 0 0 !!!
    """
    from shapely.geometry import LineString, Point
    # create Points from center & thickness
    cntr_pnt = Point(pdb.memb_res.cntr.x, pdb.memb_res.cntr.y, pdb.memb_res.cntr.z)
    thkn_m_pnt = Point(-pdb.memb_res.thkn.x, pdb.memb_res.thkn.y, pdb.memb_res.thkn.z)
    thkn_pnt = Point(pdb.memb_res.thkn.x, pdb.memb_res.thkn.y, pdb.memb_res.thkn.z)

    # define the line between center and thickness
    line = LineString([thkn_m_pnt, thkn_pnt])
    thickness = cntr_pnt.distance(thkn_pnt)

    result = set()
    # iterate over all CAs in the pdb
    for cid in sorted(pdb.chains.keys()):
        for rid in sorted(pdb[cid].residues.keys()):
            atom = pdb[cid][rid]['CA']

            # the atom as a Point
            p = Point(atom.xyz.x, atom.xyz.y, atom.xyz.z)

            # projection of the CA atom on the center-thickness line
            np = line.interpolate(line.project(p))

            # if the distance on the center-thickness line is smaller than 15, than this is in the membrane
            if cntr_pnt.distance(np) < thickness-0.1:
                result.add(atom.res_seq_num)

    return result

if __name__ == '__main__':
    main()
