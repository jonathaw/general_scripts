#!/usr/bin/env python3.5
import argparse
from MyPDB import MyPDB, parse_PDB, write_PDB


def parse_Bfactor_list(file_name: str) -> dict:
    """
    return a dict of res numbers pointing to list value normalized to 1-50
    """
    result = {}
    i = 1
    for l in open(file_name, 'r'):
        s = l.split()
        if len(s) == 2:
            try:
                result[int(s[0])] = float(s[1])
            except:
                pass
        if len(s) == 1:
            if s[0] == '#NUM!':
                s[0] = 0
            result[i] = float(s[0])
            i += 1
    talk = False
    min_ = min(result.values()) + 0.01
    print(min_)
    for k, v in result.items():
        print(k, v, v-min_)
        result[k] = v - min_
    for k, v in result.items():
        if v < 0.001:
            result[k] = 0
            talk = True

    if talk:
        print('FORCED VALUS TO BE 0')
    min_ = min(result.values())
    max_ = max(result.values())
    new_results = {}
    for k, v in result.items():
        print(k, v)
        new_results[k] = (v - min_) / ( max_ - min_ ) * 50.0
    return new_results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-pdb', type=str, help='PDB file')
    parser.add_argument('-list', type=str, help='B factor list. residue number\t number')
    parser.add_argument('-mode', default='list')

    args = vars(parser.parse_args())

    if args['mode'] == 'list':
        Bfactor_by_list(args)

    elif args['mode'] == 'z':
        Bfactor_by_z(args)

    else:
        print('no mode chosen')


def Bfactor_by_z(args):
    pdb = parse_PDB(args['pdb'])
    z_Bfactor_dict = {}
    for rsd in pdb.iter_all_res():
        if -15 <= rsd['CA'].xyz.z <= 15:
            z_Bfactor_dict[rsd.res_num] = 50
        else:
            z_Bfactor_dict[rsd.res_num] = 0
    new_file = args['pdb'].split('.pdb')[0]+'_z.pdb'
    repalce_Bfactor(pdb, z_Bfactor_dict, args, new_file)



def Bfactor_by_list(args):
    pdb = parse_PDB(args['pdb'])

    resnum_occupancy_dict = parse_Bfactor_list(args['list'])
    file_name = args['pdb'].split('.pdb')[0] + '_' + args['list'].split('.')[0] + '.pdb'
    repalce_Bfactor(pdb, resnum_occupancy_dict, args, new_name=file_name)


def repalce_Bfactor(pdb: MyPDB, resnum_occupancy_dict_: dict, args: dict, new_name: str) -> None:
    resnum_occupancy_dict = {}
    min_ = min(resnum_occupancy_dict_.values())
    max_ = max(resnum_occupancy_dict_.values())
    for k, v in resnum_occupancy_dict_.items():
        resnum_occupancy_dict[k] = (v - min_) / (max_ - min_) * 50.0

    for cid in sorted(pdb.chains.keys()):
        for rid in sorted(pdb[cid].residues.keys()):
            for aid in sorted(pdb[cid][rid].keys()):
                if pdb[cid][rid].res_num in resnum_occupancy_dict.keys():
                    pdb[cid][rid][aid].set_temp(resnum_occupancy_dict[pdb[cid][rid].res_num])
                else:
                    pdb[cid][rid][aid].set_temp(0.0)

    write_PDB(new_name, pdb)
    print('writing the pdb into %s' % new_name)


if __name__ == '__main__':
    main()
