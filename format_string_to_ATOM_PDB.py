#!/usr/bin/env python3.5
"""
a script to read a poorly formated atom list and print it as PDB format
"""
import MyPDB
import sys
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', default=None)
    parser.add_argument('-xyz', default=None, type=str)
    parser.add_argument('-num', default=100, type=int)
    parser.add_argument('-name', default='NNN')
    parser.add_argument('-res_name', default='X')
    parser.add_argument('-chain', default='X')
    parser.add_argument('-res_num',default=500 )
    args = vars(parser.parse_args())

    if args['f'] is not None:
        by_file(args)

    elif args['xyz'] is not None:
        by_xyz(args)

    else:
        print('no mode chosen')


def by_xyz(args):
    print(args['xyz'])
    i = args['num']
    j = args['res_num']
    for l in open(args['xyz'], 'r'):
        s = l.rstrip().split()
        a = MyPDB.Atom('HETATM', i, args['name'], "", args['res_name'], args['chain'], j, float(s[0]), float(s[1]), float(s[2]), "", charge=0, element="O", si="", occupancy=1, temp=0)
        print(a)
        i += 1
        j += 1

def by_file(args):
    atoms = []
    for l in open(args['f'], 'r'):
        s = l.rstrip().split()
        header = s[1]
        num = int(s[2])
        name = s[3]
        res_name = s[4]
        chain = s[5]
        res_num = int(s[6])
        x = float(s[7])
        y = float(s[8])
        z = float(s[9])

        a = MyPDB.Atom(header, num, name, "", res_name, chain, res_num, x, y, z, "", charge=0, element="O", si="", occupancy=1, temp=0)
        print(a)


if __name__ == '__main__':
    main()
