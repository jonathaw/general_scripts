#!/usr/bin/env python3.5
"""
gather all module names used in python scripts
"""
import os
import sys
import argparse

import MyPDB as mp


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-params', default=sys.argv[1])
    args = vars(parser.parse_args())

    params_atoms = parse_params_for_atoms(file_name)


def parse_params_for_atoms(file_name: str) -> dict:
    atoms = {}
    for l in open( file_name, 'r' ):
        if 'ATOM' in file_name:
            s = l.rstrip().split()
            atoms[ s[1] ] = {'RAtomNAme': s[2], 'MMAtomType': s[3], 'charge': float(s[4]), 'charge2': float(s[5])}


if __name__ == '__main__':
    main()
