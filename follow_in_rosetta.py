#!/usr/bin/env python3.5
"""
gather all module names used in python scripts
"""
import os
import re
import pathlib
import argparse
import mimetypes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-func_name')
    parser.add_argument('-root', default='./')
    parser.add_argument('-num', default=5, type=int)
    args = vars(parser.parse_args())

    walk_dict = generate_walk_dict(args['root'])
    find_occurences_of(args['func_name'], args['root'], walk_dict)


def find_occurences_of(name, path, walk_dict) -> list:
    for k, v in walk_dict.items():
        for v2 in v:
            new_names = find_func_names_in_file('%s/%s' % (k, v2), name)


def find_func_names_in_file(file_path, name) -> list:
    print('A')
    if mimetypes.guess_type(file_path)[0] == 'text/plain':
        regex = re.compile("}\n*\D*::(\D*)[(]\D*%s[(]" % name)
        print('B', regex)
        try:
            l = regex.findall(open(file_path, 'r').read())
            print('C', l)
            print(open(file_path, 'r').read())
            if l != []:
                print(l)
        except:
            pass
    return []


def generate_walk_dict(root) -> dict:
    walk = os.walk(root)
    walk_dict = {}
    for a in walk:
        walk_dict[a[0]] = a[2]
    return walk_dict


if __name__ == '__main__':
    find_func_names_in_file('/home/labs/fleishman/jonathaw/Rosetta/main/source/src/core/scoring/membrane/MPSpanInsertionEnergy.cc', 'calc_span_score')
    # main()
