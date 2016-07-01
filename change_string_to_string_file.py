#!/usr/bin/env python3.5
"""
a script to change a string in a txt file to another
"""
import re
import sys
import zipfile
import gzip
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-file_src', default=sys.argv[1], type=str, help='source file name')
    parser.add_argument('-str_src', default=sys.argv[2], type=str, help='string to replace')
    parser.add_argument('-str_dst', default=sys.argv[3], type=str, help='string to insert')
    parser.add_argument('-file_dst', default=sys.argv[4], type=str, help='file destination')
    parser.add_argument('-file_type', default='txt', type=str, help='file type')
    args = vars(parser.parse_args())

    if args['file_type'] == 'zip':
        zipped(args)

    elif args['file_type'] == 'gzip':
        gzipped(args)

    elif args['file_type'] == 'txt':
        not_zipped(args)

    elif args['file_type'] == 'gzip_rgx':
        gzipped_rgx(args)

    elif args['file_type'] == 'txt_rgx':
        txt_rgx(args)

    else:
        print('file type not recognised')


def gzipped(args):
    with gzip.open(args['file_dst'], 'wb+') as fout:
        with gzip.GzipFile(args['file_src'], 'r') as zin:
            for l in zin:
                fout.write(l.replace(str.encode(args['str_src']), str.encode(args['str_dst'])))


def txt_rgx(args):
    name = re.compile('S_[0-9]{4}_[0-9]*')
    with open(args['file_dst'], 'wb+') as fout:
        with open(args['file_src'], 'rb') as zin:
            for l in zin:
                name_ = re.findall(name, str(l))
                if name_ != []:
                    orginal_name = name_[0]
                    new_name = '%s_S_%s' % (orginal_name.split('_')[-1], orginal_name.split('_')[1])
                    fout.write(l.replace(str.encode(orginal_name), str.encode(new_name)))


def gzipped_rgx(args):
    name = re.compile('S_[0-9]{4}_[0-9]*')
    with gzip.open(args['file_dst'], 'wb+') as fout:
        with gzip.GzipFile(args['file_src'], 'r') as zin:
            for l in zin:
                name_ = re.findall(name, l)[0]
                print(name_)
                fout.write(l.replace(str.encode(args['str_src']), str.encode(args['str_dst'])))


def zipped(args):
    with zipfile.ZipFile(args['file_src'], 'r') as zin:
        print(zin, type(zin))
        # with zin.open(args['file_src']) as fin:
        #     for l in fin:
        #         print(l)


def not_zipped(args):
    with open(args['file_dst'], 'w+') as fout:
        with open(args['file_src'], 'r') as fin:
            for l in fin:
                print(l)
                fout.write(l.replace(args['str_src'], args['str_dst']))

if __name__ == '__main__':
    main()
