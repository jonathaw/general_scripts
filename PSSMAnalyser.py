#!/usr/bin/env python3.5
"""
a script to present PSSMs in a human readable way
"""


def parse_pssm(file_name):
    f = open(file_name, 'r')
    cont = f.read().split('\n')
    result = {}
    aas = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    for l in cont:
        s = l.split()
        if s == []:
            continue
        try:
            num = int(s[0])-1
        except:
            continue
        result[num+1] = {}
        result[num+1]['type'] = s[1]
        for i, sc in enumerate(s[2:22]):
            result[num+1][aas[i]] = int(sc)
    f.close()
    return result


def main(args):
    pssm = parse_pssm(args['pssm_file'])
    for k, v in pssm.items():
        print('pos %i %s: %s' % (k, v['type'], ', '.join([a+'='+str(b) for a, b in v.items() if a != 'type'
                                                          and b >= args['grade']])))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-pssm_file')
    parser.add_argument('-grade', default=0, type=int)

    args = vars(parser.parse_args())
    main(args)