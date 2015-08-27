#!/usr/bin/env python
def main():
    import argparse
    from random import shuffle
    from math import floor
    parser = argparse.ArgumentParser()
    parser.add_argument('-file')
    parser.add_argument('-num_parts', default=2, type=int)
    args = vars(parser.parse_args())
    in_list = read_list(args['file'])
    shuffle(in_list)
    part_size = int(floor(float(len(in_list)) / float(args['num_parts'])))
    print 'cutting to %i parts', part_size
    gotto = 0
    result = []
    print 'in list length', len(in_list)
    for i in range(args['num_parts']):
        res = []
        for a in range(gotto, gotto+part_size):
            try:
                res.append(in_list[a])
            except:
                continue

        gotto += part_size
        print res, len(res)
        result.append(res)
    test = [i for a in result for i in a]
    print set([x for x in test if test.count(x) > 1])
    for i, a in enumerate(result):
        write_list('split_list_'+str(i), a)


def write_list(file_o, l):
    with open(file_o, 'wr+') as fout:
        for a in l:
            fout.write('%s\n' % a)


def read_list(file_i):
    with open(file_i, 'r') as fin:
        cont = fin.read().split('\n')
    return [a for a in cont if a != '']


if __name__ == '__main__':
    main()