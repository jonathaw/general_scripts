#!/usr/bin/env python3.5
import os, re
import time
import argparse
import pickle
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import RosettaFilter as Rf
from Logger import Logger
import pandas as pd
import numpy as np
from collections import OrderedDict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', default='best_in_cluster')
    parser.add_argument('-clust_log', type=str)
    parser.add_argument('-top_n', default=1)
    parser.add_argument('-skip_empty', type=bool, default=True)

    args = vars(parser.parse_args())

    if args['mode'] == 'best_in_cluster':
        best_in_cluster(args)

    else:
        print('no mode chosen!!!')


def best_in_cluster(args):
    logger = Logger('logeer_%s.log' % time.strftime("%d.%0-m"))
    sorted_data = parse_cluster_log(args['clust_log'])
    names = []
    for cluster, data in sorted_data.items():
        print('cluster #%i' % cluster)
        i = 0
        for k, v in data.items():
            if args['skip_empty'] and 'empty' in k:
                continue
            logger.log('\t%s\t%f' % (k, v))
            names.append(k)
            i += 1
            if i == args['top_n']:
                break

    for n in names:
        logger.log(n)


def parse_cluster_log(file_name, verbose=False):
    saw_summary = 0
    results = {}
    current_cluster = 999
    summary_n = 0
    for l in open(file_name, 'r'):
        if '---------- Summary ---------------------------------' in l:
            summary_n += 1
    for l in open(file_name, 'r'):
        if '---------- Summary ---------------------------------' in l:
            saw_summary += 1
        if saw_summary == summary_n:
            s = re.split('\s+', l.rstrip())
            if len(s) == 3 and 'Clusters:' not in s and 'Structures:' not in s:
                results[ current_cluster ][s[1]] = float(s[2])
            elif len(s) == 7 and 'Cluster:' in s:
                current_cluster = int( s[2] )
                results[ current_cluster ] = {}

            if 'Structures' in l:
                break

    sorted_results = {}
    for cluster, data in results.items():
        sorted_results[cluster] = OrderedDict(sorted(data.items(), key=lambda x: x[1]))

    if verbose:
        for cluster, data in sorted_results.items():
            print(cluster)
            for k, v in data.items():
                print('\t%s\t%f' % (k, v))

        print(sorted_results.keys())

    return sorted_results

if __name__ == '__main__':
    main()
