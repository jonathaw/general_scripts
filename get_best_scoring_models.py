#!/usr/bin/env python3
"""
a script to read a score file, and print a list of X% besto scoring models
"""
import argparse
import RosettaFilter as Rf
from Logger import Logger
import pandas as pd
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-sc', type=str, help='score file')
    parser.add_argument('-percent', type=float, default=5, help='percent (1-100) best scoring to get')
    parser.add_argument('-filter', type=str, default='score', help='filter or score term to use')
    parser.add_argument('-num', default=10, type=int, help='use if you want a number of results, not better than percentile')
    parser.add_argument('-mode', default='%')
    parser.add_argument('-over_under', type=str, default='under', help='under/over score should be over/under threshold')
    parser.add_argument('-result', type=str, default=None, help='should the names be written to a file separate from the log file')
    parser.add_argument('-terms', nargs='+', default=['score', 'a_shape', 'a_pack', 'a_ddg', 'res_solv'])
    parser.add_argument('-thresholds', nargs='+', type=float)
    parser.add_argument('-percentile', default=10, type=int)
    args = vars(parser.parse_args())

    logger = Logger('top_%.1f_%s.log' % (args['percent'], args['filter']))

    # read in the score file, determine the threshold for the percentile
    sc_df = Rf.score_file2df(args['sc'])
    score = sc_df[args['filter']]


    if args['mode'] == '%':
        threshold = np.percentile(score, args['percent'])
        logger.log('found a threshold for %f for filter %s to be %.2f' % (args['percent'], args['filter'], threshold))

        # create a df for lines that pass the threshold, either over or above it...
        if args['over_under'] == 'over':
            pass_df = sc_df[sc_df[args['filter']] >= threshold]
        elif args['over_under'] == 'under':
            pass_df = sc_df[sc_df[args['filter']] <= threshold]

    if args['mode'] == 'num':
        sc_df.sort_values(args['filter'], inplace=True)
        pass_df = sc_df.head(args['num'])

    if args['mode'] == 'best_of_best':
        threshold = np.percentile(score, args['percent'])
        sc_df = sc_df[sc_df[args['filter']] <= threshold]
        pass_df = Rf.get_best_of_best(sc_df, args['terms'], args['percentile'])

    if args['mode'] == 'thresholds':
        for term, thrs in zip(args['terms'], args['thresholds']):
            if term in ['a_sasa', 'a_pack', 'a_shape', 'a_tms_span_fa',
                        'a_tms_span', 'a_span_topo']:
                sc_df = sc_df[sc_df[term] > thrs]
            elif term in ['a_mars', 'a_ddg', 'score', 'total_score',
                          'a_res_solv', 'a_span_ins']:
                sc_df = sc_df[sc_df[term] < thrs]
            threshold = np.percentile(score, args['percent'])
            pass_df = sc_df[sc_df[args['filter']] < threshold]

    # output the names (description) of models that pass the threshold, either to the logger file, or to a separate file
    if args['result'] is None:
        logger.create_header('models passing the threshold:')
        for idx, row in pass_df.iterrows():
            logger.log('%s %f' % (row['description'], row['score']), skip_stamp=True)
    else:
        with open(args['result'], 'w+') as fout:
            for name in pass_df['description']:
                fout.write(name + '\n')
    logger.close()

if __name__ == '__main__':
    main()
