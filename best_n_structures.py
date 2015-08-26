#!/usr/bin/env python
import argparse
import os
import re
from rosetta_score_files import score2dict
parser = argparse.ArgumentParser()
parser.add_argument('-score_file', default=[x for x in os.listdir('./') if re.match('^(?!purple).*\.score', x)][0])
parser.add_argument('-n', default=5, type=int)
args = vars(parser.parse_args())
scores = score2dict(args['score_file'])
passed_scores = {v['ddg']: k for k, v in scores.items() if v['purple']}
sorted_scores = sorted(passed_scores.keys())
fout = open('extract_purples', 'wr+')
try:
    print '\n'.join(passed_scores[sorted_scores[i]] for i in range(args['n']))
    fout.write('\n'.join(passed_scores[sorted_scores[i]] for i in range(args['n'])))
except:
    print 'found nothin'
fout.close()