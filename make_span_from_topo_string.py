#!/usr/bin/env python3.5

import argparse
from Logger import Logger
from retrive_natural_TMs_scores import find_topo

parser = argparse.ArgumentParser()
parser.add_argument('-ts')
parser.add_argument('-create_text', default=False)
args = vars(parser.parse_args())
args['logger'] = Logger('log_file.log')

spans = find_topo(args['ts'])
print(spans)

if args['create_text']:
    for h in spans:
        print('<Span start=%i end=%i orientation=%s/>' % (h[0], h[1], 'in2out' if h[2] == 'fwd' else 'out2in'))
