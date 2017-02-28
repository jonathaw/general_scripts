#!/usr/bin/env python3.5
"""
a script to gaet spans from a topo string
"""
import sys
from prepare_fnd import get_spans_from_topo_string

def main():
    ts = sys.argv[1]
    spans = get_spans_from_topo_string( {'topo_string': ts, 'symm_num': 1, '1st_span_orient': None} )
    for span in spans:
        print('%i %i %s' % (span['start'], span['end'], span['orientation']))
    for span in spans:
        print('<Span start="%i" end="%i" orientation="%s"/>' % ( span['start'], span['end'], span['orientation']))
    for i, span in enumerate( spans ):
        print('select span_%i, resi %i-%i' % (i, span['start'], span['end']))

if __name__ == '__main__':
    main()
