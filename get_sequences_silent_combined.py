import sys
import re

file_name = sys.argv[1]
seq_re = re.compile(r'(\[.+?\])')
with open(file_name, 'r') as f:
    for line in f:
        if line.split()[0] == 'ANNOTATED_SEQUENCE:':
            print '>' + line.split()[-1]
            print re.sub(seq_re, '', line.split()[1]) + '\n'