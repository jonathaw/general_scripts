#!/usr/bin/env python3.5
for kind in ['pending', 'processed', 'running']:
    num_lines = sum(1 for line in open('/home/labs/fleishman/jonathaw/general_lists/%s_folders.txt' % kind, 'r'))
    print('found %i lines in %s' % (num_lines, kind))
num_lines = sum(1 for line in open('/home/labs/fleishman/jonathaw/general_lists/switched_jobs.txt', 'r'))
print('found %i lines in switched_jobs' % num_lines)