#!/usr/bin/env python3.5
"""
gather all module names used in python scripts
"""
import os
import pathlib


def part(t):
    if '.' in t:
        return t.split('.')[0]
    else:
        return t


result = set()
for path, subdirs, files in os.walk('/home/labs/fleishman/jonathaw/scripts/'):
    for name in files:
        if name[-3:] == '.py':
            try:
                for l in open(str(pathlib.PurePath(path, name)), 'r'):
                    if 'import' in l:
                        s = l.split()
                        if s[0] == 'import':
                            result.add(part(s[1]))
                        if s[2] == 'import':
                            result.add(part(s[1]))
            except:
                pass
print(result)
