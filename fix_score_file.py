#!/usr/bin/env python3.5
"""
fix score files with more than one header
"""
from RosettaFilter import score_file2df
import pandas as pd
import sys
import re
from io import StringIO


def main():
    file_name = sys.argv[1]
    df = concat_multi_tables(file_name)
    print(df.to_string())


def concat_multi_tables(file_name):
    header = re.compile('SCORE:.*description.*')
    tables = []
    cont = open(file_name, 'r').read()

    all_ps = []
    for p in re.finditer(header, cont):
        all_ps.append( [ p.start(), p.end() ] )

    paras = []
    major_df = pd.DataFrame()
    for i in range(len(all_ps)-1):
        parag = cont[ all_ps[i][0] : all_ps[i+1][0] ]

        s = StringIO( parag )
        tables.append(pd.read_csv(s, sep='\s+'))

    parag = cont[ all_ps[-1][0]: ]
    s = StringIO( parag  )
    tables.append(pd.read_csv(s, sep='\s+'))
    major_df = pd.concat(tables, ignore_index=True)

    # rearrange columns so description is last
    cols = major_df.columns.tolist()
    cols.remove( 'description' )
    cols.append( 'description' )
    if 'description.1' in cols:
        cols.remove('description.1')
    major_df = major_df[ cols ]


    return major_df


if __name__ == '__main__':
    main()
