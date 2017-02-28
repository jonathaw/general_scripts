#!/usr/bin/env python3.5
"""
"""
import pandas as pd

def main():
    mut_df = parse_mutations_table()


def parse_mutations_table() -> pd.DataFrame:
    df = pd.read_csv( '/home/labs/fleishman/jonathaw/elazaridis/mut_recap_20Feb/bAR_20Feb/data/experimental_data/processed_table.tsv', sep='\t' )
    print( df )



if __name__ == '__main__':
    main()
