#!/usr/bin/env python3
"""

"""
import sys
import colorama
import subprocess
import pandas as pd

QUEUES = ['fleishman', 'new-all.q', 'new-short', 'new-medium', 'new-long']


def main():
    bj_df = get_bjobs_df()
    summerise_bjobs_df(bj_df)


def summerise_bjobs_df(df: pd.DataFrame) -> pd.DataFrame:
    user = sys.argv[1]
    g = df.groupby(['QUEUE', 'USER', 'STAT']).size()
    df = pd.DataFrame(g)[0]
    for l in df.to_string().split('\n'):
        s = l.split()
        if user in s:
            print(colorama.Fore.RED + l + colorama.Style.RESET_ALL)
        else:
            print(l)


def get_bjobs_df():
    df = pd.DataFrame()
    task = subprocess.Popen(
        "bjobs -u all |grep -E '(fleishman|new-|USER)'",
        stdout=subprocess.PIPE, shell=1)
    df = pd.read_table(task.stdout, sep='\s+', parse_dates=True,
                       encoding="utf-8", error_bad_lines=False,
                       warn_bad_lines=False, skiprows=1,
                       names=['JOBID', 'USER', 'STAT', 'QUEUE', 'FROM_HOST',
                              'EXEC_HOST', 'JOB_NAME', 'MONTH', 'DAY', 'HOUR'])
    return df

if __name__ == '__main__':
    main()
