#!/usr/bin/env python3.5
import sys
import time
import subprocess


def main(user):
    iter_num = 1
    while True:
        stat = get_status(user)
        sys.stdout.write('%i@%s RUN: %i PEND: %i\r' % (iter_num, time.strftime("%H:%M:%S"), stat['running'],
                                                       stat['pending']))
        sys.stdout.flush()
        iter_num += 1
        time.sleep(5)
        # break


def get_status(user):
    proc = subprocess.Popen(['bjobs', '-u', user], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = str(proc.stdout.read())
    stderr = str(proc.stderr.read())
    if stderr == 'No unfinished job found':
        return 'No unfinished job found'
    else:
        data = stdout.split('\\n')
        result = {'running': 0, 'pending': 0}
        for l in data[1:]:
            s = l.split()
            if len(s) >= 2:
                if s[2] == 'RUN':
                    result['running'] += 1
                elif s[2] == 'PEND':
                    result['pending'] += 1
        return result


if __name__ == '__main__':
    main('jonatha')
