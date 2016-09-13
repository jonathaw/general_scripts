#!/usr/bin/env python3.5
import os
import sys
import time
import inspect
import curses
import subprocess

HEIGHT, WIDTH = (int(a) for a in os.popen('stty size', 'r').read().split())


class User:
    def __init__(self, name):
        self.name = name
        self.run = 0
        self.pend = 0
        self.susp = 0
        self.ssusp = 0
        self.unkwn = 0

    def draw_at(self, pos: int):
        pass


def main(stdscr):
    stdscr.clear()
    user = sys.argv[1]
    stdscr.nodelay(True)
    iter_num = 1
    loc_y = HEIGHT / 2
    key = ''
    total_last = 0
    time_to_sleep = 6
    while key != ord('q'):
        # stat = get_status(user)
        stat = get_all_status()
        for i in range(1, time_to_sleep):
            try:
                key = stdscr.getch()
            except:
                pass
            if key == 'q':
                break
            if user not in stat.keys():
                # stat[user] = {'RUN': 0, 'PEND': 0}
                stat[user] = User(user)
            msg = '%i@%s RUN: %i PEND: %i [%s%s]               \r' % (iter_num,
                                                        time.strftime("%H:%M:%S"),
                                                        stat[user].run,
                                                        stat[user].pend, 
                                                        '-'*i, ' '*(time_to_sleep-1-i))
            sys.stdout.write('\33[%i;%iH%s' % (loc_y,
                                               (WIDTH-len(msg))/2,
                                               msg))
            if stat[user].run + stat[user].pend > total_last:
                time_to_sleep = 15
            else:
                time_to_sleep = 6
            total_last = stat[user].run + stat[user].pend
            
            sys.stdout.flush()
            time.sleep(1)
            iter_num += 1
            stdscr.refresh()


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


def get_all_status():
    proc = subprocess.Popen(['bjobs', '-u', 'all'], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    stdout = str(proc.stdout.read())
    stderr = str(proc.stderr.read())
    data = stdout.split('\\n')
    result = {}
    user_list = {}
    for l in data[1:]:
        s = l.split()
        if len(s) >= 2:
            if s[1] not in result.keys():
                result[s[1]] = {'RUN': 0, 'PEND': 0, 'SUSP': 0, 'UNKWN': 0,
                                'SSUSP': 0}
                user_list[s[1]] = User(s[1])
            result[s[1]][s[2]] += 1
            if s[2] == 'RUN':
                user_list[s[1]].run += 1
            elif s[2] == 'PEND':
                user_list[s[1]].pend += 1
            elif s[2] == 'SSUSP':
                user_list[s[1]].ssusp += 1
            elif s[2] == 'unkwn':
                user_list[s[1]].unkwn += 1
            elif s[2] == 'SUSP':
                user_list[s[1]].susp += 1

    return user_list


def other():
    print('in other')
    ins = inspect.getouterframes( inspect.currentframe() )
    (frame, filename, lineno, function, code_context, index) = ins[1]
    print('%s:%i' % (filename.split('/')[-1], lineno))

def yaba():
    print('calling other')
    other()
    print('back')

if __name__ == '__main__':
    yaba()
    sys.exit()
    curses.wrapper(main)
