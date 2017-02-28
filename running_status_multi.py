#!/usr/bin/env python3.5
import os
import sys
import time
import curses
import subprocess
import functools

HEIGHT, WIDTH = (int(a) for a in os.popen('stty size', 'r').read().split())

@functools.total_ordering
class User:
    def __init__(self, name):
        self.name = name
        self.run = {}
        self.pend = {}
        self.susp = {}
        self.ssusp = {}
        self.unkwn = {}

    def __repr__( self ):
        return "%s: run: %i pend: %i" % ( self.name, self.run, self.pend )

    def __eq__( self, other ):
        return self.total() == other.total()

    def __lt__( self, other ):
        return self.total() > other.total()

    def update_queues( self, q: str ) -> None:
        if q not in self.run.keys(): self.run[q] = 0
        if q not in self.pend.keys(): self.pend[q] = 0
        if q not in self.susp.keys(): self.susp[q] = 0
        if q not in self.ssusp.keys(): self.ssusp[q] = 0
        if q not in self.unkwn.keys(): self.unkwn[q] = 0

    def draw_at(self, pos: int):
        pass

    def graphic_status( self ) -> str:
        fle_run = round( self.run['fleishman'] / 100 )
        fle_pnd = round( self.pend['fleishman'] / 100 )

        new_run = round( self.run['new-all.q']  / 100 )
        new_pnd = round( self.pend['new-all.q'] / 100 )

        msg = '|' * fle_run
        msg += '/' * fle_pnd
        msg += '!' * new_run
        msg += 'i' * new_pnd
        msg += '%i, %i' % (self.run['fleishman'], self.pend['fleishman'])
        msg += '/%i, %i' % (self.run['new-all.q'], self.pend['new-all.q'])
        return msg

    def total( self ) -> int:
        return self.run['fleishman'] + self.run['new-all.q'] + self.pend['fleishman'] + self.pend['new-all.q']


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
            msg = '%s RUN: %i PEND: %i [%s%s]               \r' % (time.strftime("%H:%M:%S"),
                                                        stat[user].run['fleishman'],
                                                        stat[user].pend['fleishman'],
                                                        '-'*i, ' '*(time_to_sleep-1-i))
            sys.stdout.write('\33[%i;%iH%s' % (loc_y,
                                               (WIDTH-len(msg))/2,
                                               msg))
            y_delta = round(HEIGHT/2 - ( len(stat.keys()) -1 )/2 -  1)
            sorted_users = sorted( list(stat.values()) )
            for u in sorted_users:
                if loc_y-2 <= y_delta <= loc_y+2:
                    y_delta += 1
                    continue
                sys.stdout.write('\33[%i;%iH%s\t->\t%s' % ( y_delta, 0, u.name, u.graphic_status() ))
                y_delta += 1

            if stat[user].run['fleishman'] + stat[user].pend['fleishman'] > total_last:
                time_to_sleep = 15
            else:
                time_to_sleep = 6
            total_last = stat[user].run['fleishman'] + stat[user].pend['fleishman']

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
                                'SSUSP': 0, 'USUSP': 0}
                user_list[s[1]] = User(s[1])
            result[s[1]][s[2]] += 1
            queue = s[3]
            user_list[s[1]].update_queues( queue )
            if s[2] == 'RUN':
                user_list[s[1]].run[queue] += 1
            elif s[2] == 'PEND':
                user_list[s[1]].pend[queue] += 1
            elif s[2] == 'SSUSP':
                user_list[s[1]].ssusp[queue] += 1
            elif s[2] == 'unkwn':
                user_list[s[1]].unkwn[queue] += 1
            elif s[2] == 'SUSP':
                user_list[s[1]].sus[queue] += 1

    for k, v in user_list.items():
        v.update_queues( 'fleishman' )
        v.update_queues( 'new-all.q')
    return user_list


if __name__ == '__main__':
    curses.wrapper(main)
