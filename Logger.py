#!/usr/bin/env python3.5
"""
a class for making log files, during stdout print
"""
import os
import colorama
import inspect
import sys
import time


class Logger:
    def __init__(self, log_file, path=None, err_file=None, emblem=None):
        self.started = time.time()
        self.HEIGHT, self.WIDTH = (int(a) for a in os.popen('stty size', 'r').read().split())
        if path is None:
            path = os.getcwd()+'/'
        elif path[-1] != '/':
            path += '/'
        self.emblem = emblem if emblem is not None else '~'
        self.HOME = os.environ['HOME']
        self.log_file = open(path+log_file, 'w+')
        header = ' Starting %s ' % sys.argv[0]
        self.create_header(header)
        self.log('given command: %s' % ' '.join(sys.argv))
        self.files_logged = {}

    def log(self, string, to_print=True, emphasize=False, time_stamp=True, skip_stamp: bool=False):
        # create time stamp
        log_ts = "%s" % time.strftime("%H:%M") if time_stamp else ""
        ts = colorama.Fore.RED + log_ts + colorama.Style.RESET_ALL

        # create path stamp
        log_pt = os.path.expanduser(os.getcwd()).replace(self.HOME, '~')
        pt = colorama.Fore.GREEN + log_pt + colorama.Style.RESET_ALL

        # create calling frame stamp
        ins = inspect.getouterframes( inspect.currentframe()  )
        (frame, filename, lineno, function, code_context, index) = ins[1]
        log_lc ='%s$%s:%i' % (filename.split('/')[-1], function, lineno)
        lc = colorama.Fore.MAGENTA + log_lc + colorama.Style.RESET_ALL

        # shrink path stamp
        if pt.count('/') > 1:
            s = pt.split('/')
            log_pt = "%s" % s[-1]
            pt = colorama.Fore.GREEN + log_pt + colorama.Style.RESET_ALL

        if len('<%s@%s@%s> %s' % (pt, ts,lc, string)) > self.WIDTH-2:
            string = '\n--> %s' % string

        # print or not and log
        if to_print:
            if emphasize:
                print(colorama.Fore.RED + colorama.Back.GREEN + string + colorama.Style.RESET_ALL)
            elif not skip_stamp:
                print('<%s@%s@%s> %s' % (pt, ts,lc, string))
                self.log_file.write('<%s@%s@%s> %s\n' % (log_pt, log_ts, log_lc, string))
            else:
                print(string)
                self.log_file.write('%s\n' % string)
        self.log_file.flush()
        sys.stdout.flush()

    def log_text_file(self, file_name: str, to_print: bool = True):
        if file_name in self.files_logged.keys():
            self.log('file logged at %i' % self.files_logged[file_name])
        else:
            self.files_logged[file_name] = len(list(self.files_logged.keys()))+1
            self.log('logging the file %s at num %i' % (file_name, len(list(self.files_logged.keys()))+1))
            for l in open(file_name, 'r'):
                if to_print:
                    print(l.rstrip())
                self.log_file.write(l)
            self.log_file.write('\n\n')

    def close(self):
        time_delta = time.time() - self.started
        h, remain = divmod(time_delta, 3600)
        m, s = divmod(remain, 60)
        self.log('finished at %s, process took %i:%i:%i' % (time.strftime("%H:%M"), h, m, s))
        self.create_header('Finishing %s' % sys.argv[0])
        self.log_file.close()

    def create_header(self, msg: str='') -> None:
        print('='*self.WIDTH)
        print('%s%s%s' % ('='*int((self.WIDTH-len(msg))/2), colorama.Fore.BLUE + msg + colorama.Style.RESET_ALL, '='*int((self.WIDTH-len(msg))/2)))
        print('='*self.WIDTH)
        self.log_file.write('%s%s%s\n' % ('='*int((self.WIDTH-len(msg))/2), msg, '='*int((self.WIDTH-  len(msg))/2)))
        self.log_file.flush()
        sys.stdout.flush()
