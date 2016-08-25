"""
a class for making log files, during stdout print
"""
import os
import colorama
import sys
import time


class Logger:
    def __init__(self, log_file, path=None, err_file=None, emblem=None):
        if path is None:
            path = os.getcwd()+'/'
        elif path[-1] != '/':
            path += '/'
        self.emblem = emblem if emblem is not None else '~'
        self.HOME = os.environ['HOME']
        self.log_file = open(path+log_file, 'w+')
        self.log('given command: %s' % ' '.join(sys.argv))
        self.files_logged = {}

    def log(self, string, to_print=True, emphasize=False, time_stamp=True):
        ts = "@%s" % time.strftime("%H:%M") if time_stamp else ""
        pt = os.path.expanduser(os.getcwd()).replace(self.HOME, '~')
        if to_print:
            if emphasize:
                print(colorama.Fore.RED + colorama.Back.GREEN + string + colorama.Style.RESET_ALL)
            else:
                print('<%s%s> %s' % (pt, ts, string))
        self.log_file.write('<%s%s> %s\n' % (pt, ts, string))
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
        self.log_file.close()
