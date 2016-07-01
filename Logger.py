"""
a class for making log files, during stdout print
"""
import os
import colorama
import sys


class Logger:
    def __init__(self, log_file, path=None, err_file=None, emblem=None):
        if path is None:
            path = os.getcwd()+'/'
        elif path[-1] != '/':
            path += '/'
        self.emblem = emblem if emblem is not None else '~'
        self.log_file = open(path+log_file, 'w+')
        self.log('given command: %s' % ' '.join(sys.argv))

    def log(self, string, to_print=True, emphasize=False):
        if to_print:
            if emphasize:
                print(colorama.Fore.RED + colorama.Back.GREEN + string + colorama.Style.RESET_ALL)
            else:
                print('<%s> %s' % (self.emblem, string))
        self.log_file.write(string+'\n')
        sys.stdout.flush()

    def close(self):
        self.log_file.close()
