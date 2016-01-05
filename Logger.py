"""
a class for making log files, during stdout print
"""
import os


class Logger():
    def __init__(self, log_file, path=None, err_file=None):
        if path is None:
            path = os.getcwd()+'/'
        elif path[-1] != '/':
            path += '/'
        self.log_file = open(path+log_file, 'w+')

    def log(self, string, to_print=True):
        if to_print:
            print(string)
        self.log_file.write(string+'\n')

    def close(self):
        self.log_file.close()