from cStringIO import StringIO
import sys


class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout


def log_maker(log, log_file='log_maker.log'):
    import sys
    cmd_line = sys.argv
    with open(log_file, 'wr+') as l:
        print "%s ".join()