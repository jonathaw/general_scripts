"""
A script to align and analyse RB out results from Rafael's MD simulation
"""
class RBVector():
    def __init__(self, first, second):
        self.first = {'x':}
def rb_parser(file_name):
    vectors = []
    with open(file_name, 'r') as f:
        for line in f:
            split = line.split()