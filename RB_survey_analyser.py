#!/usr/bin/env python3.5
from MyPDB import XYZ


class RigidBody():
    def __init__(self, p1: XYZ, p2: XYZ, p3: XYZ, p4: XYZ):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4

    def __str__(self):
        return '%r\n%r\n%r\n%r' % (self.p1, self.p2, self.p3, self.p4)

    def __repr__(self):
        return self.__str__()

def read_RBs(file_name):
    results = []
    with open(file_name, 'r') as fin:
        for l in fin:
            floats = [float(a) for a in l.split()[1:]]
            rb = RigidBody(XYZ(floats[0], floats[1], floats[2]), XYZ(floats[3], floats[4], floats[5]),
                           XYZ(floats[6], floats[7], floats[8]), XYZ(floats[9], floats[10], floats[11]))
            results.append(rb)
            break
    print(results)


if __name__ == '__main__':
    read_RBs('all_dir_0_RB.db')