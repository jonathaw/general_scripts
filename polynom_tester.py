#!/usr/bin/env python3.5
import sys
c4 = float(sys.argv[1])
c3 = float(sys.argv[2])
c2 = float(sys.argv[3])
c1 = float(sys.argv[4])
c0 = float(sys.argv[5])
z = float(sys.argv[6])
print('found these c4 %f, c3 %f c2 %f, c1 %f and c0 %f. z is %f' % (c4, c3, c2, c1, c0, z))
res = c4*z**4 + c3*z**3 + c2*z**2 + c1*z +c0
print('result is %f' % res)