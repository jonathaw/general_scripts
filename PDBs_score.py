#!/usr/bin/env python3.5
from RosettaFilter import score2dict
import sys
__author__ = 'jonathan'

score_dict = score2dict(sys.argv[1])
print(score_dict[sys.argv[2]])