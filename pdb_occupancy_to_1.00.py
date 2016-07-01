"""
replace occupancy 0.00 to 1.00 so that rosetta loads them
"""
import sys, os
with open(sys.argv[1], 'r')as f:
    with open('./temp.pdb', 'wr+') as o:
        for line in f:
            edit = list(line)
            # print {k:v for k, v in enumerate(edit)}
            # break
            edit[56] = '1'
            o.write(''.join(edit))
os.rename('./temp.pdb', sys.argv[1])