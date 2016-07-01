import timeit
import math
start = timeit.default_timer()
import re
import os
import sys
import __main__
__main__.pymol_argv = ['pymol', '-qc']  # Pymol: quiet and no GUI '-qc'
pymol_argv = ['pymol', '-qc']
import pymol
from pymol import cmd
pymol.finish_launching()
        # cmd.load(PDBS_PATH + pdb, object='query')
        # cmd.load(crys_name, object='crys')
        # cmd.select('query_sele', 'query and resi '+aln_thr_seg_str)
        # cmd.select('crys_sele', 'crys and resi '+crys_seg_str)
        # cmd.align('crys_sele', 'query_sele', object='alignment')[0]
        # rmsd = cmd.rms_cur('crys_sele and name ca', 'query_sele and name ca')
        # main_dict[pdb] = rmsd

pdb_list = [x for x in os.listdir(os.getcwd()) if re.match('.*\.pdb', x)]
cwd = os.getcwd() + '/'
name1 = sys.argv[-2] if sys.argv[-1] == '-qc' else sys.argv[-1]
cmd.load(name1, 'q1')
with open(name1[:-4] + '.rmsd', 'wr+') as f:
    for name2 in pdb_list:
        cmd.load(name2, 'q2')
        rmsd = cmd.align('q1 and name CA', 'q2 and name CA', object='alignment')[0]
        cmd.delete('q2')
        f.write(name1 + ' ' + name2 + ' ' + str(rmsd) + '\n')
stop = timeit.default_timer()
print math.floor((stop - start) / 60), " minutes and ", math.floor((stop - start) % 60), "seconds"