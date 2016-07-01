#!/usr/bin/env python
"""
A script to make a .pml to align two pdbs
1st argument is the full name (path+name) to the unmobile pdb
2nd argument is the full path (path+name) to the mobile pdb
creates a .pml file with same name as the mobile with _aln.pdb
use this with pymol.
examples:
applying this script to many pdbs, with same target
for file in dir_*/1ohz_*.pdb;do python ~/scripts/jonathaw/job_makers/make_pml_to_align_pdbs.py
                                                        ~/data/PDB_all/1ohz_AB.pdb ${file};done
using the .pml files made here:
for pml in dir_*/*pml;do sh /home/labs/fleishman/jonathaw/bin/fleish_sub_pymol_new_all_q_gideon.sh ${pml} ;done
"""
import sys
target_in = sys.argv[1]
mobile_in = sys.argv[2]
path2mobile = '/'.join(mobile_in.split('/')[:-1])+'/'
mobile_name = mobile_in.split('/')[-1].split('.')[0]
target_name = target_in.split('/')[-1].split('.')[0]
with open(path2mobile+mobile_name+'_aln.pml', 'wr+') as o:
    o.write('load %s\n' % target_in)
    o.write('load %s\n' % mobile_in)
    o.write('cealign %s, %s\n' % (target_name, mobile_name))
    o.write('save %s, %s' % (path2mobile+mobile_name+'_aln.pdb', mobile_name))