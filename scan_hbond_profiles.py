#!/usr/bin/env python3.5

import argparse
import os
import matplotlib.pyplot as plt
import MyPDB as my
from CreateEMPBenchmark import POS_Z_DICT_total
import RosettaFilter as Rf

def draw_hbonds_profiles():
    parser = argparse.ArgumentParser()
    parser.add_argument('-pdb')
    parser.add_argument('-stage', type=int)
    args = vars(parser.parse_args())

    pdb = my.parse_PDB(args['pdb'])

    if args['stage'] == 1:
        seq_length = pdb.seq_length()

        command = "for i in `seq 1 %i`;do ~/bin/fleish_sub_general.sh /home/labs/fleishman/jonathaw/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease -parser:protocol ~/elazaridis/protocols/scan_hbonds.xml -s %s -mp:scoring:hbond -corrections::beta_nov15 -score:elec_memb_sig_die -score:memb_fa_sol -overwrite -out:prefix ${i}_ -script_vars energy_function=beta_nov15_elazaridis res_num=${i} s1=%i e1=%i ori1=%s s2=%i e2=%i ori2=%s ;done" % (seq_length, args['pdb'], 1, 24, 'out2in', 25, 48, 'out2in')
        print('issuing command\n%s' % command)
        os.system(command)

    if args['stage'] == 2:
        os.system("head -2 1_score.sc|tail -1 > all_score.sc")
        os.system("grep SCORE: *_score.sc|grep -v des >> all_score.sc")
        z_dict = {id: res.memb_z for id, res in pdb.res_items()}
        pos_dict = {v: k for k, v in z_dict.items()}
        sc_df = Rf.score_file2df('all_score.sc')
        zs, scs = [], []
        for d, sc in zip(sc_df['description'].values, sc_df['a_e_res']):
            zs.append(z_dict[ int( d.split('_')[0] ) ])
            scs.append(sc)
        plt.scatter(zs, scs)
        for z, sc in zip(zs, scs):
            if z is not None:
                plt.annotate(pos_dict[z], xy=(z, sc), xytext=(-20, 20),
                            textcoords = 'offset points', ha = 'right', va = 'bottom',
                             bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
                             arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
        plt.show()


if __name__ == '__main__':
    draw_hbonds_profiles()
