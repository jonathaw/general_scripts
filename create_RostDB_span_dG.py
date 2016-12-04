#!/usr/bin/env python3.5
"""
scripts for working with the Rost data set structurally
"""
import argparse
import matplotlib.pyplot as plt
import numpy as np
from retrive_natural_TMs_scores import parse_rost_db, find_topo, create_AddMembrane_xml
import shutil
import os
import re
from get_pdb import download_pdb
from MyPDB import parse_PDB, MyPDB, write_PDB, extract_seq
from AASeq import AASeq
import sys
from Logger import Logger


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode')
    args = vars(parser.parse_args())

    if args['mode'] == 'setup':
        """
        iterate the Rost data set and make a dir with each with:
        1: download each PDB, isolate the correct chain
        2: read the spans from Rost PDBTM, create a XML subroutine for add membrane
        3: create flags file to run an embedding xml
        """
        setup_db(args)

    elif args['mode'] == 'dg_dist':
        """
        go over all results from running jobs that run the flag files made by setup_db.
        gather the ∆G_insertion for every span, and show a hist and some statistics
        """
        analyse_span_dg_distribution(args)

    else:
        print('no mode found')


def analyse_span_dg_distribution(args):
    rost_db = parse_rost_db()
    all_dgs = []
    single_span_dgs, multi_span_dgs = [], []
    for k in rost_db.keys():
        if k in ['q9u6b8', 'q03nm0', 'p32851', 'p0ab98', 'p0abj3']:
            print('skipping %s becuause its wierd' % k)
            continue
        all_files = os.listdir('%s' % k)
        out_file = [a for a in all_files if 'out.' in a][0]
        this_entry_dgs = []
        for l in open('%s/%s' % (k, out_file), 'r'):
            if 'core.scoring.membrane.MPSpanInsertionEnergy' in l:
                if 'for span' in l and 'skipping' not in l:
                    # print(l)
                    all_dgs.append(float(l.split()[5]))
                    this_entry_dgs.append(float(l.split()[5]))
        if len(this_entry_dgs) == 1:
            single_span_dgs += this_entry_dgs
        else:
            multi_span_dgs += this_entry_dgs
    for lst, name in zip([single_span_dgs, multi_span_dgs, all_dgs], ['single', 'multi', 'all']):
        print('%s mean %.2f' % (name, np.mean(lst)))
        print('%s median %.2f' % (name, np.median(lst)))
        print('%s std %.2f' % (name, np.std(lst)))
        print('%s entries %i' % (name, len(lst)))
    plt.hist(all_dgs, bins=50, color='red', alpha=0.5, label='all')
    plt.hist(single_span_dgs, bins=50, stacked=True, color='b', alpha=0.5, label='single')
    plt.hist(multi_span_dgs, bins=50, stacked=True, color='g', alpha=0.5, label='multi')
    plt.legend()
    plt.title('span ∆G destributions')
    plt.xlabel('∆G')
    plt.ylabel('count')
    plt.show()


def setup_db(args):
    rost_db = parse_rost_db()
    failed = []
    logger = Logger('./db_setup.log')
    for k, v in rost_db.items():
        # if k != 'q9u6b8': continue
        logger.create_header('working on %s' % k)
        logger.log('seq: %s' % v['seq'])
        logger.log('pdb: %s' % v['pdb'])
        logger.log('chain: %s' % v['chain'])
        logger.log('ts: %s' % v['ts'])
        os.mkdir(k)
        os.chdir(k)

        # get pdb and extract chain
        download_pdb({'name': v['pdb'], 'path': './'})
        empty_pdb = MyPDB(name=v['pdb'])
        pdb = parse_PDB('pdb%s.ent' % v['pdb'])
        chain = pdb.chains[v['chain']]
        empty_pdb.add_chain(chain)
        write_PDB('%s_%s.pdb' % (k, v['chain']), empty_pdb)
        pdb_seq = extract_seq(empty_pdb)
        rdb_seq = AASeq(v['seq'])
        score, start, end = pdb_seq[v['chain']].align(rdb_seq)
        logger.log('pdb seq: %s' % pdb_seq[v['chain']].aligned)
        logger.log('rst seq: %s' % rdb_seq.aligned)


        # get spans and print xml
        spans = find_topo(v['ts'])

        new_spans = []
        for sp in spans:
            start = pdb_seq[v['chain']].aligned_position_at_non_aligned(sp[0]) + 1
            end = pdb_seq[v['chain']].aligned_position_at_non_aligned(sp[1]) + 1
            logger.log('span %i->%i %s moving to %i->%i' %(sp[0], sp[1], sp[2], start, end))
            new_spans.append([start, end, sp[2]])
        create_AddMembrane_xml(new_spans, '%s_AddMembrane.xml' % v['pdb'])

        # create flags file
        with open('embed.flags', 'w+') as fout:
            fout.write('-parser:protocol /home/labs/fleishman/jonathaw/elazaridis/protocols/embed_in_membrane.xml\n')
            fout.write('-s %s\n' % '%s_%s.pdb' % (k, v['chain']))
            fout.write('-parser:script_vars add_memb_xml=%s\n' % '%s_AddMembrane.xml' % v['pdb'])
            fout.write('-overwrite\n')
            fout.write('-score::elec_memb_sig_die\n')
            fout.write('-corrections::beta_nov15\n')
            fout.write('-score::memb_fa_sol\n')
        os.chdir('../')


if __name__ == '__main__':
    main()
