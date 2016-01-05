#!/usr/bin/env python3.5
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import pickle
import os

import MyPDB as mp

work_dir = '/home/labs/fleishman/jonathaw/pseudo_symmetry/'


class Entry:
    def __init__(self, pdb, chain, segs, axis: mp.XYZ):
        self.pdb = pdb
        self.chain = chain
        self.symm_segs = segs
        self.interface = {}
        self.coverage = {}
        self.symm_axis = axis

    def __repr__(self):
        return "%s.%s %r %r %r %r" % (self.pdb, self.chain, self.symm_segs, self.interface, self.coverage,
                                      self.symm_axis)

    def set_interface(self, interface: dict):
        """
        :param interface: {chain: [1, 2, 3...]} interface residues
        :return: None
        """
        self.interface = interface

    def calc_coverage(self):
        """
        :return: calculates the fraction of interface residues that are in the symmetry
        """
        self.coverage = {}
        rng = [a for seg in self.symm_segs for a in range(seg[0], seg[1]+1)]
        for ch, interface in self.interface.items():
            self.coverage[ch] = len([a for a in interface if a in rng]) / len(interface)


def get_ce_symm_file():
    import urllib.request
    url = 'http://source.rcsb.org/jfatcatserver/downloadCensus.jsp?type=xml'
    print('downloading data from %s' % url)
    data = urllib.request.urlopen(url).read().decode('utf-8')
    with open(work_dir+'ce_symm.xml', 'w+') as fout:
        fout.write(data)


def parse_symm_db(args):
    """
    :param args:
    :return: list of Entry instances that are significant in the CE database
    """
    results = []
    with open(work_dir+'ce_symm.xml', 'r') as fin:
        xml = fin.read().split('data')
    ABC = 'QWERTYUIOPASDFGHJKLZXCVBNM'
    for data in xml:
        axis = {}
        significant = False
        for l in data.split('\n'):
            if 'isSig' in l:
                significant = l.split('>')[1].split('<')[0] == 'true'
                if not significant:
                    break
            if '<x>' in l:
                axis['x'] = float(l.split('>')[1].split('<')[0])
            if '<y>' in l:
                axis['y'] = float(l.split('>')[1].split('<')[0])
            if '<z>' in l:
                axis['z'] = float(l.split('>')[1].split('<')[0])
            if 'protodomain' in l:
                try:
                    s = l.split('>')[1].split('<')[0].split(',')
                    pdb = s[0][:4]
                    chain = s[0][5]
                    if len(s) < 2:
                        continue
                    segs = []
                    for s_ in s:
                        if '_-' not in s_:
                            start = int(''.join([a for a in s_.split('_')[1].split('-')[0] if a not in ABC]))
                            end = int(''.join([a for a in s_.split('_')[1].split('-')[1] if a not in ABC]))
                        if '_-' in s_:
                            start = int(''.join([a for a in s_.split('_-')[1].split('-')[0] if a not in ABC]))
                            end = int(''.join([a for a in s_.split('_-')[1].split('-')[1] if a not in ABC]))
                        segs.append((start, end))
                    results.append(Entry(pdb, chain, segs, mp.XYZ(axis['x'], axis['y'], axis['z'])))
                except:
                    print('SKIPPING', s)
    return results


def get_pdb_list(args):
    symm_db = parse_symm_db(args)
    for e in symm_db:
        print(e.pdb, e.chain)


def parse_pdb_list() -> dict:
    """
    :return: {pdb: chain} from the original pdb_name_chain list
    """
    with open(work_dir+'pdb_list.txt', 'r') as fin:
        return {l.split()[0]: l.split()[1] for l in fin.read().split('\n') if len(l.split()) > 1}


def find_pdb(pdb_name: str) -> str:
    dirs = ['dir_pdb_list_unique.txt.%i' % i for i in range(16)]
    for dir in dirs:
        dir_pdbs = [a for a in os.listdir(work_dir+'database/'+dir) if a[-4:] == '.pdb']
        if pdb_name+'.pdb' in dir_pdbs:
            return work_dir+'database/'+dir+'/'+pdb_name+'.pdb'
    return None


def find_interface_analysis(pdb_name: str, chain: str) -> str:
    dirs = ['dir_pdb_list_unique.txt.%i' % i for i in range(16)]
    for dir in dirs:
        dir_txts = [a for a in os.listdir(work_dir+'database/'+dir) if a[-4:] == '.txt']
        if '%s_%s_result.txt' % (pdb_name, chain) in dir_txts:
            return work_dir+'database/'+dir+'/'+'%s_%s_result.txt' % (pdb_name, chain)
    return None


def pdb2interface(pdb_name: str, chain: str):
    """
    :param pdb: a pdb name
    :param chain: a chain
    :return: creates a list of residues for any interface of chain with other chains with seq isentity < 0.95
    """
    results = {}
    pdb_path = find_pdb(pdb_name)
    pdb = mp.parse_PDB(file_in=pdb_path, name=pdb_name, with_non_residue=False)
    if len(list(pdb.seqs.keys())) < 2:
        return None
    wanted_seq = pdb.seqs[chain]
    for ch_id, ch_seq in pdb.seqs.items():
        if ch_id != chain:
            score = wanted_seq.align(ch_seq)
            if wanted_seq.aligned_identity(ch_seq) < 0.95:
                interface = mp.interface_residues(pdb[chain], pdb[ch_id])
                results[ch_id] = sorted([a.res_num for a in interface])
    out_path = '/'.join(pdb_path.split('/')[:-1])+'/'+pdb_name+'_'+chain+'_result.txt'
    with open(out_path, 'w+') as fout:
        fout.write('name %s\n' % pdb_name)
        fout.write('chain %s\n' % chain)
        for ch, inter in results.items():
            fout.write('%s %s\n' % (ch, ' '.join([str(a) for a in inter])))
            fout.write('seq %s %s\n' % (ch, pdb.seqs[ch].get_seq()))


def parse_unwanted_list(list_name) -> list:
    """
    :param list_name: name of unwanted list
    :return:
    """
    with open(work_dir+'unwanted_lists/'+list_name+'.txt', 'r') as fin:
        return fin.read().split()


def parse_interface_analysis(pdb: str, chain: str) -> dict:
    """
    :param pdb: pdb name
    :param chain: chain name
    :return: {chian: [1, 2, 3, ...]} interface residues
    """
    result = {}
    analisys_file = find_interface_analysis(pdb, chain)
    if analisys_file is not None:
        with open(analisys_file, 'r') as fin:
            cont = fin.read().split('\n')
        if len(cont) < 3:
            return None
        for l in cont:
            s = l.split()
            if len(s) < 2:
                continue
            if s[0] == 'name':
                if s[1] != pdb:
                    print('SKIPPING', cont)
                    return None
            if s[0] == 'chain':
                if s[1].upper() != chain.upper():
                    print('SKIPPING', cont)
                    return None
            if len(s[0]) == 1:
                result[s[0]] = [int(a) for a in s[1:]]
        return result
    # else:
    #     print('no interface analysis found')


def analyse_interface(args):
    """
    :param args:
    :return: {pdb: {chain: [1, 2, 3]}} interface residues
    """
    results = {}

    all_unwanted = []
    [all_unwanted.extend(parse_unwanted_list(k)) for k in ['antibodies']]
    pdb_list = parse_pdb_list()

    i = 0
    for pdb, chain in pdb_list.items():
        if pdb not in all_unwanted:
            i += 1
            interface = parse_interface_analysis(pdb, chain)
            if interface is None:
                continue
            if len(list(interface.keys())) == 0:
                continue
            # print(pdb, chain, interface)
            results[pdb] = {chain: interface}
    return results


def integrate_symm_interface(args):
    symm_list = parse_symm_db(args)
    interface_dict = analyse_interface(args)
    missing_interface = []
    for e in symm_list:
        if e.pdb in interface_dict.keys():
            if e.chain in interface_dict[e.pdb].keys():
                e.set_interface(interface_dict[e.pdb][e.chain])
                print('found', e)
        else:
            missing_interface.append(e)
            # print('missing interface info', e)
    print('a total of %i entries had no interface data out of a total of %i entries' %
          (len(missing_interface), len(symm_list)))
    with open(work_dir+'symm_interface_analysis.obj', 'wb') as fout:
        pickle.dump(symm_list, fout)


def draw_coverage_distribution(e_list: list):
    plt.figure()
    coverages = [c for e in e_list for c in e.coverage.values()]
    plt.hist(coverages, bins=50)
    plt.xlabel('symmetric interface residues')
    plt.ylabel('#interfaces')
    plt.show()


def integrate_analysis(args):
    with open(work_dir+'symm_interface_analysis.obj', 'rb') as fin:
        entry_list = pickle.load(fin)
    i = 0
    e_list = []
    pass_threshold = []
    for e in entry_list:
        if e.interface != []:
            e.calc_coverage()
            e_list.append(e)
            # print(e)
            i += 1
            for ch, cover in e.coverage.items():
                if cover >= 0.95:
                    if e.pdb in ['1ohz', '2ccl', '2b59', '2ozn', '2vn5', '2y3n', '3ul4', '4fl4', '4fl5', '4dh2', '4uyp',
                                 '5new', '4uyq', '3kcp', '2vn6', '4iu2']:
                        print('\nAAAAAAAAAAAAAA\n%s.%s\nBBBBBBBBBBBB\n' % (e.pdb, e.chain))
                        if test_symm_axis_interface(e, ch):
                            print('THRESHOLD #$@$@#%@#$%')
                            pass_threshold.append((e, ch))
        else:
            entry_list.remove(e)
    # draw_coverage_distribution(e_list)
    print('we have %i matches' % i)
    print('%i passed the threshold' % len(pass_threshold))
    with open('symm_interface_passed.obj', 'wb') as fin:
        pickle.dump(pass_threshold, fin)


def test_symm_axis_interface(e: Entry, ch: str) -> bool:
    print(e)
    inter_resi = e.interface[ch]
    print(e.symm_segs)
    symm_resi = [a for seg in e.symm_segs for a in range(seg[0], seg[1]+1)]
    print('symm_resi', symm_resi)
    print(inter_resi)
    pdb = mp.parse_PDB(find_pdb(e.pdb), with_non_residue=False)
    symm_com = mp.com_residues(chain=pdb[e.chain], residues=symm_resi)
    print(symm_resi, symm_com)
    pdb.translate_xyz(mp.XYZ(0, 0, 0)-symm_com)
    print('pdb translated')
    print('writing pdb %s_%s_translated.pdb' % (e.pdb, e.chain))
    mp.write_PDB('%s_%s_translated.pdb' % (e.pdb, e.chain), pdb)
    print('now symm com is at', mp.com_residues(chain=pdb[e.chain], residues=symm_resi))
    print('symm axis XYZ %r' % (e.symm_axis))
    print_pymol_select(e)


def print_pymol_select(e: Entry) -> None:
    print(e.pdb, e.chain)
    print('select symm, chain %s and resi %s' % (e.chain, '+'.join([str(a) for seg in e.symm_segs for a in
                                                                    range(seg[0], seg[1]+1)])))
    print('select interface, chain %s and resi %s' % (e.chain, '+'.join([str(r) for ch, inter in e.interface.items()
                                                                         for r in inter])))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode')
    parser.add_argument('-pdb')
    parser.add_argument('-chain')

    args = vars(parser.parse_args())

    if args['mode'] == 'parse_symm_db':
        parse_symm_db(args)

    elif args['mode'] == 'get_ce_symm_file':
        get_ce_symm_file()

    elif args['mode'] == 'get_pdb_list':
        get_pdb_list(args)

    elif args['mode'] == 'pdb2interface':
        pdb2interface(args['pdb'], args['chain'])

    elif args['mode'] == 'test':
        print(find_pdb(args['pdb']))

    elif args['mode'] == 'analyse_interface':
        print(analyse_interface(args))

    elif args['mode'] == 'integrate_symm_interface':
        integrate_symm_interface(args)

    elif args['mode'] == 'integrate_analysis':
        integrate_analysis(args)

    else:
        print('mode unknown')


if __name__ == '__main__':
    main()