import subprocess
import re
import operator
import sys
import os
import argparse
import shutil
from rosetta_score_files import score_passes_thresholds
from DoCohJobMaker_options import design_prediction_by_list
from sequence_funcitons import seq_from_pdb


def main():
    args = parse_args()
    best_scores_dict = best_scores(args.score_file)
    best_scores_sort = sorted(best_scores_dict.items(), key=operator.itemgetter(1))
    best_desc = [best_scores_sort[x][0] for x in range(50)]
    new_names = new_names_maker(best_desc, args.short_name, args.names_list)
    extract_silent_files(best_desc, args.silent_file)
    # make_seq_files(new_names, 'coh', cwd + args.short_name + '_coh.fasta')
    # make_seq_files(new_names, 'doc', cwd + args.short_name + '_doc.fasta')
    # design_prediction_by_list(args)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-score_file', type=str)
    parser.add_argument('-silent_file', type=str)
    parser.add_argument('-design_num', default=50, type=int)
    parser.add_argument('-short_name', type=str)
    parser.add_argument('-names_list', default=cwd + '/top_names', type=str)
    parser.add_argument('-coh_name', type=str)
    parser.add_argument('-doc_name', type=str)

    return parser.parse_args()


def best_scores(score_file):
    thresholds = {'ddg': -20, 'sasa': 1300, 'pack': 0.6, 'shape': 0.5, 'buried': 0.0}
    passing_scores = {}
    with open(score_file) as f:
        for score_line in f:
            line_split = score_line.split()
            if line_split[0] == 'SCORE:' and line_split[1] == 'score':
                score_fields = {'ddg': line_split.index('a_ddg'), 'sasa': line_split.index('a_sasa'),
                                'pack': line_split.index('a_packstat'), 'shape': line_split.index('a_shape'),
                                'buried': line_split.index('a_buried_2'),
                                'description': line_split.index('description')}
            elif line_split[0] == 'SCORE:' and line_split[1] != 'score':
                if score_passes_thresholds(score_line, score_fields, thresholds):
                    passing_scores[line_split[score_fields['description']]] = float(line_split[score_fields['ddg']])
    return passing_scores


def new_names_maker(original_names, short_name, file_name):
    original_names.sort()
    new_names = []
    with open(file_name, 'wr+') as f:
        i = 11
        for name in original_names:
            f.write(name + '\t' + short_name + str(i) + '\n')
            new_names.append(short_name + str(i))
            i += 1
    return new_names


def extract_silent_files(names, silent_file):
    with open(cwd+'/pdb_temp_list', 'wr+') as f:
        [f.writelines(name+'\n') for name in names]
        subprocess.call(['/home/labs/fleishman/jonathaw/Rosetta/main/source/bin/extract_pdbs.default.linuxgccrelease',
                         '-in:file:silent *.cmb -in:file:tags `cat pdb_temp_list`',
                         '-database /home/labs/fleishman/jonathaw/Rosetta/main/database/'])


def make_seq_files(names_list, coh_doc, out_file_name):
    chain_type = {'coh': 'A', 'doc': 'B'}
    with open(out_file_name, 'wr+') as f:
        print names_list
        for name in names_list:
            f.write('>' + name + '.' + chain_type[coh_doc] + '\n')
            f.write(seq_from_pdb(name+'.pdb', chain_type[coh_doc]))


if __name__ == '__main__':
    cwd = os.getcwd()
    main()