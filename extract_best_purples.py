#!/usr/bin/env python
from rosetta_score_files import score_passes_thresholds


def main(args):
    import os
    import re
    import subprocess
    score_file = [x for x in os.listdir('./') if re.match('all.*\.score', x)][0]
    cmb_file = os.getcwd() +[x for x in os.listdir('./') if re.match('.*\.cmb', x)][0]
    print score_file, cmb_file
    thresholds = {'ddg': -20.0, 'sasa': 1300.0, 'pack': 0.6, 'shape': 0.5, 'buried': 1.0}
    filter_fields = {}
    scores = []
    with open(score_file, 'r') as f:
        for file_line in f:
            line_split = file_line.split()
            if 'SCORE:' in line_split:
                if line_split[1] == 'score':
                    filter_fields['ddg'] = line_split.index('a_ddg')
                    filter_fields['sasa'] = line_split.index('a_sasa')
                    filter_fields['pack'] = line_split.index('a_packstat')
                    filter_fields['shape'] = line_split.index('a_shape')
                    filter_fields['buried'] = line_split.index('a_buried_2')
                    filter_fields['description'] = line_split.index('description')
                else:
                    if score_passes_thresholds(file_line, filter_fields, thresholds):
                        scores.append({float(line_split[filter_fields['ddg']]): str(line_split[filter_fields['description']])})
    scores = sorted(scores, key=lambda x: x.keys()[0])

    # with open(os.getcwd()+'/purples_extract', 'wr+') as f:
    #     [f.writelines(name.values()[0]+'\n') for name in scores]
    #     subprocess.call(['/home/labs/fleishman/jonathaw/Rosetta/main/source/bin/extract_pdbs.default.linuxgccrelease',
    #                      '-in:file:silent %s -in:file:tags `cat %s/purples_extract`' % (cmb_file, os.getcwd()),
    #                      '-database /home/labs/fleishman/jonathaw/Rosetta/main/database/'])
    with open('purples_names_extract', 'wr+') as f:
        i = 1
        for purple in scores:
            print 'writing %s to file' % purple.values()[0]
            f.writelines(purple.values()[0]+'\n')
            if i == args['num_structures']:
                break
    subprocess.call(['/home/labs/fleishman/jonathaw/Rosetta/main/source/bin/extract_pdbs.default.linuxgccrelease',
                         '-in:file:silent %s -in:file:tags `cat purples_names_extract`' % cmb_file,
                         '-database /home/labs/fleishman/jonathaw/Rosetta/main/database/'])
    # subprocess.call('/home/labs/fleishman/jonathaw/Rosetta/main/source/bin/extract_pdbs.default.linuxgccrelease -in:file:silent *cmb -in:file:tags `cat purples_names_extract` -database /home/labs/fleishman/jonathaw/Rosetta/main/database/')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-num_structures', default=10, type=int)
    args = vars(parser.parse_args())
    main(args)