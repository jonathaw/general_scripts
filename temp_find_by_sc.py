#!/usr/bin/env python3.5
import os, shutil
from DoCohResultProcessor import generate_run_filters, all_who_pass_run_filters
from RosettaFilter import score2dict


topath = '/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/recliques_4Nov/clique_6_pdbs/mini_diagonal_11Nov/'
run_filters = generate_run_filters(args={'ddg': 24.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6, 'buried_2': 3,
                                         'hbonds': 10})
all_docs = [a for a in os.listdir(topath) if 'doc_' in a]
fout = open(topath+'minidiagonal_full_names.txt', 'w')

for doc in all_docs:
    all_dirs = os.listdir(topath+doc)
    for dir in all_dirs:
        try:
            sc_file = [a for a in os.listdir(topath+doc+'/'+dir) if a[-3:] == '.sc']
            if sc_file:
                sc_dict = score2dict(topath+doc+'/'+dir+'/'+sc_file[0])
                passed, failed = all_who_pass_run_filters({}, sc_dict, run_filters)
                if len(passed) > 5:
                    fout.write('%s\t%i\n' % (dir, len(passed)))
                    shutil.copy(topath+doc+'/'+dir+'/'+list(sc_dict.keys())[0]+'.pdb.gz', topath+'minidiagonal_pdbs')
        except:
            print('no folder', dir)
fout.close()

# def analyse_minidiagonal(args):
#     with open('../minidiagonal.txt', 'w') as fout:
#         run_filters = generate_run_filters(args={'ddg': 24.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6, 'buried_2': 3,
#                                                  'hbonds': 10})
#         counter = 0
#         score_files = [a for a in os.listdir('./') if a[-3:] == '.sc']
#         for sc in score_files:
#             score_dict = score2dict(sc)
#             passed, failed = all_who_pass_run_filters(args, score_dict, run_filters)
#             if len(passed) > 5:
#                 fout.write('%s\t%i\n' % (sc, len(passed)))
#                 counter += 1
#     print('%i passed minidiagonal' % counter)


"""
import os, sys

topath = '/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/recliques_4Nov/clique_6_pdbs/mini_diagonal_11Nov/'
doc_list = ['1ohz', '2vn5', '2y3n', '3ul4', '4dh2', '4fl4', '4fl5', '4uyp', '5new']
fold_lists = {doc: {fold: [a for a in os.listdir(topath+'doc_'+doc+'/'+fold) if a[-3:] == '.sc']
                    for fold in os.listdir(topath+'doc_'+doc) if fold[:2] == 'dz'} for doc in doc_list[:1]}

# for k, v in fold_lists.items():
#     print(k)
#     for k1, v1 in v.items():
#         print(k1, v1)
#         break

not_found = []
for n in open(topath+'minidiagonal.txt', 'r'):
    name = n.split()[0]
    s = name.split('_dz_')[0]
    doc_name = s.split('1ohz_A_')[1].split('_')[0]
    found = False
    for entry in fold_lists[doc_name].keys():
        if s in entry:
            if name in fold_lists[doc_name][entry]:
                print(entry)
                if found:
                    print(name, entry, s)
                    sys.exit()
                found = True
    if not found:
        not_found.append(s)
print("did not find these:\n", not_found)
"""