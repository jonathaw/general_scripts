#!/usr/bin/env python3.5
import shutil
import os
import subprocess
import re
import random
import string
import sys
from AASeq import AASeq
import seq_funcs
import MyPDB

path_seqs_rb = '/home/labs/fleishman/jonathaw/data/pdbs_final/chosen_chains_rearranged/'
rb_path = '/home/labs/fleishman/jonathaw/data/RBouts/'
COH_TEMPLATES = ['1OHZ', '2B59', '2CCL', '2OZN', '2VN5', '2VN6', '2Y3N', '3KCP', '3UL4', '4DH2', '4IU2', '4UYP', '4UYQ',
                 '4FL4', '1ANU', '1AOH', '1G1K', '1QZN', '1TYJ', '2BM3', '2JH2', '2VO8', '2W1N', '2W5F', '2XBT', '2XDH',
                 '2ZF9', '3BWZ', '3FNK', '3L8Q', '4JO5', '4UMS', '4FL5', '5NEW']
DOC_TEMPLATES = ['1OHZ', '2B59', '2CCL', '2OZN', '2VN5', '2VN6', '2Y3N', '3KCP', '3UL4', '4DH2', '4IU2', '4UYP', '4UYQ',
                 '4FL4', '4FL5', '5NEW']
real_names = {'1oha': '1ohz', '2b5a': '2b59', '2cca': '2ccl', '2oza': '2ozn', '2v5a': '2vn5', '2v6a': '2vn6',
              '2y3a': '2y3n', '3kca': '3kcp', '3ula': '3ul4', '4dha': '4dh2', '4iua': '4iu2', '4upa': '4uyp',
              '4uqa': '4uyq', '4f4a': '4fl4', '4fla': '4fl4', '4f5a': '4fl5', '5nea': '5new'}
COH_CLUSTER = {'1ohz': {'temp': ['1ohz', '1anu', '2ccl', '4fl4', '1aoh', '3ul4'],
                        'strands': ['GIANCDFVFRY', 'FDTAIYPDR', 'KIIVFLFAEDSG', 'KSFDTAIYP']},
               '2vn5': {'temp': ['2vn5', '2vn6', '1g1k'], 'strands': ['VNFSSSASN', 'GTISFLFL', 'NVGTCNFYLGY',
                                                                      'FKDGGAFGDGTMSKIASV']},
               '2b59': {'temp': ['2b59', '3kcp', '4fl5', '2bm3'],
                        'strands': ['GPIQIADND', 'GILNFALAYS', 'NFAGFQVNIVY', 'FQDTLSMPGAISGTQLFDWDGEVITGYE']},
               '1qzn': {'temp': ['1qzn', '3l8q', '3fnk', '3bwz'],
                        'strands': ['FNLRQVADN', 'KGILNFSKAYV', 'NFSGCQLNMKY', 'FEDTTSVPNAIDGTVLFDWNGDRI']},
               '2y3n': {'temp': ['2y3n', '1tyj'],
                        'strands': ['NATDMSKHN', 'GVLNFGRLYMNL', 'NFSGYQFNIKY', 'TFENGSSMNNAVDGTMLFDWDGNMY']},
               '5new': {'temp': ['5new'],
                        'strands': ['SDQTASLLPTFDVSIQ', 'YSSVIWSTA', 'QGIDFAVT',
                                    'IKLEAVKRETYVGSGTDNSSISAGYSANDKAVKYTV']},
               '2w5f': {'temp': ['2w5f'],
                        'strands': ['WTEISAKYK', 'EYKYSVFVK', 'FIFDDVTITRK', 'YEVVHD']},
               '4jo5': {'temp': ['4jo5'], 'strands': ['VKGTFVKMSSST', 'ADTYLEISF', 'LTLRYYY', 'VEWDQVTAYLNGVLVWGKE']},
               '2zf9': {'temp': ['2zf9'], 'strands': ['LKKAEADT', 'SFFTATG', 'ATGLHIQF', 'KYDVQVAYQSRTTNEDLFTNVKKD']},
               '4iu2': {'temp': ['4iu2'], 'strands': ['LAKAENNG', 'GVFVASG', 'TTGYHIYW', 'VYPIDVAYQWDPSKGDLFTDNKDSA']},
               '4dh2': {'temp': ['4dh2'], 'strands': ['NNFDYNIV', 'EIVFLFD', 'SFNLSLY', 'LITFGESNFC']},
               '4ums': {'temp': ['4ums'], 'strands': ['NFSNNKSE', 'KISFLFN', 'NCDFKLV', 'RKDLVGSFSGLKDNKMTSIG']},
               '4uyp': {'temp': ['4uyp', '4uyq'], 'strands': ['SDFTTYYN', 'FASMTFE', 'ALSFRTN', 'TNSAYTSFYYSGTD']},
               '2ozn': {'temp': ['2ozn', '2jh2', '2vo8', '2w1n'],
                        'strands': ['VFVNAKKIE', 'VRVLVSS', 'AYDFTLN', 'LSVTNSSVGDGEGLVHEIAGT']},
               '2xdh': {'temp': ['2xdh'], 'strands': ['LFDYQVE', 'QIKVGIADS', 'SINLILS', 'TLQGIEIYDIDGNS']}}
# no 2xbt. it's a weird cohesin.
DOC_CLUSTER = {'2b59': {'temp': ['2b59', '2ozn', '3kcp'], 'strands': ['DIVKDNSINLLDVAEVIRCF', 'NRNGAINMQDIMIVHKH']},
               '4iu2': {'temp': ['4iu2'], 'strands': ['DTDLNNIVDGRDATATLTYYAATS', 'DGRDASSILTFYTKSSV']},
               '4fl5': {'temp': ['4fl5'], 'strands': ['NSINLLDVAEVIRC', 'NRNGAINMQDIMIVHK']},
               '1ohz': {'temp': ['1ohz', '2ccl', '2y3n', '3ul4', '4dh2', '4fl4', '4uyp', '4uyq', '5new'],
                        'strands': ['DVNGDGTINSTDLTMLKRSVL', 'DVDKNGSINSTDVLLLSRYLL']},
               '2vn5': {'temp': ['2vn5', '2vn6'], 'strands': ['DYNNDGNVDALDFAGLKKYIM', 'DVNLDNEVNSTDLAILKKYLLGM']}}

RB_NAMES = {1: 'RB_1ohz_1.db', 2: 'RB_1ohz_2.db', 3: 'RB_2b59_1.db', 4: 'RB_2b59_2.db', 5: 'RB_2ccl_1.db',
            6: 'RB_2ccl_2.db', 7: 'RB_2ozn_1.db', 8: 'RB_2ozn_2.db', 9: 'RB_2vn5_1.db', 10: 'RB_2vn5_2.db',
            11: 'RB_2vn6_1.db', 12: 'RB_2vn6_2.db', 13: 'RB_2y3n_1.db', 14: 'RB_2y3n_2.db', 15: 'RB_3kcp_1.db',
            16: 'RB_3kcp_2.db', 17: 'RB_3ul4_1.db', 18: 'RB_3ul4_2.db', 19: 'RB_4dh2_1.db', 20: 'RB_4dh2_2.db',
            21: 'RB_4fl4_1.db', 22: 'RB_4fl4_2.db', 23: 'RB_4iu2_1.db', 24: 'RB_4iu2_2.db', 25: 'RB_4uyp_1.db',
            26: 'RB_4uyp_2.db', 27: 'RB_4uyq_1.db', 28: 'RB_4uyq_2.db', 29: 'RB_4fl5_1.db', 30: 'RB_4fl5_2.db',
            31: 'RB_5new_1.db', 32: 'RB_5new_2.db'}

doc_design_pos = {'1ohz': {'seq': AASeq('GDVNGDGTINSTDLTMLKRSVLRAITLTDDAKARADVDKNGSINSTDVLLLSRYLL', '1ohz.B'),
                           'design': {11: 'S', 12: 'T', 15: 'T', 18: 'K', 19: 'R', 21: 'V', 22: 'L', 23: 'R', 48: 'V',
                                      49: 'L', 52: 'S', 53: 'R', 56: 'L', 46: 'T', 44: 'N', 50: 'L', 14: 'L', 45: 'S',
                                      24: 'A'}},
                  '2b59': {'seq': AASeq(
                      'GYKVSGYILPDFSFDATVAPLVKAGFKVEIVGTELYAVTDANGYFEITGVPANASGYTLKISRATYLDRVIANVVVTGDTSVSTSQAPIMMWVGDIVKDNSINLLDVAEVIRCFNATKGSANYVEELDINRNGAINMQDIMIVHKHFGATSSDYDA',
                      '2b59.B'),
                           'design': {104: 'L', 111: 'I', 114: 'F', 115: 'N', 136: 'N', 137: 'M', 140: 'I', 141: 'M',
                                      144: 'H', 147: 'F', 13: 'S', 112: 'R', 117: 'T', 120: 'S', 122: 'N', 138: 'Q',
                                      145: 'K', 107: 'V'}},
                  '2ozn': {'seq': AASeq(
                      'DKTNLGELINQGKSLLDESVEGFNVGEYHKGAKDGLTVEINKAEEVFNKEDATEEEINLAKESLEGAIARFNSLLIEESTGDFNGNGKIDIGDLAMVSKNIGSTTNTSLDLNKDGSIDEYEISFINHRILN',
                      '2ozn.B'),
                           'design': {23: 'F', 30: 'K', 91: 'I', 94: 'L', 95: 'A', 98: 'S', 99: 'K', 101: 'I', 102: 'G',
                                      104: 'T', 116: 'S', 119: 'E', 120: 'Y', 122: 'I', 123: 'S', 126: 'N', 127: 'H',
                                      130: 'L', 131: 'N'}}
                  }

for doc, data in doc_design_pos.items():
    for i, aa in data['design'].items():
        assert data['seq'][i] == aa, 'in %s found %s in seq at %i but %s in design' % (doc, data['seq'][i], i, aa)


def RBIn_SASA_test(args):
    args['protocol'] = 'RBIn_SASA_test.xml'
    args['flags_file'] = 'docoh_modeling.flags'
    shutil.copyfile(path_seqs_rb + 'all_docs.fasta', os.getcwd() + '/all_doc.fasta')
    doc_fasta = read_multi_fastas('all_doc.fasta')

    coh_seq_1ohz = AASeq(
        'GVVVEIGKVTGSVGTTVEIPVYFRGVPSKGIANCDFVFRYDPNVLEIIGIDPGDIIVDPNPTKSFDTAIYPDRKIIVFLFAEDSGTGAYAITKDGVFAKIRATVKSSAPGYITFDEVGGFADNDLVEQKVSFIDGGVNVG')
    coh_name = '1ohz'
    cohs_dict = {'rb_coh': RBposCoh(coh_seq_1ohz.get_seq(), coh_name), 'coh_len': len(coh_seq_1ohz)}
    for doc in ['1ohz', '2vn5', '2y3n', '3ul4', '4dh2', '4uyp', '5new', '4fl4', '4fl5']:
        os.chdir(args['path'])
        os.mkdir(doc + '_RB_test')
        os.chdir(doc + '_RB_test')
        doc_rb_pos = RBposDoc1st(doc_fasta[doc]['seq'], doc)
        for rb in RB_NAMES.keys():
            job_args = {
                '-parser:protocol': '/home/labs/fleishman/jonathaw/protocols/%s' % args['protocol'],
                '@/home/labs/fleishman/jonathaw/protocols/flags/%s' % args['flags_file']: None,
                '-database': '/home/labs/fleishman/jonathaw/Rosetta/main/database/',
                '-parser:script_vars': 'RB_file=%s%s' % (rb_path, RB_NAMES[rb]),
                '-overwrite': None,
                '-s': '%s%s_A.pdb' % (path_seqs_rb, '1ohz'),
                '-out:suffix': '_%s_%i.pdb' % (doc, rb),
            }
            script_vars = {
                'doc_pdb=': '%s%s_B.pdb' % (path_seqs_rb, doc.lower()),
                'coh_rb_pos=': str(cohs_dict['rb_coh']) + 'A',
                'doc_rb_pos=': str(doc_rb_pos + cohs_dict['coh_len']) + 'B',
                'doc_start_res=': str(cohs_dict['coh_len'] + 1),
            }
            make_jobs(args, job_args, script_vars)


def minidiagonal(args):
    pdb = MyPDB.parse_PDB(args['pdb'], args['coh_name'])
    seqs = MyPDB.extract_seq(pdb)
    args['coh_seq'] = seqs['A'].get_seq
    args['doc_seq'] = seqs['B'].get_seq
    args['coh_name'] = args['pdb'].split('.pdb')[0]
    args['doc_name'] = args['pdb'].split('.pdb')[0]

    args['flags_file'] = 'docoh_modeling.flags'
    print("copying data files")
    shutil.copyfile(path_seqs_rb + 'all_cohs.fasta', os.getcwd() + '/all_coh.fasta')
    shutil.copyfile(path_seqs_rb + 'all_docs.fasta', os.getcwd() + '/all_doc.fasta')

    print("adding query seqs and aligning:")
    add_align(args, 'coh')
    add_align(args, 'doc')

    coh_fasta = read_multi_fastas(args['coh_aln_file'] + '_aln')
    doc_fasta = read_multi_fastas(args['doc_aln_file'] + '_aln')

    rb = int(args['pdb'].split('_')[4])
    coh_bb = args['pdb'].split('_')[1]
    doc_bb = args['pdb'].split('_')[3]

    cohs_passed = [coh_bb]
    docs_passed = [doc_bb]

    cohs_dict = {k: {'rb_coh': None, 'coh_len': None} for k in cohs_passed}
    docs_dict = {k: {'rb_doc': None} for k in docs_passed}
    for t_coh in cohs_passed:
        cohs_dict[t_coh]['rb_coh'] = RBposCoh(coh_fasta[t_coh]['seq'], t_coh)
        cohs_dict[t_coh]['coh_len'] = len(coh_fasta[t_coh]['seq'].replace('-', ''))
    for t_doc in docs_passed:
        docs_dict[t_doc]['rb_doc'] = RBposDoc1st(doc_fasta[t_doc]['seq'], t_doc)

    job_args = {'-parser:protocol': '/home/labs/fleishman/jonathaw/protocols/%s' % args['protocol'],
                '@/home/labs/fleishman/jonathaw/protocols/flags/%s' % args['flags_file']: None,
                '-database': '/home/labs/fleishman/jonathaw/Rosetta/main/database/',
                '-parser:script_vars': 'RB_file=%s%s' % (rb_path, RB_NAMES[rb]), '-overwrite': None,
                '-s': '%s%s_A.pdb' % (path_seqs_rb, coh_bb.lower()), '-out:suffix': '_%s_%i.pdb' % (doc_bb, rb)
                }
    job_args['-out:prefix'] = '%s_%s_on_' % (args['coh_name'], args['doc_name'])
    script_vars = {
        'doc_pdb=': '%s%s_B.pdb' % (path_seqs_rb, doc_bb.lower()),
        'coh_aln_file=': args['coh_aln_file'] + '_aln',
        'doc_aln_file=': args['doc_aln_file'] + '_aln',
        'coh_query_name=': args['coh_name'].lower()[:4],
        'doc_query_name=': args['doc_name'].lower()[:4],
        'coh_template_name=': coh_bb + '.A',
        'doc_template_name=': doc_bb + '.B',
        'coh_rb_pos=': str(cohs_dict[coh_bb]['rb_coh']) + 'A',
        'doc_rb_pos=': str(docs_dict[doc_bb]['rb_doc'] + cohs_dict[coh_bb]['coh_len']) + 'B',
        'doc_start_res=': str(cohs_dict[coh_bb]['coh_len'] + 1),
        # 'hot_spot=': str(cohs_dict[t_coh]['rb_coh']+2) + 'A',
    }  # hot_spot is coh_rb_pos +2...
    make_jobs(args, job_args, script_vars)


def predict_pair(args):
    args['flags_file'] = 'docoh_modeling.flags'
    print("copying data files")
    shutil.copyfile(path_seqs_rb + 'all_cohs.fasta', os.getcwd() + '/all_coh.fasta')
    shutil.copyfile(path_seqs_rb + 'all_docs.fasta', os.getcwd() + '/all_doc.fasta')

    print("adding query seqs and aligning:")
    add_align(args, 'coh')
    add_align(args, 'doc')

    # added some docs without actually testing, SA 2b59, 2ozn, 3kcp, 4uyq, 2ccl. 8Nov2015
    doc_rb_select = {'1ohz': [10, 11, 17, 18, 19, 1, 20, 22, 25, 27, 2, 32, 3, 5, 6, 8, 9],
                     '2b59': range(1, 33),
                     '2ccl': [10, 11, 17, 18, 19, 1, 20, 22, 25, 27, 2, 32, 3, 5, 6, 8, 9],
                     '2vn5': [10, 16, 17, 18, 19, 1, 21, 25, 26, 27, 28, 2, 30, 32, 4, 5, 6, 9],
                     '2vn6': [10, 16, 17, 18, 19, 1, 21, 25, 26, 27, 28, 2, 30, 32, 4, 5, 6, 9],
                     '2ozn': range(1, 33),
                     '2y3n': [11, 12, 13, 18, 19, 1, 29, 31, 32, 9],
                     '3kcp': range(1, 33),
                     '3ul4': [10, 11, 12, 13, 17, 1, 20, 26, 28, 2, 31, 32, 5, 6, 9],
                     '4dh2': [11, 12, 14, 16, 18, 19, 1, 26, 29, 2, 30, 31, 32, 3, 5, 9],
                     '4fl4': [10, 11, 12, 16, 18, 19, 1, 21, 22, 25, 26, 27, 2, 30, 31, 32, 4, 5, 6, 9],
                     '4fl5': [10, 19, 20, 21, 22, 26, 32],
                     '4uyp': [11, 19, 20, 25, 26, 27, 28, 2, 30, 9],
                     '4uyq': [11, 19, 20, 25, 26, 27, 28, 2, 30, 9],
                     '5new': [10, 11, 12, 14, 16, 17, 19, 1, 20, 22, 29, 2, 31, 32, 3, 6, 9]}

    coh_fasta = read_multi_fastas(args['coh_aln_file'] + '_aln')
    doc_fasta = read_multi_fastas(args['doc_aln_file'] + '_aln')
    # print(coh_fasta)
    cohs_passed = cohs_pass(args, coh_fasta)
    args['logger'].log('chosen cohesins: %s' % ', '.join(cohs_passed))
    docs_passed = docs_pass(args, doc_fasta)
    if False:
        docs_passed = [a.lower() for a in DOC_TEMPLATES]

    args['logger'].log('chosen dockerins: %s' % ', '.join(docs_passed))
    cohs_dict = {k: {'rb_coh': None, 'coh_len': None} for k in cohs_passed}
    docs_dict = {k: {'rb_doc': None} for k in docs_passed}
    for t_coh in cohs_passed:
        cohs_dict[t_coh]['rb_coh'] = RBposCoh(coh_fasta[t_coh]['seq'], t_coh)
        cohs_dict[t_coh]['coh_len'] = len(coh_fasta[t_coh]['seq'].replace('-', ''))
    for t_doc in docs_passed:
        docs_dict[t_doc]['rb_doc'] = RBposDoc1st(doc_fasta[t_doc]['seq'], t_doc)

    if args['mode'] == 'jack':
        cohs_names = list(cohs_dict.keys())
        for t_coh in cohs_names:
            if seq_iden(coh_fasta[t_coh]['seq'], coh_fasta[args['coh_name'].lower()[:4]]['seq']) >= 0.9:
                args['logger'].log('%s coh is rejected %f' %
                                   (t_coh, seq_iden(coh_fasta[t_coh]['seq'],
                                                    coh_fasta[args['coh_name'].lower()[:4]]['seq'])))
                cohs_dict.pop(t_coh)
                if t_coh.upper() in DOC_TEMPLATES:
                    RB_NAMES.pop([k for k, v in RB_NAMES.items() if t_coh in v][0])
        docs_names = list(docs_dict.keys())
        for t_doc in docs_names:
            if seq_iden(doc_fasta[t_doc]['seq'], doc_fasta[args['doc_name'].lower()[:4]]['seq']) >= 0.9:
                args['logger'].log('%s doc is rejected %f' %
                                   (t_doc, seq_iden(doc_fasta[t_doc]['seq'],
                                                    doc_fasta[args['doc_name'].lower()[:4]]['seq'])))
                docs_dict.pop(t_doc)
                if t_doc.upper() in DOC_TEMPLATES:
                    RB_NAMES.pop([k for k, v in RB_NAMES.items() if t_doc in v][0])
        args['logger'].log('seq identity left cohs %s' % ', '.join(cohs_dict.keys()))
        args['logger'].log('seq identity left docs %s' % ', '.join(docs_dict.keys()))
        args['logger'].log('seq identity left RBs %s' % ', '.join(RB_NAMES.values()))

    args['logger'].log('there will be %i jobs made' % (len(list(cohs_dict.keys())) * len(list(docs_dict.keys())) *
                                                       len(list(RB_NAMES.keys()))))
    for t_coh in cohs_dict.keys():
        for t_doc in docs_dict.keys():
            for rb in RB_NAMES.keys():
                if args['specific_RBs']:
                    # if chosen, will only create jobs for RBs that passed the test, used for diagonal
                    if rb not in doc_rb_select[t_doc]:
                        continue
                job_args = {
                    '-parser:protocol': '/home/labs/fleishman/jonathaw/protocols/%s' % args['protocol'],
                    '@/home/labs/fleishman/jonathaw/protocols/flags/%s' % args['flags_file']: None,
                    '-database': '/home/labs/fleishman/jonathaw/Rosetta/main/database/',
                    '-parser:script_vars': 'RB_file=%s%s' % (rb_path, RB_NAMES[rb]),
                    '-overwrite': None,
                    '-s': '%s%s_A.pdb' % (path_seqs_rb, t_coh.lower()),
                    '-out:suffix': '_%s_%i.pdb' % (t_doc, rb),
                }
                if args['coh_name'] == args['doc_name']:
                    job_args['-out:prefix'] = '%s_on_' % args['coh_name']
                else:
                    job_args['-out:prefix'] = '%s_%s_on_' % (args['coh_name'], args['doc_name'])
                if args['mode'] == 'jack':
                    job_args['-native'] = '%s%s_AB.pdb' % (path_seqs_rb, real_names[args['coh_name']])
                if args['silent']:
                    job_args['-out:file:silent'] = '%s_%s_on_%s_%s_%i.out' % (args['coh_name'], args['doc_name'], t_coh,
                                                                              t_doc, rb)
                    job_args['-out:file:silent_struct_type'] = 'binary'
                script_vars = {
                    'doc_pdb=': '%s%s_B.pdb' % (path_seqs_rb, t_doc.lower()),
                    'coh_aln_file=': args['coh_aln_file'] + '_aln',
                    'doc_aln_file=': args['doc_aln_file'] + '_aln',
                    'coh_query_name=': args['coh_name'].lower()[:4],
                    'doc_query_name=': args['doc_name'].lower()[:4],
                    'coh_template_name=': t_coh + '.A',
                    'doc_template_name=': t_doc + '.B',
                    'coh_rb_pos=': str(cohs_dict[t_coh]['rb_coh']) + 'A',
                    'doc_rb_pos=': str(docs_dict[t_doc]['rb_doc'] + cohs_dict[t_coh]['coh_len']) + 'B',
                    'doc_start_res=': str(cohs_dict[t_coh]['coh_len'] + 1),
                    # 'hot_spot=': str(cohs_dict[t_coh]['rb_coh']+2) + 'A',
                }  # hot_spot is coh_rb_pos +2...
                make_jobs(args, job_args, script_vars)
    append_to_pending(os.getcwd() + '/')


def append_to_pending(folder: str):
    with open('/home/labs/fleishman/jonathaw/general_lists/pending_folders.txt', 'a') as fout:
        fout.write('%s\n' % folder)


def design_parts(args):
    from PSSMAnalyser import parse_pssm
    args['protocol'] = 'DoCohDesigner_3_13Oct.xml'
    shutil.copyfile(path_seqs_rb + 'all_docs.fasta', os.getcwd() + '/all_doc.fasta')
    coh_seq_1ohz = AASeq(
        'GVVVEIGKVTGSVGTTVEIPVYFRGVPSKGIANCDFVFRYDPNVLEIIGIDPGDIIVDPNPTKSFDTAIYPDRKIIVFLFAEDSGTGAYAITKDGVFAKIRATVKSSAPGYITFDEVGGFADNDLVEQKVSFIDGGVNVG')
    coh_name = args['coh_name']
    cohs_dict = {'rb_coh': RBposCoh(coh_seq_1ohz.get_seq(), coh_name), 'coh_len': len(coh_seq_1ohz)}
    coh_pssm = parse_pssm('/home/labs/fleishman/jonathaw/data/pssm/cohs/%s.pssm' % coh_name)

    # docs = [a.lower() for a in DOC_TEMPLATES if a not in ['2CCL', '3KCP']]
    docs_fasta = seq_funcs.read_multi_fastas(args['doc_aln_file'])

    doc_rb_select = {'1ohz': [10, 11, 17, 18, 19, 1, 20, 22, 25, 27, 2, 32, 3, 5, 6, 8, 9],
                     '2vn5': [10, 16, 17, 18, 19, 1, 21, 25, 26, 27, 28, 2, 30, 32, 4, 5, 6, 9],
                     '2y3n': [11, 12, 13, 18, 19, 1, 29, 31, 32, 9],
                     '3ul4': [10, 11, 12, 13, 17, 1, 20, 26, 28, 2, 31, 32, 5, 6, 9],
                     '4dh2': [11, 12, 14, 16, 18, 19, 1, 26, 29, 2, 30, 31, 32, 3, 5, 9],
                     '4fl4': [10, 11, 12, 16, 18, 19, 1, 21, 22, 25, 26, 27, 2, 30, 31, 32, 4, 5, 6, 9],
                     '4fl5': [10, 19, 20, 21, 22, 26, 32],
                     '4uyp': [11, 19, 20, 25, 26, 27, 28, 2, 30, 9],
                     '5new': [10, 11, 12, 14, 16, 17, 19, 1, 20, 22, 29, 2, 31, 32, 3, 6, 9]}
    docs = list(doc_rb_select.keys())
    res_type = {'"KR"': 'p', '"DE"': 'n', '"ACDEFGHIKLMNPQRSTVWY"': '-'}

    docs_dict = {k: {'rb_doc': None} for k in docs}
    for t_doc in docs:
        docs_dict[t_doc]['rb_doc'] = RBposDoc1st(docs_fasta[t_doc + '.B'].get_seq, t_doc.upper())


    # all_coh_poses = [33, 35, 37, 39, 63, 66, 68, 70, 73, 75, 77, 79, 81, 83, 85, 121, 125, 127, 118, 119, 115, 116, 32, 123,
    #                  82, 87]
    coh_switch_poses = [33, 35, 37, 63, 66, 70, 73, 75, 83, 127]
    new_combos = create_combos(res_type.keys(), coh_pssm, res_type)
    for c in new_combos:
        args['logger'].log(''.join(res_type[c[i]] for i in coh_switch_poses))
    args['logger'].log('combo num', len(new_combos))
    args['logger'].log('there are %i combos' % len(new_combos))

    # coh_no_design = design_to_no_design(1, cohs_dict['coh_len'], all_coh_poses)

    args['logger'].log('every do will have %i jobs' % (len(new_combos) * len(list(RB_NAMES.keys()))))
    args['logger'].log('there will be %i jobs created' % (len(docs) * len(list(RB_NAMES.keys())) * len(new_combos)))

    for t_doc in docs:
        if os.path.isdir('%s%s_on_%s' % (args['path'], args['coh_name'], t_doc)):
            continue
        os.chdir(args['path'])
        os.mkdir('%s%s_on_%s' % (args['path'], args['coh_name'], t_doc))
        os.chdir('%s%s_on_%s' % (args['path'], args['coh_name'], t_doc))

        # doc_desgin_pos = sorted([a for a in doc_design_pos[t_doc]['design'].keys()])
        # doc_desgin_pos_corrected = [a+cohs_dict['coh_len'] for a in doc_desgin_pos]

        # doc_no_design = [i for i in range(cohs_dict['coh_len']+1, cohs_dict['coh_len']+
        #                                   len(docs_fasta[t_doc+'.B'].get_seq())+1) if i not in doc_desgin_pos_corrected]

        for combo in new_combos:
            comb_str = ''.join(res_type[combo[i]] for i in coh_switch_poses)
            # if comb_str.count('n') < 3 or comb_str.count('p') < 3:
            #     continue
            os.chdir('%s%s_on_%s' % (args['path'], args['coh_name'], t_doc))
            os.mkdir('%s%s_on_%s/dir%s' % (args['path'], args['coh_name'], t_doc, comb_str))
            os.chdir('%s%s_on_%s/dir%s' % (args['path'], args['coh_name'], t_doc, comb_str))
            for rb in RB_NAMES.keys():
                # exclude RBIns that were found to be pointless... 15.10.15
                if rb not in doc_rb_select[t_doc]:
                    continue

                job_args = {
                    '-parser:protocol': '/home/labs/fleishman/jonathaw/protocols/%s' % args['protocol'],
                    '@/home/labs/fleishman/jonathaw/protocols/flags/%s' % args['flags_file']: None,
                    '-database': '/home/labs/fleishman/jonathaw/Rosetta/main/database/',
                    '-overwrite': None,
                    '-s': '%s%s_A.pdb' % (path_seqs_rb, coh_name.lower()),
                    '-out:prefix': 'dz_',
                    '-out:suffix': '_%s_%i_%s.pdb' % (t_doc, rb, comb_str),
                }
                script_vars = {
                    'RB_file=': '%s%s' % (rb_path, RB_NAMES[rb]),
                    'doc_pdb=': '%s%s_B.pdb' % (path_seqs_rb, t_doc.lower()),
                    'coh_rb_pos=': str(cohs_dict['rb_coh']) + 'A',
                    'doc_rb_pos=': str(docs_dict[t_doc]['rb_doc'] + cohs_dict['coh_len']) + 'B',
                    'doc_start_res=': str(cohs_dict['coh_len'] + 1),
                    # 'coh_no_design=': ','.join(str(i) for i in coh_no_design),
                    # 'doc_no_design=': ','.join(str(i) for i in doc_no_design),
                    'coh_pssm=': '/home/labs/fleishman/jonathaw/data/pssm/cohs/%s.pssm' % coh_name,
                    'doc_pssm=': '/home/labs/fleishman/jonathaw/data/pssm/docs/%s.pssm' % t_doc,
                }
                for k, v in combo.items():
                    script_vars['res_types_%i=' % k] = v
                make_jobs(args, job_args, script_vars)
                # sys.exit()


def create_combos(res_types, pssm, res_type):
    combinations = []
    for res33 in res_types:
        if not can_pos_by_pssm(res33, pssm[33]):
            continue
        for res35 in res_types:
            if not can_pos_by_pssm(res35, pssm[35]):
                continue
            for res37 in res_types:
                if not can_pos_by_pssm(res37, pssm[37]):
                    continue
                for res63 in res_types:
                    if not can_pos_by_pssm(res63, pssm[63]):
                        continue
                    for res66 in res_types:
                        if not can_pos_by_pssm(res66, pssm[66]):
                            continue
                        for res68 in res_types:
                            if not can_pos_by_pssm(res68, pssm[68]):
                                continue
                            for res70 in res_types:
                                if not can_pos_by_pssm(res70, pssm[70]):
                                    continue
                                for res73 in res_types:
                                    if not can_pos_by_pssm(res73, pssm[73]):
                                        continue
                                    for res75 in res_types:
                                        if not can_pos_by_pssm(res75, pssm[75]):
                                            continue
                                        for res83 in res_types:
                                            if not can_pos_by_pssm(res83, pssm[83]):
                                                continue
                                            for res121 in res_types:
                                                if not can_pos_by_pssm(res121, pssm[121]):
                                                    continue
                                                for res125 in res_types:
                                                    if not can_pos_by_pssm(res125, pssm[125]):
                                                        continue
                                                    for res127 in res_types:
                                                        if not can_pos_by_pssm(res127, pssm[127]):
                                                            continue
                                                        combo = {33: res33, 35: res35,
                                                                 37: res37,
                                                                 63: res63, 66: res66,
                                                                 68: res68, 70: res70,
                                                                 73: res73, 75: res75,
                                                                 83: res83,
                                                                 121: res121, 125: res125,
                                                                 127: res127}
                                                        s = ''.join(res_type[i] for i in combo.values())
                                                        if s.count('n') + s.count('p') >= 5:
                                                            combinations.append(combo)
    return combinations


def can_pos_by_pssm(reses, pssm_dict):
    return any(pssm_dict[res] >= 0 for res in reses if res != '"')


def can_comb_by_pssm(comb, pssm):
    for pos, aas in comb.items():
        if not any([pssm[pos][aa] >= 0 for aa in aas if aa != '"']):
            # print('didnt find!!!', pos, aas, ' '.join('%s %i' % (aa, pssm[pos][aa]) for aa in aas if aa != '"')
            return False
    return True


def design_to_no_design(start, end, design_poses):
    """
    :param start: start pos, the first residue for this chain (1 or len(coh)+1)
    :param end: the length of the chain
    :param design_poses: positions to be designed
    :return: postions to be not designed
    >>> start, end = 101, 110
    >>> design_poses = {3: 'A', 8: 'B'}
    >>> design_to_no_design(start, end, design_poses)
    [101, 102, 104, 105, 106, 107, 109, 110]
    """
    normed = [a + start - 1 for a in design_poses]
    result = []
    for i in range(start, end + 1):
        if i not in normed:
            result.append(i)
    return result


def make_jobs(args, job_args, script_vars=None):
    pwd = os.getcwd() + '/'
    num = ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(8))
    jobname = '%sjob.%s' % (pwd, num)
    outname = '%sout.%s' % (pwd, num)
    errname = '%serr.%s' % (pwd, num)
    cmdname = '%scommand' % pwd
    with open(jobname, 'w+') as job:
        job.write('#!/bin/bash\n')
        job.write('. /usr/share/lsf/conf/profile.lsf\n')
        job.write('cd ' + pwd + '\n')
        if args['queue'] == 'new-all.q':
            job.write(
                '/apps/RH6U4/blcr/0.8.5/bin/cr_run /home/labs/fleishman/jonathaw/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease ')
        else:
            job.write('/home/labs/fleishman/jonathaw/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease ')
        for k, v in job_args.items():
            if v is None:
                job.write('%s ' % k)
            else:
                job.write('%s %s ' % (k, v))
        if script_vars is not None:
            job.write('-parser:script_vars ')
            for k, v in script_vars.items():
                job.write('%s%s ' % (k, v))
        job.write('\n')
    with open(cmdname, 'a+') as cmd:
        if args['queue'] == 'new-all.q':
            cmd.write(str(
                'bsub -C 1024 -u /dev/null -N -u /dev/null -R rusage[mem=1024] -L /bin/bash -G fleishman-wx-grp-lsf -q ' +
                args['queue'] + ' -o ' + outname + ' -e ' + errname + ' /apps/RH6U4/blcr/0.8.5/bin/cr_run ' +
                jobname + '\n'))
        else:
            cmd.write(str('bsub -L /bin/bash -N -u /dev/null -G fleishman-wx-grp-lsf -q ' + args['queue'] + ' -o ' +
                          outname + ' -e ' + errname + ' ' + jobname + '\n'))
    subprocess.call(['chmod', '+x', jobname])


def add_align(args, coh_doc) -> None:
    aln_file = args['coh_aln_file'] if coh_doc == 'coh' else args['doc_aln_file']
    with open(aln_file, 'a+') as fout:
        fout.write('>%s\n%s\n' % (args['coh_name'] if coh_doc == 'coh' else args['doc_name'],
                                  args['coh_seq'] if coh_doc == 'coh' else args['doc_seq']))
    if coh_doc == 'coh':
        subprocess.call(['muscle', '-in', aln_file, '-out', aln_file + '_aln'])  # , '2>&1'])
    elif coh_doc == 'doc':
        subprocess.call(['muscle', '-in', aln_file, '-out', aln_file + '_aln'])  # , '2>&1'])
    with open(aln_file + '_aln', 'a') as fout:
        fout.write('>mock\nAAAAAAAAAA\n')


def read_multi_fastas(file):
    """
    :param file: fasta with gaps MSA
    :return: {name: {name:, seq:}} dictionary. sequences whose name are in the TEMPLATE lists will be named with the
    short name only
    """
    with open(file, 'r') as f:
        cont = f.read().split('>')
    result = {}
    for entry in cont:
        split_entry = entry.split('\n')
        if len(split_entry) < 2:
            continue
        name = '_'.join(split_entry[0].rstrip().split())
        seq = ''.join(a.rstrip() for a in split_entry[1:])
        if name[:4].upper() not in DOC_TEMPLATES and name[:4].upper() not in COH_TEMPLATES:
            result[name.lower()] = {'name': name.lower(), 'seq': seq}
        else:
            result[name[:4].lower()] = {'name': name[:4].lower(), 'seq': seq}
    return result


def docs_pass(args, doc_fasta):
    """
    :param args: user arguments
    :param cohs_fasta: dockerins MSA in dict fasta
    :return: list of doc names that are good for threading
    """
    passed = []
    query_seq = doc_fasta[args['doc_name'].lower()]['seq']
    for clstr_rep, temp_strands in DOC_CLUSTER.items():
        clst_seq = doc_fasta[clstr_rep]['seq']
        clst_seg = temp_strands['strands']
        clst_tmp = temp_strands['temp']
        cruc_segs = poses_with_gaps_of_seg(clst_seq, clst_seg)
        passes, seg_fail = same_gap_pattern(query_seq, clst_seq, cruc_segs, args['passing_problem'])
        if passes:
            passed.extend(clst_tmp)
        else:
            args['logger'].log('%s DOC is failed by seg %i %i' % (clstr_rep, seg_fail[0], seg_fail[1]))
            args['logger'].log('cluster seq %s' % clst_seq[seg_fail[0]:seg_fail[1] + 1])
            args['logger'].log('query   seq %s\n' % query_seq[seg_fail[0]:seg_fail[1] + 1])
    return passed


def cohs_pass(args, cohs_fasta):
    """
    :param args: user arguments
    :param cohs_fasta: cohesin MSA in dict fasta
    :return: list of coh names that are good for threading
    """
    passed = []
    try:
        query_seq = cohs_fasta[args['coh_name'].split('.')[0].lower()]['seq']
    except:
        query_seq = cohs_fasta[args['coh_name'].lower()]['seq']
    for clst_rep, temps_strands in COH_CLUSTER.items():
        clst_seq = cohs_fasta[clst_rep]['seq']
        clst_seg = [a[2:-2] for a in temps_strands['strands']]
        clst_tmp = temps_strands['temp']
        cruc_segs = poses_with_gaps_of_seg(clst_seq, clst_seg)
        passes, seg_fail = same_gap_pattern(query_seq, clst_seq, cruc_segs, args['passing_problem'])
        if passes:
            passed.extend(clst_tmp)
        else:
            args['logger'].log('%s COH is failed by seg %i %i' % (clst_rep, seg_fail[0], seg_fail[1]))
            args['logger'].log('cluster seq %s' % clst_seq[seg_fail[0]:seg_fail[1] + 1])
            args['logger'].log('query   seq %s\n' % query_seq[seg_fail[0]:seg_fail[1] + 1])
    return passed


def same_gap_pattern(q_seq, t_seq, segs, passing_problem=False):
    """
    :param q_seq: query seq (with gaps)
    :param t_seq: target seq (with gaps)
    :param segs: list of segments ot test
    :return: True iff query and target have the same gap pattern
    >>> q_seq = 'AAAAAAAAAAA-AAAAAAA'
    >>> t_seq = 'AAAAACCCAAA-BBBBBBB'
    >>> segs = [[0, 3], [9, 13]]
    >>> same_gap_pattern(q_seq, t_seq, segs)
    (True, None)
    >>> t_seq = 'AAAAACCCAAA--BBBBBB'
    >>> same_gap_pattern(q_seq, t_seq, segs)
    (False, [9, 13])
    >>> t_seq = 'AAAAACCCAAA--BBBBB-'
    >>> segs = [[16, 18]]
    >>> same_gap_pattern(q_seq, t_seq, segs)
    (False, [16, 18])
    >>> q_seq = 'AAAAAAAAAAA-AAAAAA-'
    >>> same_gap_pattern(q_seq, t_seq, segs)
    (True, None)
    """
    if not passing_problem:
        for seg in segs:
            for i in range(seg[0], seg[1] + 1):
                if (t_seq[i] == '-' and q_seq[i] != '-') or (t_seq[i] != '-' and q_seq[i] == '-'):
                    return False, seg
        return True, None
    else:
        for seg in segs:
            counter = 0
            for i in range(seg[0], seg[1] + 1):
                if (t_seq[i] == '-' and q_seq[i] != '-') or (t_seq[i] != '-' and q_seq[i] == '-'):
                    counter += 1
            if counter > 1:
                return False, seg
        return True, None


def target_coverage(t_seq, q_seq):
    """
    :param t_seq: target sequence from alignment
    :param q_seq: query sequence from alignment
    :return: coverage calculated by the number of positions in target that have an aligned AA (non gap) in query,
    divided by the length of the target
    >>> t_seq = '123456-789-112'
    >>> q_seq = '-236---5678---'
    >>> target_coverage(t_seq, q_seq)
    0.5
    """
    counter = 0.
    for i in range(len(t_seq)):
        if t_seq[i] != '-' and q_seq[i] != '-':
            counter += 1.
    return counter / len(t_seq.replace('-', ''))


def poses_with_gaps_of_seg(seq_gaps, segs):
    """
	:param seq_gaps: a sequence with gaps
	:param seg: a list of AA segments
	:return: a list of the positions (i, j) of where segs are within seq_gaps
	>>> seq_gaps = 'ABC-DEF-GHIJKLMNOP-QRST'
	>>> seg = ['EFG', 'OPQ']
	>>> poses_with_gaps_of_seg(seq_gaps, seg)
	[[5, 8], [16, 19]]
	"""
    result = []
    for seg in segs:
        string = '-*'.join(a for a in seg)  # +'.*'
        patt = re.compile(string)
        result.append([[m.start(0), m.end(0) - 1] for m in re.finditer(patt, seq_gaps)][0])
    return result


def RBposCoh(coh_target_seq, coh_model):
    NAME = coh_model.upper()

    if NAME == '1OHZ' or NAME == '1OHA':
        pat = re.compile("V-*F-*L-*F-*A.*")
    elif NAME == '2B59':
        pat = re.compile("N-*F-*A-*L-*A.*")
    elif NAME == '2CCL':
        pat = re.compile("V-*F-*L-*F-*A.*")
    elif NAME == '2OZN':
        pat = re.compile("R-*V-*L-*V-*S.*")
    elif NAME == '2VN5':
        pat = re.compile("S-*F-*L-*F-*L.*")
    elif NAME == '2Y3N':
        pat = re.compile("N-*F-*G-*R-*L.*")
    elif NAME == '3UL4':
        pat = re.compile("A-*V-*L-*Y-*L.*")
    elif NAME == '3KCP':
        pat = re.compile("N-*F-*A-*L-*A.*")
    elif NAME == '4DH2':
        pat = re.compile("V-*F-*L-*F-*D.*")
    elif NAME == '4IU2':
        pat = re.compile("F-*V-*A-*S-*G.*")
    elif NAME == '4UYP':
        pat = re.compile("S-*M-*T-*F-*E.*")
    elif NAME == '4UYQ':
        pat = re.compile("S-*M-*T-*F-*E.*")
    elif NAME == '4UMS':
        pat = re.compile("S-*F-*L-*F-*N.*")

    elif NAME == '1ANU':
        pat = re.compile("V-*F-*L-*F-*A.*")
    elif NAME == '1AOH':
        pat = re.compile("V-*F-*L-*F-*A.*")
    elif NAME == '1G1K':
        pat = re.compile("S-*F-*L-*F-*L.*")
    elif NAME == '1QZN':
        pat = re.compile("N-*F-*S-*K-*A.*")
    elif NAME == '1TYJ':
        pat = re.compile("N-*F-*G-*R-*L.*")
    elif NAME == '2BM3':
        pat = re.compile("N-*F-*A-*L-*A.*")
    elif NAME == '2JH2':
        pat = re.compile("R-*V-*L-*V-*S.*")
    elif NAME == '2VN6':
        pat = re.compile("S-*F-*L-*F-*L.*")
    elif NAME == '2VO8':
        pat = re.compile("R-*I-*L-*V-*A.*")
    elif NAME == '2W1N':
        pat = re.compile("R-*V-*L-*V-*S.*")
    elif NAME == '2W5F':
        pat = re.compile("A-*K-*Y-*K-*A.*")
    elif NAME == '2XBT':
        pat = re.compile("E-*I-*K-*L-*A.*")
    elif NAME == '2XDH':
        pat = re.compile("K-*V-*G-*I-*A.*")
    elif NAME == '2ZF9':
        pat = re.compile("F-*F-*T-*A-*T.*")
    elif NAME == '3BWZ':
        pat = re.compile("N-*F-*S-*K-*A.*")
    elif NAME == '3FNK':
        pat = re.compile("L-*S-*R-*S-*Y.*")
    elif NAME == '3L8Q':
        pat = re.compile("N-*F-*S-*K-*A.*")
    elif NAME == '4FL4':
        pat = re.compile("V-*F-*L-*F-*A.*")
    elif NAME == '4JO5':
        pat = re.compile("E-*I-*S-*F-*T.*")
    elif NAME == '4FL5':
        pat = re.compile("N-*F-*A-*L-*A.*")
    elif NAME == '5NEW':
        pat = re.compile("S-*V-*I-*W-*S.*")
    subseq = coh_target_seq[0:pat.search(coh_target_seq).start() + 1]
    # cleaned = subseq.translate(None, '-')
    cleaned = subseq.replace('-', '')
    return len(cleaned)


def RBposDoc1st(doc_target_seq, doc_model):
    NAME = doc_model.upper()
    if NAME == '1OHZ' or NAME == '1OHA':
        pat = re.compile("S-*T-*D-*L-*T.*")
    # elif NAME == '2B59':
    #     pat = re.compile("L-*D-*V-*A-*E.*")
    elif NAME == '2B59':
        pat = re.compile("L-*L-*D-*V-*A.*")
    elif NAME == '2CCL':
        pat = re.compile("S-*T-*D-*L-*T.*")
    elif NAME == '2OZN':
        pat = re.compile("I-*G-*D-*L-*A.*")
    elif NAME == '2VN5':
        pat = re.compile("A-*L-*D-*F-*A.*")
    elif NAME == '2Y3N':
        pat = re.compile("M-*A-*D-*V-*M.*")
    elif NAME == '3UL4':
        pat = re.compile("D-*E-*D-*Y-*I.*")
    elif NAME == '3KCP':
        pat = re.compile("L-*L-*D-*V-*A.*")
    elif NAME == '4DH2':
        pat = re.compile("I-*S-*D-*Y-*V.*")
    elif NAME == '4IU2':
        pat = re.compile("G-*R-*D-*A-*T.*")
    elif NAME == '4UYP':
        pat = re.compile("S-*I-*D-*A-*V.*")
    elif NAME == '4UYQ':
        pat = re.compile("I-*N-*D-*A-*V.*")
    elif NAME == '2VN6':
        pat = re.compile("S-*T-*D-*F-*A.*")
    elif NAME == '4FL4':
        pat = re.compile("S-*T-*D-*L-*T.*")
    elif NAME == '4FL5':
        pat = re.compile("L-*L-*D-*V-*A.*")
    elif NAME == '5NEW':
        pat = re.compile("S-*D-*D-*L-*T.*")
    subseq = doc_target_seq[0:pat.search(doc_target_seq).start() + 1]
    cleaned = subseq.replace('-', '')
    return len(cleaned)


def seq_iden(s1, s2):
    """
    :param s1: sequence 1
    :param s2: sequence 1
    :return: sequence identity, calculated as the number of identical residues, devided by the length of s1
    >>> s1 = '123456789-'
    >>> s2 = '123456-890'
    >>> seq_iden(s1, s2)
    0.8888888888888888
    """
    result = 0
    for i in range(len(s1)):
        if s1[i] == s2[i] and s1[i] != '-':
            result += 1
    return float(result) / float(len(s1.replace('-', '')))


def jack_simple(args):
    """
    :param args: run arguments 
    :return: crates jons and command file to predict a structure with known structure using a simple docking procedure
    """
    args['flags_file'] = 'docoh_modeling.flags'
    args['doc_name'] = args['coh_name']
    args['protocol'] = 'DoCohModeller_simple_docking_19Nov.xml'
    print("copying data files")
    shutil.copyfile(path_seqs_rb + 'all_cohs.fasta', os.getcwd() + '/all_coh.fasta')
    shutil.copyfile(path_seqs_rb + 'all_docs.fasta', os.getcwd() + '/all_doc.fasta')

    print("adding query seqs and aligning:")
    add_align(args, 'coh')
    add_align(args, 'doc')

    coh_fasta = read_multi_fastas(args['coh_aln_file'] + '_aln')
    doc_fasta = read_multi_fastas(args['doc_aln_file'] + '_aln')

    coh_iden = {}
    for coh in coh_fasta.values():
        if coh['seq'] == 'AAAAAAAAAA'or coh['name'].upper() not in DOC_TEMPLATES:
            continue
        iden = seq_iden(coh_fasta[args['coh_name']]['seq'], coh['seq'])
        if iden < 0.9:
            coh_iden[iden] = coh['name']

    doc_iden = {}
    for doc in doc_fasta.values():
        if doc['seq'] == 'AAAAAAAAAA':
            continue
        iden = seq_iden(doc_fasta[args['doc_name']]['seq'], doc['seq'])
        if iden < 0.9 and not (args['coh_name'] == '2v5a' and doc['name'] == '2vn6'):
            doc_iden[iden] = doc['name']
    max_coh = coh_iden[max(coh_iden.keys())]
    max_doc = doc_iden[max(doc_iden.keys())]
    print(max_coh, max(coh_iden.keys()))
    print(max_doc, max(doc_iden.keys()))

    for chosen in {max_coh, max_doc}:
        os.mkdir('on_%s' % chosen)
        os.chdir('on_%s' % chosen)
        shutil.copy('../'+args['coh_aln_file'] + '_aln', '.')
        shutil.copy('../'+args['doc_aln_file'] + '_aln', '.')

        coh_len = len(coh_fasta[chosen]['seq'].replace('-', ''))
        for itr in range(5000):
                job_args = {
                    '-parser:protocol': '/home/labs/fleishman/jonathaw/protocols/%s' % args['protocol'],
                    '@/home/labs/fleishman/jonathaw/protocols/flags/%s' % args['flags_file']: None,
                    '-database': '/home/labs/fleishman/jonathaw/Rosetta/main/database/',
                    '-overwrite': None,
                    '-s': '%sno_calcium/%s_AB.pdb' % (path_seqs_rb, chosen),
                    '-out:prefix': '%s_on_' % args['coh_name'],
                    '-out:suffix': '_%i.pdb' % itr,
                    '-native': '%sno_calcium/%s_AB.pdb' % (path_seqs_rb, real_names[args['coh_name']]),
                    '-spin': 1
                }
                script_vars = {
                    'coh_aln_file=': args['coh_aln_file'] + '_aln',
                    'doc_aln_file=': args['doc_aln_file'] + '_aln',
                    'coh_query_name=': args['coh_name'].lower()[:4],
                    'doc_query_name=': args['doc_name'].lower()[:4],
                    'coh_template_name=': chosen + '.A',
                    'doc_template_name=': chosen + '.B',
                    'doc_start_res=': coh_len + 1,
                }
                make_jobs(args, job_args, script_vars)
        append_to_pending(os.getcwd() + '/')
        os.chdir('../')


def jack_no_docking(args):
    """
    :param args: run arguments
    :return: crates jons and command file to predict a structure with known structure using a simple docking procedure
    """
    args['flags_file'] = 'docoh_modeling.flags'
    args['doc_name'] = args['coh_name']
    args['protocol'] = 'DoCohModeller_NOdocking_22Dec.xml'
    print("copying data files")
    shutil.copyfile(path_seqs_rb + 'all_cohs.fasta', os.getcwd() + '/all_coh.fasta')
    shutil.copyfile(path_seqs_rb + 'all_docs.fasta', os.getcwd() + '/all_doc.fasta')

    print("adding query seqs and aligning:")
    add_align(args, 'coh')
    add_align(args, 'doc')

    coh_fasta = read_multi_fastas(args['coh_aln_file'] + '_aln')
    doc_fasta = read_multi_fastas(args['doc_aln_file'] + '_aln')

    coh_iden = {}
    for coh in coh_fasta.values():
        if coh['seq'] == 'AAAAAAAAAA'or coh['name'].upper() not in DOC_TEMPLATES:
            continue
        iden = seq_iden(coh_fasta[args['coh_name']]['seq'], coh['seq'])
        print('iden for coh %s %f' % (coh['name'], iden))
        if iden < 0.9:
            coh_iden[iden] = coh['name']

    doc_iden = {}
    for doc in doc_fasta.values():
        if doc['seq'] == 'AAAAAAAAAA':
            continue
        iden = seq_iden(doc_fasta[args['doc_name']]['seq'], doc['seq'])
        if iden < 0.9 and not (args['coh_name'] == '2v5a' and doc['name'] == '2vn6'):
            doc_iden[iden] = doc['name']
    max_coh = coh_iden[max(coh_iden.keys())]
    max_doc = doc_iden[max(doc_iden.keys())]
    print(max_coh, max(coh_iden.keys()))
    print(max_doc, max(doc_iden.keys()))

    for chosen in {max_coh, max_doc}:
        os.mkdir('on_%s' % chosen)
        os.chdir('on_%s' % chosen)
        shutil.copy('../'+args['coh_aln_file'] + '_aln', '.')
        shutil.copy('../'+args['doc_aln_file'] + '_aln', '.')

        coh_len = len(coh_fasta[chosen]['seq'].replace('-', ''))
        for itr in range(1):
                job_args = {
                    '-parser:protocol': '/home/labs/fleishman/jonathaw/protocols/%s' % args['protocol'],
                    '@/home/labs/fleishman/jonathaw/protocols/flags/%s' % args['flags_file']: None,
                    '-database': '/home/labs/fleishman/jonathaw/Rosetta/main/database/',
                    '-overwrite': None,
                    '-s': '%sno_calcium/%s_AB.pdb' % (path_seqs_rb, chosen),
                    '-out:prefix': '%s_on_' % args['coh_name'],
                    '-out:suffix': '_%i.pdb' % itr,
                    '-native': '%sno_calcium/%s_AB.pdb' % (path_seqs_rb, real_names[args['coh_name']]),
                    '-spin': 1
                }
                script_vars = {
                    'coh_aln_file=': args['coh_aln_file'] + '_aln',
                    'doc_aln_file=': args['doc_aln_file'] + '_aln',
                    'coh_query_name=': args['coh_name'].lower()[:4],
                    'doc_query_name=': args['doc_name'].lower()[:4],
                    'coh_template_name=': chosen + '.A',
                    'doc_template_name=': chosen + '.B',
                    'doc_start_res=': coh_len + 1,
                }
                make_jobs(args, job_args, script_vars)
        append_to_pending(os.getcwd() + '/')
        os.chdir('../')

if __name__ == '__main__':
    import argparse
    from Logger import Logger as lg

    parser = argparse.ArgumentParser()
    parser.add_argument('-coh_name', type=str)
    parser.add_argument('-doc_name', type=str)
    parser.add_argument('-coh_seq', type=str)
    parser.add_argument('-doc_seq', type=str)
    parser.add_argument('-coh_aln_file', default="./all_coh.fasta")
    parser.add_argument('-doc_aln_file', default="./all_doc.fasta")
    parser.add_argument('-coh_seqs_file')
    parser.add_argument('-doc_seqs_file')
    parser.add_argument('-mode', type=str, default='predict_pair')
    parser.add_argument('-job_num', type=int)
    parser.add_argument('-protocol', type=str)
    parser.add_argument('-flags_file', default='docoh_modeling.flags')
    parser.add_argument('-queue', default='fleishman')
    parser.add_argument('-log_file', default='jobber.log')
    parser.add_argument('-make_log', default=True, type=bool)
    parser.add_argument('-path', default=os.getcwd() + '/')
    parser.add_argument('-silent', type=bool, default=False)
    parser.add_argument('-passing_problem', type=bool, default=False, help='if a jack knifing procedure has nothing to '
                                                                           'be modeled on, i reduce the contraints to 1'
                                                                           ' gap in segment')
    parser.add_argument('-specific_RBs', type=bool, default=False, help="in predicting pairs, whehter to use only RBs "
                                                                        "that passed a test. to use for diagonal "
                                                                        "prediction")
    parser.add_argument('-pdb', type=str)
    args = vars(parser.parse_args())

    args['logger'] = lg(log_file=args['log_file'])

    if args['mode'] == 'predict_pair':
        # args['protocol'] = 'DoCohModeller.xml'
        args['protocol'] = 'DoCohModeller_normsd_9Nov_hbonds.xml'
        predict_pair(args)

    elif args['mode'] == 'jack':
        # args['protocol'] = 'DoCohModeller_rmsd.xml'
        predict_pair(args)

    elif args['mode'] == 'post':
        # args['protocol'] = 'DoCohModeller_normsd_8Oct.xml'
        args['protocol'] = 'DoCohModeller_normsd_9Nov_hbonds.xml'
        cohs = read_multi_fastas('/home/labs/fleishman/jonathaw/data/seqs_all/cohesins_from_rachel.fasta')
        print(cohs.keys())
        args['coh_seq'] = cohs[args['coh_name']]['seq']
        docs = read_multi_fastas('/home/labs/fleishman/jonathaw/data/seqs_all/dockerins_from_rachel.fasta')
        args['doc_seq'] = docs[args['doc_name']]['seq']
        predict_pair(args)

    elif args['mode'] == 'all_coh_post':
        assert os.getcwd().split('/')[-1] == args['coh_name'], 'not in the right folder'
        args['protocol'] = 'DoCohModeller_normsd_9Nov_hbonds.xml'
        cohs = read_multi_fastas('/home/labs/fleishman/jonathaw/data/seqs_all/cohesins_from_rachel.fasta')
        args['coh_seq'] = cohs[args['coh_name']]['seq']
        # args['protocol'] = 'DoCohModeller_normsd_8Oct.xml'
        docs = read_multi_fastas('/home/labs/fleishman/jonathaw/data/seqs_all/dockerins_from_rachel.fasta')
        for doc in docs.keys():
            if doc.upper() not in DOC_TEMPLATES and doc not in real_names.keys() and doc.upper() not in COH_TEMPLATES:
                # print('found doc', doc)
                os.chdir(args['path'])
                os.mkdir('%s%s_on_%s' % (args['path'], args['coh_name'], doc))
                os.chdir('%s%s_on_%s' % (args['path'], args['coh_name'], doc))
                args['logger'] = lg('jobber.log')
                args['doc_name'] = doc
                args['doc_seq'] = docs[doc]['seq']

                predict_pair(args)

    elif args['mode'] == 'design_parts':
        design_parts(args)

    elif args['mode'] == 'survey':
        cohs = read_multi_fastas(
            '/home/labs/fleishman/jonathaw/data/pdbs_final/chosen_chains_rearranged/all_cohs.fasta_aln')
        docs = read_multi_fastas(
            '/home/labs/fleishman/jonathaw/data/pdbs_final/chosen_chains_rearranged/all_docs.fasta_aln')
        for qk, qv in cohs.items():
            print('query is', qk)
            args['coh_name'] = qk
            coh_passed = cohs_pass(args, cohs)
            print('by strands these pass:', coh_passed)
            for tk, tv in cohs.items():
                if qk != tk:
                    print('for targes %s the covergae is %f' % (tk, target_coverage(qv['seq'], tv['seq'])))

    elif args['mode'] == 'predict_pair_from_pdb':
        pdb = MyPDB.parse_PDB(args['pdb'], args['coh_name'])
        seqs = MyPDB.extract_seq(pdb)
        args['coh_seq'] = seqs['A'].get_seq
        args['doc_seq'] = seqs['B'].get_seq
        args['coh_name'] = args['pdb'].split('.pdb')[0]
        args['doc_name'] = args['pdb'].split('.pdb')[0]
        # os.mkdir(args['path']+args['coh_name']+'/')
        # os.chdir(args['path']+args['coh_name']+'/')

        args['protocol'] = 'DoCohModeller_normsd_9Nov_hbonds.xml'
        predict_pair(args)
        # os.chdir(args['path']+'/')

    elif args['mode'] == 'RBIn_SASA_test':
        RBIn_SASA_test(args)

    elif args['mode'] == 'predict_seqs':
        args['protocol'] = 'DoCohModeller_normsd_9Nov_hbonds.xml'
        cohs_seqs = read_multi_fastas(args['coh_seqs_file'])
        docs_seqs = read_multi_fastas(args['doc_seqs_file'])
        args['coh_seq'] = cohs_seqs[args['coh_name'].lower()]['seq'] # ' + '.pdb.gz.a']['seq']
        args['doc_seq'] = docs_seqs[args['doc_name'].lower()]['seq'] #   + '.pdb.gz.b']['seq']
        print(args['coh_seq'])
        print(args['doc_seq'])

        predict_pair(args)

    elif args['mode'] == 'minidiagonal':
        minidiagonal(args)

    elif args['mode'] == 'jack_simple':
        jack_simple(args)

    elif args['mode'] == 'clq':
        args['protocol'] = 'DoCohModeller_normsd_9Nov_hbonds.xml'

        args['coh_seqs_file'] = '/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/reclique_18Nov/cliques_prediction/all_j_st_cohs.fasta'
        args['doc_seqs_file'] = '/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/reclique_18Nov/cliques_prediction/all_j_st_docs.fasta'

        path_here = os.getcwd()
        my_folder = path_here.split('/')[-1]
        args['coh_name'] = my_folder.split('_on_')[0]
        args['doc_name'] = my_folder.split('_on_')[1]
        cohs_seqs = read_multi_fastas(args['coh_seqs_file'])
        docs_seqs = read_multi_fastas(args['doc_seqs_file'])
        args['coh_seq'] = cohs_seqs[args['coh_name'].lower() + '_st.a']['seq']
        args['doc_seq'] = docs_seqs[args['doc_name'].lower() + '_st.b']['seq']

        print(args['coh_name'], args['doc_name'])

        predict_pair(args)

    elif args['mode'] == 'full_clq':
        from choose_cliques import all_i_know
        args['protocol'] = 'DoCohModeller_normsd_9Nov_hbonds.xml'
        already_know = all_i_know()
        args['coh_seqs_file'] = '/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/reclique_18Nov/cliques_prediction/all_j_st_cohs.fasta'
        args['doc_seqs_file'] = '/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/reclique_18Nov/cliques_prediction/all_j_st_docs.fasta'
        cohs_seqs = read_multi_fastas(args['coh_seqs_file'])
        docs_seqs = read_multi_fastas(args['doc_seqs_file'])

        path_here = os.getcwd()
        my_folder = path_here.split('/')[-1]

        clq_names = my_folder.split('_')
        for n1 in clq_names:
            for n2 in clq_names:
                if [n1, n2] not in already_know:
                    os.mkdir('%s_on_%s' % (n1, n2))
                    os.chdir('%s_on_%s' % (n1, n2))
                    args['coh_name'] = n1
                    args['doc_name'] = n2
                    args['coh_seq'] = cohs_seqs[args['coh_name'].lower() + '_st.a']['seq']
                    args['doc_seq'] = docs_seqs[args['doc_name'].lower() + '_st.b']['seq']

                    print("preparing %s Vs. %s" % (args['coh_name'], args['doc_name']))

                    predict_pair(args)
                    os.chdir('../')
                else:
                    print('skipping %s %s' % (n1, n2))

    elif args['mode'] == 'jack_no_docking':
        jack_no_docking(args)

    else:
        print('no mode was chosen')
    args['logger'].close()
