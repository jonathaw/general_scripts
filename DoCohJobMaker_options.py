#!/usr/bin/env python
import subprocess
import re
import sys
import os
import argparse
import shutil

coh_seq_file = '/home/labs/fleishman/jonathaw/data/seqs_all/cohesins_from_rachel.fasta'
doc_seq_file = '/home/labs/fleishman/jonathaw/data/seqs_all/dockerins_from_rachel.fasta'
path_seqs_rb = '/home/labs/fleishman/jonathaw/data/crys_seqs_and_RB_for_modeling/'
coh_pssm_path = '/home/labs/fleishman/jonathaw/pssm/cohs/'
doc_pssm_path = '/home/labs/fleishman/jonathaw/pssm/docs/'
COH_TEMPLATES = ['1OHZ', '2B59', '2CCL', '2OZN', '2VN5', '2VN6', '2Y3N', '3KCP', '3UL4', '4DH2', '4IU2', '4UYP', '4UYQ',
                 '4FL4', '1ANU', '1AOH', '1G1K', '1QZN', '1TYJ', '2BM3', '2JH2', '2VO8', '2W1N', '2W5F', '2XBT', '2XDH',
                 '2ZF9', '3BWZ', '3FNK', '3L8Q', '4JO5', '4UMS', '4FL5', '5NEW']
DOC_TEMPLATES = ['1OHZ', '2B59', '2CCL', '2OZN', '2VN5', '2VN6', '2Y3N', '3KCP', '3UL4', '4DH2', '4IU2', '4UYP', '4UYQ',
                 '4FL4', '4FL5', '5NEW']
is_first_line_in_cmd = True
names_Ad = {'1ohz': '1oha', '2b59': '2b5a', '2ccl':'2cca', '2ozn':'2oza', '2vn5': '2v5a', '2vn6': '2v6a',
            '2y3n': '2y3a', '3kcp': '3kca', '3ul4': '3ula', '4dh2': '4dha', '4iu2': '4iua', '4uyp': '4upa',
             '4uyq': '4uqa', '4fl4': '4f4a', '1anu': '1ana', '1aoh': '1aoa', '1g1k': '1g1a', '1qzn': '1qza',
             '1tyj': '1tya', '2bm3': '2bma', '2jh2': '2jha', '2vo8': '2voa', '2w1n': '2w1a', '2w5f': '2w5a',
             '2xbt': '2xba', '2xdh': '2xda', '2zf9': '2zfa', '3bwz': '3bwa', '3fnk': '3fna', '3l8q': '3l8a',
             '4jo5': '4joa', '4ums': '4uma', '4fl5': '4f5a', '5new': '5nea'}


def GetSequence(name, seq_file):
    seq = ''
    with open(seq_file, 'r') as file_temp:
        for line in file_temp:
            if line[1:1 + len(name)].upper() == name.upper():
                for line_1 in file_temp:
                    if line_1[0] == '>' or len(line_1.rstrip()) == 0:
                        break
                    seq += line_1.rstrip()
                    # break
        if len(seq) == 0:
            print 'could not find', name, 'breaking', seq, 'in file', seq_file
            sys.exit(0)
        return seq.split()[0]


def RBposCoh(aln_file, coh_model):
    IN = open(aln_file, 'r')
    NAME = coh_model.upper()

    seq = ''

    for line in IN:
        if line[1:5] == NAME:
            for i_line in IN:
                if i_line[0] != '>':
                    seq += i_line.rstrip()
                else:
                    break

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
        pat = re.compile("N-*L-*S-*R-*S.*")
    elif NAME == '4FL4':
        pat = re.compile("V-*F-*L-*F-*A.*")
    elif NAME == '4JO5':
        pat = re.compile("E-*I-*S-*F-*T.*")
    elif NAME == '4FL5':
        pat = re.compile("N-*F-*A-*L-*A.*")
    elif NAME == '5NEW':
        pat = re.compile("S-*V-*I-*W-*S.*")
    subseq = seq[0:pat.search(seq).start() + 1]
    cleaned = subseq.translate(None, '-')
    IN.close()
    return len(cleaned)


def RBposDoc1st(aln_file, doc_model):
    IN = open(aln_file, 'r')
    NAME = doc_model.upper()
    seq = ''
    for line in IN:
        if line[1:5] == NAME:
            for i_line in IN:
                if i_line[0] != '>':
                    seq += i_line.rstrip()
                else:
                    break

    if NAME == '1OHZ' or NAME == '1OHA':
        pat = re.compile("S-*T-*D-*L-*T.*")
    elif NAME == '2B59':
        pat = re.compile("L-*D-*V-*A-*E.*")
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
    subseq = seq[0:pat.search(seq).start() + 1]
    cleaned = subseq.translate(None, '-')
    IN.close()
    return len(cleaned)


def RBposDoc2nd(aln_file, doc_model):
    IN = open(aln_file, 'r')
    NAME = doc_model.upper()
    seq = ''

    for line in IN:
        if line[1:5] == doc_model:
            for i_line in IN:
                if i_line[0] != '>':
                    seq += i_line.rstrip()
                else:
                    break

    if NAME == '1OHZ':
        pat = re.compile("S-*T-*D-*V-*L.*")
    elif NAME == '2B59':
        pat = re.compile("M-*Q-*D-*I-*M.*")
    elif NAME == '2CCL':
        pat = re.compile("A-*A-*D-*V-*L.*")
    elif NAME == '2OZN':
        pat = re.compile("E-*Y-*E-*I-*S.*")
    elif NAME == '2VN5':
        pat = re.compile("S-*T-*D-*L-*A.*")
    elif NAME == '2Y3N':
        pat = re.compile("S-*D-*D-*A-*I.*")
    elif NAME == '3UL4':
        pat = re.compile("S-*T-*D-*C-*L.*")
    elif NAME == '3KCP':
        pat = re.compile("M-*Q-*D-*I-*M.*")
    elif NAME == '4DH2':
        pat = re.compile("D-*I-*D-*C-*N.*")
    elif NAME == '4IU2':
        pat = re.compile("G-*R-*D-*A-*S.*")
    elif NAME == '4UYP':
        pat = re.compile("I-*n-*d-*A-*V.*")
    elif NAME == '4UYQ':
        pat = re.compile("S-*I-*D-*A-*V.*")
    elif NAME == '2VN6':
        pat = re.compile("A-*F-*D-*L-*A.*")
    elif NAME == '4FL4':
        pat = re.compile("S-*S-*D-*V-*T.*")
    elif NAME == '4FL5':
        pat = re.compile("M-*Q-*D-*I-*M.*")
    elif NAME == '5NEW':
        pat = re.compile("V-*F-*D-*L-*I.*")
    subseq = seq[0:pat.search(seq).start() + 1]
    cleaned = subseq.translate(None, '-')
    IN.close()
    print len(cleaned)


def LengthSeqInAln(aln_file, name):
    IN = open(aln_file, 'r')
    seq = ''
    for line in IN:
        if line[1:1 + len(name)].upper() == name.upper():
            for i_line in IN:
                if i_line[0] != '>':
                    seq += i_line.rstrip()
                else:
                    break
    IN.close()
    cleaned = seq.translate(None, '-')
    return len(cleaned)


def AddAndAlign(name, seq, aln_file, coh_or_doc):
    print 'aligning', aln_file
    OUT = open(aln_file, 'a+')
    OUT.write('>' + name.upper() + '\n' + seq + '\n')
    OUT.close()
    with open(os.getcwd() + '/t_coffee.out', 'w') as stdout_file:
        with open(os.getcwd() + '/t_coffee.err', 'w') as stderr_file:
            if coh_or_doc == 'coh':
                subprocess.call(['t_coffee', aln_file, '-output', 'fasta_aln', '-template_file',
                                 '~/data/t_coffee_data/t_coffee_coh.template', '-method', 'TMalign_pair',
                                 '-mode', '3dcoffee', '-email', 'blah@blah'], stderr=stderr_file, stdout=stdout_file)
            elif coh_or_doc == 'doc':
                subprocess.call(['t_coffee', aln_file, '-output', 'fasta_aln', '-template_file',
                                 '~/data/t_coffee_data/t_coffee_doc.template', '-method', 'TMalign_pair',
                                 '-mode', '3dcoffee', '-email', 'blah@blah'], stderr=stderr_file, stdout=stdout_file)


def MakeJobs(arguments, args_outer, queue='fleishman'):
    import random
    import string
    global is_first_line_in_cmd
    pwd = os.getcwd()
    num = ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(8))
    jobname = pwd + '/job.' + num
    commandname = pwd + '/command'
    outname = pwd + '/out.' + num
    errname = pwd + '/err.' + num
    with open(jobname, 'wr+') as job:
        job.write('#!/bin/bash\n')
        job.write('. /usr/share/lsf/conf/profile.lsf\n')
        job.write('cd ' + pwd + '\n')
        if queue == 'new-all.q':
            job.write('/apps/RH6U4/blcr/0.8.5/bin/cr_run /home/labs/fleishman/jonathaw/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease ')
        else:
            job.write('/home/labs/fleishman/jonathaw/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease ')
        job.write(' '.join(arguments))
        job.write('\n')
    with open(commandname, 'a+') as command:
        if is_first_line_in_cmd:
            command.write('echo was submitted at `date` > ' + os.getcwd() + '/was_submitted\n')
            is_first_line_in_cmd = False
            # if args_outer['stop_at_purples']:
            #     command.write('echo ' + os.getcwd() + ' >> /home/labs/fleishman/jonathaw/general_lists/dirs_running_jobs\n')
        if queue == 'new-all.q':
            command.write(str('bsub -C 1024 -u /dev/null -N -u /dev/null -R rusage[mem=1024] -G fleishman-wx-grp-lsf -q ' +
                              queue + ' -o ' + outname + ' -e ' + errname + ' /apps/RH6U4/blcr/0.8.5/bin/cr_run ' +
                              jobname + '\n'))
        else:
            command.write(str('bsub -N -u /dev/null -G fleishman-wx-grp-lsf -q ' + queue + ' -o ' + outname + ' -e ' + errname
                          + ' ' + jobname + '\n'))
    subprocess.call(['chmod', '+x', jobname])


def seq_identity_by_seq(seq1, seq2):
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist

    matrix = matlist.blosum62
    gap_open = -8
    gap_extend = -4

    alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)
    top_aln = alns[0]
    aln_1, aln_2, score, begin, end = top_aln

    identical_pos = 0
    aligned_pos = 0
    internal_gap = 0
    for i in range(0, len(aln_1)):
        if aln_1[i] == aln_2[i] and aln_1[i] != '-':
            identical_pos += 1
            aligned_pos += 1
        elif aln_1[i] != '-' and aln_2[i] != '-':
            aligned_pos += 1
        elif aln_1[i] != aln_2[i] and aln_1[i] == '-' and aln_2[i] != '-':
            internal_gap += 1
        elif aln_1[i] != aln_2[i] and aln_2[i] == '-' and aln_1[i] != '-':
            internal_gap += 1
    return (100 * identical_pos) / (aligned_pos + internal_gap)


def original_names_from_list(names_list_file):
    results = {}
    with open(names_list_file, 'r') as file_list:
        for line in file_list:
            split_line = line.split()
            results[split_line[1]] = {'coh_wt': split_line[0].split('_')[0], 'doc_wt': split_line[0].split('_')[2]}
    return results


def ModelingBBs(args):
    coh_name = args['coh_name']
    doc_name = args['doc_name']
    if not args['coh_seq']:
        coh_seq = GetSequence(coh_name, args['coh_seq_db'])
    else:
        coh_seq = args['coh_seq']
    if not args['doc_seq']:
        doc_seq = GetSequence(doc_name, args['doc_seq_db'])
    else:
        doc_seq = args['doc_seq']
    coh_aln_file = args['coh_aln_file']
    doc_aln_file = args['doc_aln_file']
    JOB_NUM = args['job_num']
    shutil.copyfile(path_seqs_rb + 'all_coh.fasta', os.getcwd() + '/all_coh.fasta')
    shutil.copyfile(path_seqs_rb + 'all_doc.fasta', os.getcwd() + '/all_doc.fasta')
    shutil.copyfile(path_seqs_rb + 'all_RB.db', os.getcwd() + '/all_RB.db')
    print 'coh', coh_name, '\tseq is:', coh_seq
    print 'doc', doc_name, '\tseq is:', doc_seq
    AddAndAlign(coh_name, coh_seq, coh_aln_file, 'coh')
    AddAndAlign(doc_name, doc_seq, doc_aln_file, 'doc')
    chosen_cohs = 0
    for coh_i in COH_TEMPLATES:
        coh_iden = seq_identity_by_seq(coh_seq, GetSequence(coh_i, args['coh_aln_file']))
        if coh_iden < int(args['coh_iden_threshold']):
            print coh_i, 'rejected due to low sequecne identity', coh_iden
        else:
            chosen_cohs += 1
            print coh_i, 'chosen with sequence identity', coh_iden
            coh_RB_pos = RBposCoh(coh_aln_file, coh_i)
            coh_length = LengthSeqInAln(coh_aln_file, coh_i)
            for doc_i in DOC_TEMPLATES:
                doc_RB_pos = str(RBposDoc1st(doc_aln_file, doc_i) + coh_length)
                # doc_RB_pos_2nd = RBposDoc2nd(doc_aln_file, doc_i)
                for i in range(1, JOB_NUM):
                    if args['check_names']:
                        if coh_name.upper()[0:4] == coh_i.upper()[0:4] or doc_name.upper()[0:4] == doc_i.upper()[0:4]:
                            print 'query and template names identity, breaking!!!'
                            sys.exit(0)
                    job_arguments = ['-parser:protocol', '/home/labs/fleishman/jonathaw/protocols/' + args['protocol'],
                        '@/home/labs/fleishman/jonathaw/bin/modeling.flags',
                        '-database', '/home/labs/fleishman/jonathaw/Rosetta/main/database/',
                        '-parser:script_vars ' + 'RB_file=all_RB.db',
                        '-overwrite',
                        '-s', '/home/labs/fleishman/jonathaw/data/PDB_all_parts/' + coh_i.lower() + '_A.pdb',
                        '-parser:script_vars',
                        'doc_pdb=/home/labs/fleishman/jonathaw/data/PDB_all_parts/' + doc_i.lower() + '_B.pdb',
                        '-parser:script_vars', 'coh_aln_file=' + coh_aln_file + '_aln',
                        '-parser:script_vars', 'doc_aln_file=' + doc_aln_file + '_aln',
                        '-parser:script_vars', 'coh_query_name=' + coh_name.upper(),
                        '-parser:script_vars', 'doc_query_name=' + doc_name.upper(),
                        '-parser:script_vars', 'coh_template_name=' + coh_i + '.A',
                        '-parser:script_vars', 'doc_template_name=' + doc_i + '.B',
                        '-parser:script_vars', 'coh_rb_pos=' + str(coh_RB_pos) + 'A',
                        '-parser:script_vars', 'doc_rb_pos=' + str(doc_RB_pos) + 'B',
                        '-parser:script_vars', 'doc_start_res=' + str(coh_length + 1),
                        '-out:prefix', coh_name.lower() + '_' + doc_name.lower() + '_on_',
                        '-out:suffix', '_' + doc_i.lower() + '_' + str(i) + '.pdb'] #,
                        # '-out:file:scorefile',
                        # coh_name + '_' + doc_name + '_on_' + coh_i.lower() + '_' + doc_i.lower() + '_' + str(i) + '.sc']
                    if args['silent']:
                        job_arguments.extend(('-out:file:silent', coh_name.lower() + '_' + doc_name.lower() + '_on_' +
                        coh_name.lower() + '_A_' + doc_name.lower() + '_' + str(i) + '.out',
                                             '-out:file:silent_struct_type', 'binary'))
                    else:
                        job_arguments.extend(('-out:file:scorefile', coh_name + '_' + doc_name + '_on_' + coh_i.lower() +
                                             '_' + doc_i.lower() + '_' + str(i) + '.sc'))
                    MakeJobs(job_arguments, args_outer=args, queue=args['queue'])
    print "created", chosen_cohs * len(DOC_TEMPLATES) * JOB_NUM, "jobs, for", chosen_cohs, "cohesins"
    with open(os.getcwd() + '/name.for.files', 'wr+') as name_file:
        name_file.writelines(coh_name + '\t' + doc_name)
        print "created name.for.files with the names", coh_name, doc_name


def CohRBPosBySS(seq, NAME):
    if NAME.upper() == '1OHZ' or NAME.upper() == '1OHA':
        pat = re.compile("V-*F-*L-*F-*A.*")
    elif NAME.upper() == '2B59' or NAME.upper() == '2B5A':
        pat = re.compile("N-*F-*A-*L-*A.*")
    elif NAME.upper() == '2CCL' or NAME.upper().upper() == '2CCA':
        pat = re.compile("V-*F-*L-*F-*A.*")
    elif NAME.upper() == '2OZN' or NAME.upper() == '2OZA':
        pat = re.compile("R-*V-*L-*V-*S.*")
    elif NAME.upper() == '2VN5' or NAME.upper() == '2V5A':
        pat = re.compile("S-*F-*L-*F-*L.*")
    elif NAME.upper() == '2Y3N' or NAME.upper() == '2V6A':
        pat = re.compile("N-*F-*G-*R-*L.*")
    elif NAME.upper() == '3UL4' or NAME.upper() == '3ULA':
        pat = re.compile("A-*V-*L-*Y-*L.*")
    elif NAME.upper() == '3KCP' or NAME.upper() == '3KCA':
        pat = re.compile("N-*F-*A-*L-*A.*")
    elif NAME.upper() == '4DH2' or NAME.upper() == '4DHA':
        pat = re.compile("V-*F-*L-*F-*D.*")
    elif NAME.upper() == '4IU2' or NAME.upper() == '4IUA':
        pat = re.compile("F-*V-*A-*S-*G.*")
    elif NAME.upper() == '4UYP' or NAME.upper() == '4UPA':
        pat = re.compile("S-*M-*T-*F-*E.*")
    elif NAME.upper() == '4UYQ' or NAME.upper() == '4UQA':
        pat = re.compile("S-*M-*T-*F-*E.*")
    elif NAME.upper() == '4UMS':
        pat = re.compile("S-*F-*L-*F-*N.*")

    elif NAME.upper() == '1ANU':
        pat = re.compile("V-*F-*L-*F-*A.*")
    elif NAME.upper() == '1AOH':
        pat = re.compile("V-*F-*L-*F-*A.*")
    elif NAME.upper() == '1G1K':
        pat = re.compile("S-*F-*L-*F-*L.*")
    elif NAME.upper() == '1QZN':
        pat = re.compile("N-*F-*S-*K-*A.*")
    elif NAME.upper() == '1TYJ':
        pat = re.compile("N-*F-*G-*R-*L.*")
    elif NAME.upper() == '2BM3':
        pat = re.compile("N-*F-*A-*L-*A.*")
    elif NAME.upper() == '2JH2':
        pat = re.compile("R-*V-*L-*V-*S.*")
    elif NAME.upper() == '2VN6':
        pat = re.compile("S-*F-*L-*F-*L.*")
    elif NAME.upper() == '2VO8':
        pat = re.compile("R-*I-*L-*V-*A.*")
    elif NAME.upper() == '2W1N':
        pat = re.compile("R-*V-*L-*V-*S.*")
    elif NAME.upper() == '2W5F':
        pat = re.compile("A-*K-*Y-*K-*A.*")
    elif NAME.upper() == '2XBT':
        pat = re.compile("E-*I-*K-*L-*A.*")
    elif NAME.upper() == '2XDH':
        pat = re.compile("K-*V-*G-*I-*A.*")
    elif NAME.upper() == '2ZF9':
        pat = re.compile("F-*F-*T-*A-*T.*")
    elif NAME.upper() == '3BWZ':
        pat = re.compile("N-*F-*S-*K-*A.*")
    elif NAME.upper() == '3FNK':
        pat = re.compile("L-*S-*R-*S-*Y.*")
    elif NAME.upper() == '3L8Q':
        pat = re.compile("N-*L-*S-*R-*S.*")
    elif NAME.upper() == '4FL4':
        pat = re.compile("V-*F-*L-*F-*A.*")
    elif NAME.upper() == '4JO5':
        pat = re.compile("E-*I-*S-*F-*T.")
    elif NAME.upper() == '4FL5':
        pat = re.compile("N-*F-*A-*L-*A.*")
    elif NAME.upper() == '5NEW':
        pat = re.compile("S-*V-*I-*W-*S.*")
    return len(seq[0:pat.search(seq).start() + 1])


def DocRBPosBySS(seq, NAME):
    if NAME.upper() == '1OHZ' or NAME.upper() == '1OHA':
        pat = re.compile("S-*T-*D-*L-*T.*")
    elif NAME.upper() == '2B59':
        pat = re.compile("L-*D-*V-*A-*E.*")
    elif NAME.upper() == '2CCL':
        pat = re.compile("S-*T-*D-*L-*T.*")
    elif NAME.upper() == '2OZN':
        pat = re.compile("I-*G-*D-*L-*A.*")
    elif NAME.upper() == '2VN5':
        pat = re.compile("A-*L-*D-*F-*A.*")
    elif NAME.upper() == '2Y3N':
        pat = re.compile("M-*A-*D-*V-*M.*")
    elif NAME.upper() == '3UL4':
        pat = re.compile("D-*E-*D-*Y-*I.*")
    elif NAME.upper() == '3KCP':
        pat = re.compile("L-*L-*D-*V-*A.*")
    elif NAME.upper() == '4DH2':
        pat = re.compile("I-*S-*D-*Y-*V.*")
    elif NAME.upper() == '4IU2':
        pat = re.compile("G-*R-*D-*A-*T.*")
    elif NAME.upper() == '4UYP':
        pat = re.compile("S-*I-*D-*A-*V.*")
    elif NAME.upper() == '4UYQ':
        pat = re.compile("I-*N-*D-*A-*V.*")
    elif NAME.upper() == '2VN6':
        pat = re.compile("S-*T-*D-*F-*A.*")
    elif NAME.upper() == '4FL4':
        pat = re.compile("S-*T-*D-*L-*T.*")
    elif NAME.upper() == '4FL5':
        pat = re.compile("L-*L-*D-*V-*A.*")
    elif NAME.upper() == '5NEW':
        pat = re.compile("S-*D-*D-*L-*T.*")
    return len(seq[0:pat.search(seq).start() + 1])


def DesignJobs(args):
    coh_pssm_path = '/home/labs/fleishman/jonathaw/pssm/cohs/'
    doc_pssm_path = '/home/labs/fleishman/jonathaw/pssm/docs/'
    shutil.copyfile(path_seqs_rb + 'all_coh.fasta', os.getcwd() + '/all_coh.fasta')
    shutil.copyfile(path_seqs_rb + 'all_doc.fasta', os.getcwd() + '/all_doc.fasta')
    shutil.copyfile(path_seqs_rb + 'all_RB.db', os.getcwd() + '/all_RB.db')

    JOB_NUM = args['job_num']
    if args['design_mode'] == 'normal':
        protocol = 'DoCohDesigner_v3.xml'
    elif args['design_mode'] == 'weak_pssm':
        protocol = 'DoCohDesigner_weakPSSM.xml'
    elif args['design_mode'] == 'parts':
        protocol = 'DoCohDesigner_parts.xml'
    for coh_i in COH_TEMPLATES:
        coh_name = coh_i.lower()
        coh_seq = GetSequence(coh_name, args['coh_seq_db'])
        coh_RB_pos = CohRBPosBySS(coh_seq, coh_name)
        coh_length = LengthSeqInAln('all_coh.fasta', coh_name)
        for doc_i in DOC_TEMPLATES:
            doc_name = doc_i.lower()
            doc_seq = GetSequence(doc_name, args['doc_seq_db'])
            doc_RB_pos = str(DocRBPosBySS(doc_seq, doc_name) + coh_length)
            for i in range(0, JOB_NUM):
                job_arguments = ['-parser:protocol', '/home/labs/fleishman/jonathaw/protocols/' + protocol,
                        '@/home/labs/fleishman/jonathaw/bin/designer.flags',
                        '-database', '/home/labs/fleishman/jonathaw/Rosetta/main/database/',
                        '-parser:script_vars ' + 'RB_file=' + path_seqs_rb + 'all_RB.db',
                        '-overwrite',
                        '-s', '/home/labs/fleishman/jonathaw/data/PDB_all_parts/' + coh_name.lower() + '_A.pdb',
                        '-parser:script_vars',
                        'doc_pdb=/home/labs/fleishman/jonathaw/data/PDB_all_parts/' + doc_name.lower() + '_B.pdb',
                        '-parser:script_vars', 'coh_query_name=' + coh_name.upper(),
                        '-parser:script_vars', 'doc_query_name=' + doc_name.upper(),
                        '-parser:script_vars', 'coh_template_name=' + coh_name + '.A',
                        '-parser:script_vars', 'doc_template_name=' + doc_name + '.B',
                        '-parser:script_vars', 'coh_rb_pos=' + str(coh_RB_pos) + 'A',
                        '-parser:script_vars', 'doc_rb_pos=' + str(doc_RB_pos) + 'B',
                        '-parser:script_vars', 'coh_pssm=' + coh_pssm_path + coh_name + '_A.pssm',
                        '-parser:script_vars', 'doc_pssm=' + doc_pssm_path + doc_name + '_B.pssm',
                        '-out:suffix', '_' + doc_name.lower() + '_' + str(i) + '.pdb'] #,
                        # '-out:file:scorefile', coh_name + '_' + doc_name + '_on_' + coh_name.lower() + '_'
                        # + doc_name.lower() + '_' + str(i) + '.sc']
                if args['design_mode'] == 'parts':
                    res_type = {'pos': '"KR"', 'neg': '"DE"', 'no': '"ACDEFGHIKLMNPQRSTVWY"'}
                    pose_parts = {'no-loop': {'part_one': '115,116,37,39', 'part_two': '118,121,125,127'},
                                  'loop': {'part_one': '36,38,90,135,144', 'part_two': '33,151,152'}}
                    additional_args = []
                    if coh_name in ['1ohz', '1anu', '1aoh', '1g1k', '2ccl', '2jh2', '2ozn', '2vn5', '2vn6', '2vo8',
                                    '2w1n', '2w5f', '2xbt', '2xdh', '3ul4', '4dh2', '4fl4', '4jo5', '4ums', '4uyp',
                                    '4uyq']:
                        if args['part_one'] == 'pos':
                            additional_args.extend(['-parser:script_vars', 'poses_to_mute_1=115,116,37,39'])
                            additional_args.extend(['-parser:script_vars', 'res_types_1="KR"'])
                        if args['part_one'] == 'neg':
                            additional_args.extend(['-parser:script_vars', 'poses_to_mute_1=115,116,37,39'])
                            additional_args.extend(['-parser:script_vars', 'res_types_1="DE"'])
                        if args['part_one'] == 'no':
                            additional_args.extend(['-parser:script_vars', 'poses_to_mute_1=1'])
                            additional_args.extend(['-parser:script_vars', 'res_types_1="ACDEFGHIKLMNPQRSTVWY"'])
                        if args['part_two'] == 'pos':
                            additional_args.extend(['-parser:script_vars', 'poses_to_mute_2=118,121,125,127'])
                            additional_args.extend(['-parser:script_vars', 'res_types_2="KR"'])
                        if args['part_two'] == 'neg':
                            additional_args.extend(['-parser:script_vars', 'poses_to_mute_2=118,121,125,127'])
                            additional_args.extend(['-parser:script_vars', 'res_types_2="DE"'])
                        if args['part_two'] == 'no':
                            additional_args.extend(['-parser:script_vars', 'poses_to_mute_2=1'])
                            additional_args.extend(['-parser:script_vars', 'res_types_2="ACDEFGHIKLMNPQRSTVWY"'])
                        additional_args.extend(['-parser:script_vars', 'coh_type=/home/labs/fleishman/jonathaw/data/PDB_all_parts/1ohz_A.pdb'])
                    if coh_name in ['1qzn', '2b59', '2bm3', '2y3n', '2zf9', '3bwz', '3fnk', '3kcp', '3l8q', '4fl5',
                                    '4iu2', '5new', '1tyj']:
                        if args['part_one'] == 'pos':
                            additional_args.extend(['-parser:script_vars', 'poses_to_mute_1=36,38,90,135,144'])
                            additional_args.extend(['-parser:script_vars', 'res_types_1=K,L'])
                        if args['part_one'] == 'neg':
                            additional_args.extend(['-parser:script_vars', 'poses_to_mute_1=36,38,90,135,144'])
                            additional_args.extend(['-parser:script_vars', 'res_types_1=D,E'])
                        if args['part_one'] == 'no':
                            additional_args.extend(['-parser:script_vars', 'poses_to_mute_1='])
                            additional_args.extend(['-parser:script_vars', 'res_types_1='])
                        if args['part_two'] == 'pos':
                            additional_args.extend(['-parser:script_vars', 'poses_to_mute_2=33,151,152'])
                            additional_args.extend(['-parser:script_vars', 'res_types_2=K,L'])
                        if args['part_two'] == 'neg':
                            additional_args.extend(['-parser:script_vars', 'poses_to_mute_2=33,151,152'])
                            additional_args.extend(['-parser:script_vars', 'res_types_2=D,E'])
                        if args['part_two'] == 'no':
                            additional_args.extend(['-parser:script_vars', 'poses_to_mute_2='])
                            additional_args.extend(['-parser:script_vars', 'res_types_2='])
                        additional_args.extend(['-parser:script_vars', 'coh_type=/home/labs/fleishman/jonathaw/data/PDB_all_parts/2b59_A.pdb'])
                    job_arguments.extend(additional_args)
                if args['silent']:
                        job_arguments.extend(('-out:file:silent', coh_name.lower() + '_' + doc_name.lower() + '_on_' +
                        coh_name.lower() + '_A_' + doc_name.lower() + '_' + str(i) + '.out',
                                        '-out:file:silent_struct_type', 'binary'))
                else:
                    job_arguments.extend(('-out:file:scorefile', coh_name + '_' + doc_name + '_on_' + coh_i.lower() +
                                        '_' + doc_i.lower() + '_' + str(i) + '.sc'))
                MakeJobs(job_arguments, args, queue=args['queue'])
    with open(os.getcwd() + '/name.for.files', 'wr+') as name_file:
        name_file.writelines('design\tdesign')
        print "creating name.for.files with the names design design"
    add_to_design_list()


def design_ct_single_part(args):
    import random
    coh_pssm_path = '/home/labs/fleishman/jonathaw/pssm/cohs/'
    doc_pssm_path = '/home/labs/fleishman/jonathaw/pssm/docs/'
    shutil.copyfile(path_seqs_rb + 'all_coh.fasta', os.getcwd() + '/all_coh.fasta')
    shutil.copyfile(path_seqs_rb + 'all_doc.fasta', os.getcwd() + '/all_doc.fasta')
    shutil.copyfile(path_seqs_rb + 'all_RB.db', os.getcwd() + '/all_RB.db')

    res_type = {'pos': '"KR"', 'neg': '"DE"', 'pssm': 'ACDEFGHIKLMNPQRSTVWY'}
    # pose_parts = [37, 39, 115, 116, 118, 121, 125, 127]
    COH_TEMPLATES = ['1OHZ', '2CCL', '1ANU', '1AOH']

    JOB_NUM = 1
   # for coh_i in COH_TEMPLATES:
   #    coh_name = coh_i.lower()
   #     coh_seq = GetSequence(coh_name, args.coh_seq_db)
   #     coh_RB_pos = CohRBPosBySS(coh_seq, coh_name)
   #     coh_length = LengthSeqInAln('all_coh.fasta', coh_name)

       # for doc_i in DOC_TEMPLATES:
       #     doc_name = doc_i.lower()
       #     doc_seq = GetSequence(doc_name, args.doc_seq_db)
       #     doc_RB_pos = str(DocRBPosBySS(doc_seq, doc_name) + coh_length)
    for i in range(0, JOB_NUM):
        for pos33 in res_type.values():
            for pos35 in res_type.values():
                for pos37 in res_type.values():
                    for pos63 in res_type.values():
                        for pos66 in res_type.values():
                            for pos70 in res_type.values():
                                for pos73 in res_type.values():
                                    for pos83 in res_type.values():
                                        for pos116 in res_type.values():
                                            for pos127 in res_type.values():
                                                coh_i = random.choice(COH_TEMPLATES)
                                                coh_name = coh_i.lower()
                                                coh_seq = GetSequence(coh_name, args['coh_seq_db'])
                                                coh_RB_pos = CohRBPosBySS(coh_seq, coh_name)
                                                coh_length = LengthSeqInAln('all_coh.fasta', coh_name)

                                                doc_i = random.choice(DOC_TEMPLATES)
                                                doc_name = doc_i.lower()
                                                doc_seq = GetSequence(doc_name, args['doc_seq_db'])
                                                doc_RB_pos = str(DocRBPosBySS(doc_seq, doc_name) + coh_length)
                                                job_arguments = ['-parser:protocol', '/home/labs/fleishman/jonathaw/protocols/DoCohDesigner_Ct_byParts_2nd.xml',
                                                    '@/home/labs/fleishman/jonathaw/bin/designer.flags',
                                                    '-database', '/home/labs/fleishman/jonathaw/Rosetta/main/database/',
                                                    '-parser:script_vars ' + 'RB_file=' + path_seqs_rb + 'all_RB.db',
                                                    '-overwrite',
                                                    '-s', '/home/labs/fleishman/jonathaw/data/PDB_all_parts/' + coh_name.lower() + '_A.pdb',
                                                    '-parser:script_vars',
                                                    'doc_pdb=/home/labs/fleishman/jonathaw/data/PDB_all_parts/' + doc_name.lower() + '_B.pdb',
                                                    '-parser:script_vars', 'coh_query_name=' + coh_name.upper(),
                                                    '-parser:script_vars', 'doc_query_name=' + doc_name.upper(),
                                                    '-parser:script_vars', 'coh_template_name=' + coh_name + '.A',
                                                    '-parser:script_vars', 'doc_template_name=' + doc_name + '.B',
                                                    '-parser:script_vars', 'coh_rb_pos=' + str(coh_RB_pos) + 'A',
                                                    '-parser:script_vars', 'doc_rb_pos=' + str(doc_RB_pos) + 'B',
                                                    '-parser:script_vars', 'coh_pssm=' + coh_pssm_path + coh_name + '_A.pssm',
                                                    '-parser:script_vars', 'doc_pssm=' + doc_pssm_path + doc_name + '_B.pssm',
                                                    '-out:suffix', '_' + doc_name.lower() + '_' + str(i) + '.pdb',
                                                    '-parser:script_vars', 'res_types_33=' + pos33,
                                                    '-parser:script_vars', 'res_types_35=' + pos35,
                                                    '-parser:script_vars', 'res_types_37=' + pos37,
                                                    '-parser:script_vars', 'res_types_63=' + pos63,
                                                    '-parser:script_vars', 'res_types_66=' + pos66,
                                                    '-parser:script_vars', 'res_types_70=' + pos70,
                                                    '-parser:script_vars', 'res_types_73=' + pos73,
                                                    '-parser:script_vars', 'res_types_83=' + pos83,
                                                    '-parser:script_vars', 'res_types_116=' + pos116,
                                                    '-parser:script_vars', 'res_types_127=' + pos127,
                                                    '-out:file:silent', coh_name.lower() + '_' + doc_name.lower() + '_on_' + coh_name.lower() + '_A_' + doc_name.lower() + '_' + str(i) + '.out',
                                                    '-out:file:silent_struct_type', 'binary']
                                                MakeJobs(job_arguments, args, queue=args['queue'])
    with open(os.getcwd() + '/name.for.files', 'wr+') as name_file:
        name_file.writelines('design\tdesign')
        print "creating name.for.files with the names design design"
    add_to_design_list()


def Design_sigmoids_jobs(args):
    wt_VS_ds = {'dsn5': '1aoa_dsn5_on_1anu_A_2vn6_3.pdb_0001.pdb', 'dsn7': '1qza_dsn7_on_2b59_A_2y3n_7.pdb_0005.pdb',
                'dsn9': '2v6a_dsn9_on_2vn6_A_2vn5_2.pdb_0004.pdb', 'ds10': '3kca_ds10_on_4fl5_A_4fl5_6.pdb_0001.pdb',
                'ds11': '3kca_ds11_on_2b59_A_2y3n_2.pdb_0002.pdb', 'ds14': '2bma_ds14_on_2bm3_A_3ul4_11.pdb_0003.pdb',
                'ds15': '2bma_ds15_on_2bm3_A_3ul4_8.pdb_0003.pdb', 'ds17': '2y3a_ds17_on_2y3n_A_2y3n_3.pdb_0002.pdb',
                'ds18': '4f5a_ds18_on_3kcp_A_2y3n_11.pdb_0003.pdb', 'ds20': '4uqa_ds20_on_4uyp_A_2vn6_8.pdb_0001.pdb',
                'ds21': '4uqa_ds21_on_4uyp_A_4fl4_13.pdb_0001.pdb', 'ds22': '5nea_ds22_on_5new_A_4dh2_11.pdb_0003.pdb',
                'ds23': '5nea_ds23_on_5new_A_2vn5_3.pdb_0005.pdb', 'ds24': '5nea_ds24_on_5new_A_1ohz_6.pdb_0002.pdb',
                'ds25': '5nea_ds25_on_5new_A_2vn6_10.pdb_0001.pdb', 'ds26': '5nea_ds26_on_5new_A_2vn6_13.pdb_0005.pdb',
                'dz34': '4uqa_dz34_on_4uyp_A_5new_8.pdb_0002.pdb', 'dz35': '2b5a_dz35_on_3kcp_A_4iu2_13.pdb_0002.pdb',
                'dz42': '4f4a_dz42_on_2ccl_A_4fl4_3.pdb_0002.pdb'}
    ds_VS_wt = {'dsn5': 'dsn5_2v6a_on_2ccl_A_3ul4_7.pdb_0005.pdb', 'dsn7': 'dsn7_2y3a_on_2y3n_A_3kcp_11.pdb_0004.pdb',
                'dsn9': 'dsn9_2v5a_on_2ccl_A_4dh2_1.pdb_0002.pdb', 'ds10': 'ds10_4f5a_on_3bwz_A_3kcp_6.pdb_0005.pdb',
                'ds11': 'ds11_4f5a_on_2b59_A_3kcp_10.pdb_0003.pdb', 'ds14': 'ds14_5nea_on_1qzn_A_5new_9.pdb_0001.pdb',
                'ds15': 'ds15_5nea_on_3kcp_A_3kcp_4.pdb_0005.pdb', 'ds17': 'ds17_2y3a_on_1tyj_A_2y3n_8.pdb_0003.pdb',
                'ds18': 'ds18_5nea_on_4fl5_A_5new_5.pdb_0002.pdb', 'ds20': 'ds20_4uqa_on_4uyp_A_4uyq_11.pdb_0004.pdb',
                'ds21': 'ds21_4uqa_on_4uyp_A_4uyq_4.pdb_0001.pdb', 'ds22': 'ds22_3ula_on_5new_A_2ccl_10.pdb_0003.pdb',
                'ds23': 'ds23_4dha_on_5new_A_2ccl_8.pdb_0002.pdb', 'ds24': 'ds24_4upa_on_5new_A_4dh2_6.pdb_0005.pdb',
                'ds25': 'ds25_5nea_on_5new_A_2ccl_1.pdb_0002.pdb', 'ds26': 'ds26_5nea_on_5new_A_5new_10.pdb_0002.pdb',
                'dz34': 'dz34_4dha_on_4uyq_A_4uyq_6.pdb_0005.pdb', 'dz35': 'dz35_4iua_on_3bwz_A_4iu2_13.pdb_0004.pdb',
                'dz42': 'dz42_4f4a_on_2ccl_A_2ccl_3.pdb_0003.pdb'}
    for i in range(0, args['job_num']):
        job_arguments = ['-parser:protocol', '/home/labs/fleishman/jonathaw/protocols/DoCohDesigner_sigmoids.xml',
                        '@/home/labs/fleishman/jonathaw/bin/designer.flags',
                        '-database', '/home/labs/fleishman/jonathaw/Rosetta/main/database/',
                        '-overwrite',
                        '-s', '/home/labs/fleishman/jonathaw/designs_21.1/' + args['design_name'] + '.pdb',
                        '-parser:script_vars', 'coh_pssm=' + coh_pssm_path + args['coh_name'] + '_A.pssm',
                        '-parser:script_vars', 'doc_pssm=' + doc_pssm_path + args['doc_name'] + '_B.pssm',
                        '-parser:script_vars', 'wt_VS_ds=' +
                        '/home/labs/fleishman/jonathaw/designs_21.1/reference_poses_new/' + wt_VS_ds[args['design_name']],
                        '-parser:script_vars', 'ds_VS_wt=' +
                        '/home/labs/fleishman/jonathaw/designs_21.1/reference_poses_new/' + ds_VS_wt[args['design_name']],
                        '-out:file:silent', args['design_name'] + '_sigmoids_' + str(i) + '.out',
                        '-out:file:silent_struct_type', 'binary']
        if args['dump_pose']:
            job_arguments.extend(('-parser:script_vars DP_wt_VS_ds=' + os.getcwd() + '/' + wt_VS_ds[args['design_name']] +
            '_DP_' + str(i) + '.pdb',
            '-parser:script_vars DP_ds_VS_wt=' + os.getcwd() + '/' + ds_VS_wt[args['design_name']] +
            '_DP_' + str(i) + '.pdb'))
        MakeJobs(job_arguments, args, queue=args['queue'])
    with open(os.getcwd() + '/name.for.files', 'wr+') as name_file:
        name_file.writelines(args['design_name'] + '\t' + 'sigmoids')
        print "creating name.for.files with the names", args['design_name'], 'sigmoids'
    add_to_design_list()


def AppendToCommandList(queue, priority):
    cmd_in_list = False
    if queue == 'new-all.q':
        command_list = '/home/labs/fleishman/jonathaw/bin/command_list_new_all_q'
    else:
        command_list = '/home/labs/fleishman/jonathaw/bin/command_list'
    with open(command_list, 'r') as r:
        for line in r:
            if line == os.getcwd() + '/command\n':
                print 'found command in command list, will not append'
                cmd_in_list = True
                return
    if not cmd_in_list:
        if not priority:
            with open(command_list, 'a') as f:
                f.write(os.getcwd() + '/command\n')
                print 'added to the', queue, 'command list queue'
        else:
            with open(command_list, 'r') as f:
                origin = []
                for line in f:
                    origin.append(line)
            with open(command_list, 'w') as f:
                f.write(os.getcwd() + '/command\n')
                for line in origin:
                    f.write(line)
            print "at user's request, added command to", command_list, "as priority"


def JackKnifing(args):
    coh_name = args['coh_name']
    doc_name = args['doc_name']
    if not args['coh_seq']:
        coh_seq = GetSequence(coh_name, args['coh_seq_db'])
    else:
        coh_seq = args['coh_seq']
    if not args['doc_seq']:
        doc_seq = GetSequence(doc_name, args['doc_seq_db'])
    else:
        doc_seq = args['doc_seq']
    coh_aln_file = args['coh_aln_file']
    doc_aln_file = args['doc_aln_file']
    JOB_NUM = args['job_num']
    shutil.copyfile(path_seqs_rb + 'all_coh.fasta', os.getcwd() + '/all_coh.fasta')
    shutil.copyfile(path_seqs_rb + 'all_doc.fasta', os.getcwd() + '/all_doc.fasta')
    shutil.copyfile(path_seqs_rb + 'all_RB.db', os.getcwd() + '/all_RB.db')
    print 'coh', coh_name, '\tseq is:', coh_seq
    print 'doc', doc_name, '\tseq is:', doc_seq
    AddAndAlign(coh_name, coh_seq, coh_aln_file)
    AddAndAlign(doc_name, doc_seq, doc_aln_file)
    chosen_cohs = 0
    for coh_i in COH_TEMPLATES:
        coh_iden = seq_identity_by_seq(coh_seq, GetSequence(coh_i, args['coh_aln_file']))
        if coh_iden < int(args['coh_iden_threshold']):
            print coh_i, 'rejected due to low sequecne identity', coh_iden
        else:
            chosen_cohs += 1
            print coh_i, 'chosen with sequence identity', coh_iden
            coh_RB_pos = RBposCoh(coh_aln_file, coh_i)
            coh_length = LengthSeqInAln(coh_aln_file, coh_i)
            for doc_i in DOC_TEMPLATES:
                doc_RB_pos = str(RBposDoc1st(doc_aln_file, doc_i) + coh_length)
                # doc_RB_pos_2nd = RBposDoc2nd(doc_aln_file, doc_i)
                for i in range(1, JOB_NUM):
                    job_arguments = ['-parser:protocol', '/home/labs/fleishman/jonathaw/protocols/' +
                                     'DoCohModeller_JackKnifing.xml',
                        '@/home/labs/fleishman/jonathaw/bin/modeling.flags',
                        '-database', '/home/labs/fleishman/jonathaw/Rosetta/main/database/',
                        '-parser:script_vars ' + 'RB_file=all_RB.db',
                        '-overwrite',
                        '-s', '/home/labs/fleishman/jonathaw/data/PDB_all_parts/' + coh_i.lower() + '_A.pdb',
                        '-parser:script_vars',
                        'doc_pdb=/home/labs/fleishman/jonathaw/data/PDB_all_parts/' + doc_i.lower() + '_B.pdb',
                        '-parser:script_vars', 'coh_aln_file=' + coh_aln_file + '_aln',
                        '-parser:script_vars', 'doc_aln_file=' + doc_aln_file + '_aln',
                        '-parser:script_vars', 'coh_query_name=' + coh_name.upper(),
                        '-parser:script_vars', 'doc_query_name=' + doc_name.upper(),
                        '-parser:script_vars', 'coh_template_name=' + coh_i + '.A',
                        '-parser:script_vars', 'doc_template_name=' + doc_i + '.B',
                        '-parser:script_vars', 'coh_rb_pos=' + str(coh_RB_pos) + 'A',
                        '-parser:script_vars', 'doc_rb_pos=' + str(doc_RB_pos) + 'B',
                        '-parser:script_vars', 'doc_start_res=' + str(coh_length + 1),
                        '-out:prefix', coh_name.lower() + '_' + doc_name.lower() + '_on_',
                        '-out:suffix', '_' + doc_i.lower() + '_' + str(i) + '.pdb',
                        '-native', '/home/labs/fleishman/jonathaw/data/PDB_all/' + coh_name[0:4] + '_AB.pdb']
                    if args['silent']:
                        job_arguments.extend(('-out:file:silent', coh_name.lower() + '_' + doc_name.lower() + '_on_' +
                        coh_name.lower() + '_A_' + doc_name.lower() + '_' + i + '.out',
                                             '-out:file:silent_struct_type', 'binary'))
                    else:
                        job_arguments.extend(('-out:file:scorefile', coh_name + '_' + doc_name + '_on_' + coh_i.lower() +
                                             '_' + doc_i.lower() + '_' + str(i) + '.sc'))
                    MakeJobs(job_arguments, args, queue=args['queue'])
    print "created", chosen_cohs * len(DOC_TEMPLATES) * JOB_NUM, "jobs, for", chosen_cohs, "cohesins"


def design_prediction_by_list(args):
    wt_names = original_names_from_list(args['names_list'])
    pwd_parent = os.getcwd()
    for design in wt_names:
        os.makedirs(pwd_parent + '/' + design)
        os.chdir(pwd_parent + '/' + design)
        pwd = os.getcwd()
        coh_wt = names_Ad[wt_names[design]['coh_wt']]
        doc_wt = names_Ad[wt_names[design]['doc_wt']]
        ### makes the design_VS_design dir, and produces its jobs
        os.makedirs(pwd + '/' + design + '_VS_' + design)
        os.chdir(pwd + '/' + design + '_VS_' + design)
        args['coh_name'] = design
        args['doc_name'] = design
        ModelingBBs(args)
        if args['add_queue'] == 1:
            AppendToCommandList(queue=args['queue'], priority=args['priority'])
        else:
            print "at user's request, did not add to queue", args['queue']
        ### makes the design_VS_wt dir, and produces its jobs
        os.makedirs(pwd + '/' + design + '_VS_' + doc_wt)
        os.chdir(pwd + '/' + design + '_VS_' + doc_wt)
        is_first_line_in_cmd = True
        args['coh_name'] = design
        args['doc_name'] = doc_wt
        ModelingBBs(args)
        if args['add_queue'] == 1:
            AppendToCommandList(queue=args['queue'], priority=args['priority'])
        else:
            print "at user's request, did not add to queue", args['queue']
        ### makes the wt_VS_design dir, and produces its jobs
        os.makedirs(pwd + '/' + coh_wt + '_VS_' + design)
        os.chdir(pwd + '/' + coh_wt + '_VS_' + design)
        is_first_line_in_cmd = True
        args['coh_name'] = coh_wt
        args['doc_name'] = design
        ModelingBBs(args)
        if args['add_queue'] == 1:
            AppendToCommandList(queue=args['queue'], priority=args['priority'])
        else:
            print "at user's request, did not add to queue", args['queue']
        os.chdir(pwd_parent)
    args['add_queue'] = 0


def design_prediciton(args):
    global is_first_line_in_cmd
    coh_wt = {'dz31': '5nea', 'dz32': '5nea', 'dz33': '5nea', 'dz34': '4uqa', 'dz35': '2b5a', 'dz36': '1ana',
              'dz37': '1oha', 'dz38': '2y3a', 'dz39': '1ana', 'dz40': '2cca', 'dz41': '1qza', 'dz42': '4f4a',
              'dz43': '4dha', 'dz44': '4uqa', 'dz45': '2bma', 'dz46': '1qza', 'dz51': '1aoa', 'dz52': '1aoa',
              'dz53': '1aoa', 'dz54': '1qza', 'dz55': '2b5a', 'dz56': '2b5a', 'dz57': '2b5a', 'dz58': '2b5a',
              'dz59': '2bma', 'dz60': '2bma', 'dz61': '2bma', 'dz62': '2bma', 'dz63': '2v6a', 'dz64': '3kca',
              'dz65': '3kca', 'dz66': '3kca', 'dz67': '4f5a', 'dz68': '4upa', 'dz69': '4upa', 'dz70': '4uqa',
              'dz71': '5nea', 'dz72': '5nea', 'dz73': '5nea'}
    doc_wt = {'dz31': '2cca', 'dz32': '2cca', 'dz33': '4upa', 'dz34': '4dha', 'dz35': '4iua', 'dz36': '2v5a',
              'dz37': '2v5a', 'dz38': '4f5a', 'dz39': '2v6a', 'dz40': '4f4a', 'dz41': '2cca', 'dz42': '4f4a',
              'dz43': '2v5a', 'dz44': '4f4a', 'dz45': '5nea', 'dz46': '3kca', 'dz51': '2v6a', 'dz52': '3ula',
              'dz53': '4dha', 'dz54': '4f5a', 'dz55': '2y3a', 'dz56': '3kca', 'dz57': '4f5a', 'dz58': '4f5a',
              'dz59': '2y3a', 'dz60': '2y3a', 'dz61': '2y3a', 'dz62': '4f5a', 'dz63': '2v6a', 'dz64': '2y3a',
              'dz65': '2y3a', 'dz66': '3kca', 'dz67': '4uqa', 'dz68': '1oha', 'dz69': '4uqa', 'dz70': '2cca',
              'dz71': '1oha', 'dz72': '4iua', 'dz73': '4upa'
                                                      ''}
    pwd = os.getcwd()
    ### makes the design_VS_design dir, and produces its jobs
    os.makedirs(pwd + '/' + args['design_name'] + '_VS_' + args['design_name'])
    os.chdir(pwd + '/' + args['design_name'] + '_VS_' + args['design_name'])
    args['coh_name'] = args['design_name']
    args['doc_name'] = args['design_name']
    ModelingBBs(args)
    if args['add_queue'] == 1:
        AppendToCommandList(queue=args['queue'], priority=args['priority'])
    else:
        print "at user's request, did not add to queue", args['queue']
    ### makes the design_VS_wt dir, and produces its jobs
    os.makedirs(pwd + '/' + args['design_name'] + '_VS_' + doc_wt[args['design_name']])
    os.chdir(pwd + '/' + args['design_name'] + '_VS_' + doc_wt[args['design_name']])
    is_first_line_in_cmd = True
    args['coh_name'] = args['design_name']
    args['doc_name'] = doc_wt[args['design_name']]
    ModelingBBs(args)
    if args['add_queue'] == 1:
        AppendToCommandList(queue=args['queue'], priority=args['priority'])
    else:
        print "at user's request, did not add to queue", args['queue']
    ### makes the wt_VS_design dir, and produces its jobs
    os.makedirs(pwd + '/' + coh_wt[args['design_name']] + '_VS_' + args['design_name'])
    os.chdir(pwd + '/' + coh_wt[args['design_name']] + '_VS_' + args['design_name'])
    is_first_line_in_cmd = True
    args['coh_name'] = coh_wt[args['design_name']]
    args['doc_name'] = args['design_name']
    ModelingBBs(args)
    if args['add_queue'] == 1:
        AppendToCommandList(queue=args['queue'], priority=args['priority'])
    else:
        print "at user's request, did not add to queue", args['queue']
    args['add_queue'] = 0


def add_to_design_list():
    with open('/home/labs/fleishman/jonathaw/general_lists/design_filders', 'wr+') as f:
        f.write(os.getcwd())


def main():
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument('-coh_name', default=None)
    parser.add_argument('-doc_name', default=None)
    parser.add_argument('-coh_seq', default=False)
    parser.add_argument('-doc_seq', default=False)
    parser.add_argument('-coh_aln_file', default="./all_coh.fasta")
    parser.add_argument('-doc_aln_file', default="./all_doc.fasta")
    parser.add_argument('-add_queue', default=1, type=int)
    parser.add_argument('-job_num', default=14, type=int)
    parser.add_argument('-queue', default='fleishman', type=str)
    parser.add_argument('-mode', default='predict', type=str)
    parser.add_argument('-protocol', default='DoCohModeller.xml', type=str)
    parser.add_argument('-coh_seq_db', default='/home/labs/fleishman/jonathaw/data/seqs_all/cohesins_from_rachel.fasta',
                        type=str)
    parser.add_argument('-doc_seq_db', default='/home/labs/fleishman/jonathaw/data/seqs_all/dockerins_from_rachel.fasta',
                        type=str)
    parser.add_argument('-coh_iden_threshold', default=27, type=int)
    parser.add_argument('-priority', default=False, type=bool)
    parser.add_argument('-silent', default=True, type=bool)
    parser.add_argument('-design_name', default=None, type=str)
    parser.add_argument('-dump_pose', default=False, type=bool)
    #used in MakeJobs to test if the query and template names are identical:
    parser.add_argument('-check_names', default=True, type=bool)
    parser.add_argument('-stop_at_purples', default=False, type=bool)
    parser.add_argument('-design_mode', default='normal')
    parser.add_argument('-part_one', type=str)
    parser.add_argument('-part_two', type=str)
    parser.add_argument('-names_list', type=str)

    args = parser.parse_args()
    args_dict = vars(args)

    if args.mode == 'predict':
        ModelingBBs(args_dict)
    elif args.mode == 'design':
        DesignJobs(args_dict)
    elif args.mode == 'sigmoids':
        Design_sigmoids_jobs(args_dict)
    elif args.mode == 'jack_knifing':
        JackKnifing(args_dict)
    elif args.mode == 'design_prediction':
        design_prediciton(args_dict)
    elif args.mode == 'design_prediction_by_list':
        design_prediction_by_list(args_dict)
    elif args.mode == 'desing_Ct_parts':
        design_ct_single_part(args_dict)

    # if args['add_queue'] == 1:
    #     AppendToCommandList(queue=args['queue'], priority=args['priority'])
    # else:
    #     print "at user's request, did not add to queue", args['queue']


if __name__ == '__main__':
    main()
