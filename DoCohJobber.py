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




def main():
    import argparse
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
