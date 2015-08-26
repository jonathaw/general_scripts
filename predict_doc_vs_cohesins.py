def doc_vs_cohesins(args):
    import os
    cwd = os.getcwd()
    # names_of_desings = names_of_cohesins(args)
    for coh_name in ['ct11', 'ct12', 'ct13', 'ct15', 'ct16', 'ct17', 'ct29', 'ct31', 'ct33', 'ct36', 'ct38', 'ct41', 'ct44',
                 'ct45', 'ct46', 'ct47', 'ct48', 'ct49', 'ct50', 'ct51', 'ct52', 'ct53', 'ct54', 'ct55', 'ct59']:
        os.chdir(cwd)
        make_prediction(args, coh_name, doc_name=args['doc_name'])


def cohs_vs_docs(args):
    import os
    cwd = os.getcwd()
    for coh in args['coh_list']:
        for doc in args['doc_list']:
            # olga's condition:
            type_i = ['4212_c6', '4212_c1', '4212_c1_ext', '1560_c1', '1560_c3', '5583_c1']
            type_ii = ['1973_c7', '1973_c10', '5233_c4', '2768_c3', '2768_c1', '4497_c2']
            no_type_i_ghs = ['GH8_0751_doc', 'GH43_3654_doc', '4717_GH30_Xdoc', 'GH9_2036_Xdoc', 'GH48_doc', 'GH5_0447_doc',
                   'GH10_4354_doc', 'GH5_2998_doc', 'GH30_1029_doc', 'GH5_2375_Xdoc']
            no_type_ii = ['1973_doc', '2768_doc', 'ScaA_doc']
            if (coh in type_i and doc in no_type_i_ghs) or (coh in type_ii and doc in no_type_ii):
                continue
            os.chdir(cwd)
            print 'making prediction for', coh, doc
            make_prediction(args, coh, doc)


def predict_diagonal(args):
    import os
    cwd = os.getcwd()
    for coh in args['coh_list']:
        os.chdir(cwd)
        print 'making prediction for ', coh, coh
        make_prediction(args, coh, coh)


def names_of_cohesins(args):
    names = []
    with open(args['coh_seq_db'], 'r') as f:
        for line in f:
            if line[0] == '>' and line[1:5].upper() not in ['1OHZ', '2B59', '2CCL', '2OZN', '2VN5', '2VN6', '2Y3N', '3KCP', '3UL4', '4DH2', '4IU2', '4UYP', '4UYQ',
                 '4FL4', '1ANU', '1AOH', '1G1K', '1QZN', '1TYJ', '2BM3', '2JH2', '2VO8', '2W1N', '2W5F', '2XBT', '2XDH',
                 '2ZF9', '3BWZ', '3FNK', '3L8Q', '4JO5', '4UMS', '4FL5', '5NEW']:
                names.append(line[1:].rstrip().split('.')[0])
    return names


def make_prediction(args, coh_name, doc_name):
    import os
    from DoCohJobMaker_options import ModelingBBs
    from DoCohJobMaker_options import AppendToCommandList
    cwd = os.getcwd()
    os.mkdir(cwd + '/' + coh_name + '_VS_' + doc_name)
    os.chdir(cwd + '/' + coh_name + '_VS_' + doc_name)
    modeling_args = {'coh_name': coh_name, 'doc_name': doc_name, 'job_num': 14, 'coh_seq': False,
                     'doc_seq': False, 'coh_seq_db': args['coh_seq_db'], 'doc_seq_db': args['doc_seq_db'],
                     'coh_aln_file': './all_coh.fasta', 'doc_aln_file': './all_doc.fasta', 'coh_iden_threshold': 27,
                     'check_names': True, 'protocol': 'DoCohModeller.xml', 'silent': True, 'queue': args['queue'],
                     'stop_at_purples': True}
    ModelingBBs(modeling_args)
    AppendToCommandList(modeling_args['queue'], False)


if __name__ == '__main__':
    import argparse
    import sys
    saveout = sys.stdout
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', default='doc_vs_cohesins', type=str)
    parser.add_argument('-doc_name', type=str)
    parser.add_argument('-doc_seq_db', type=str)
    parser.add_argument('-coh_seq_db', type=str)
    parser.add_argument('-queue', default='fleishman', type=str)
    parser.add_argument('-coh_list', nargs='+', type=str)
    parser.add_argument('-doc_list', nargs='+', type=str)
    parser.add_argument('-log_file', default='auto_job_maker.log', type=str)
    args_main = vars(parser.parse_args())
    logger = open(args_main['log_file'], 'wr+')
    sys.stdout = logger
    if args_main['mode'] == 'doc_vs_cohesins':
        doc_vs_cohesins(args_main)
    if args_main['mode'] == 'cohs_vs_docs':
        cohs_vs_docs(args_main)
    if args_main['mode'] == 'diagonal':
        predict_diagonal(args_main)
    sys.stdout = saveout
    logger.close()