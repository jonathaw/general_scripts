#!/usr/bin/env python
def analyse_doc(args):
    import seq_manager as sm
    import os
    query = sm.WorkSeq(args['name'][:4], args['name'], sm.read_multi_fastas('all_docs.fasta')[args['name']+'.B']['seq'])
    helices = get_helices_seq(args['name'][:4], query.fasta)
    all_doc_fastas = sm.read_multi_fastas('all_dockerins_May2015.txt')


    pass_thresh = [{'name': args['name'], 'seq': query.fasta}]
    for hit_name, hit_val in all_doc_fastas.items():
        pw_aln = sm.pair_wise_aln_from_seqs(query.fasta, hit_val['seq'].upper())
        if sm.aln_identity(pw_aln[0], pw_aln[1]) >= 0.4 and not gap_in_essential(helices, pw_aln[1]) \
                and not gap_in_essential(helices, pw_aln[0]):
            pass_thresh.append(hit_val)

    pass_thresh_fasta_name = args['path']+args['name']+'_passed_thresholds.fasta'
    sm.write_multi_seqs_to_file({a['name']: a for a in pass_thresh}, pass_thresh_fasta_name)
    os.system('muscle -in ' + pass_thresh_fasta_name + ' -out ' + args['path']+args['name']+'_passed_thresholds.aln')
    aln_fastas = sm.read_multi_fastas(args['path']+args['name']+'_passed_thresholds.aln')
    os.system('weblogo -f %s -o %s --annotate %s' % (args['path']+args['name']+'_passed_thresholds.aln',
                                                      args['path']+args['name']+'_passed_thresholds.eps',
                                                      ','.join(a for a in aln_fastas[args['name']]['seq'].replace('-', '-'))))
    print 'aln_fastas', aln_fastas
    sm.write_seq_from_seq(args['name'], args['path'], query.fasta)
    sm.write_multi_seqs_to_file(aln_fastas, args['path']+args['name']+'_passed_thresholds.aln', args['name'])
    sm.run_psi_blast_for_pssm(args['path']+args['name']+'.fasta', args['path']+args['name']+'_passed_thresholds.aln',
                              args['path']+args['name']+'.pssm')


def analyse_coh(args):
    import seq_manager as sm
    import os
    query = sm.WorkSeq(args['name'][:4], args['name'], sm.read_multi_fastas('all_cohs.txt')[args['name']+'.A']['seq'])
    print args
    print query
    strands = get_strands_seq(query.name, query.fasta)
    all_coh_fastas = sm.read_multi_fastas('Coh_All_Setup-Oct2014-d-for-lizi.txt')

    pass_threshold = [{'name': query.name, 'seq': query.fasta}]
    for hit_name, hit_val in all_coh_fastas.items():
        # print hit_val
        pw_aln = sm.pair_wise_aln_from_seqs(query.fasta, hit_val['seq'].upper())
        if sm.aln_identity(pw_aln[0], pw_aln[1]) >= 0.4 and not gap_in_essential(strands, pw_aln[1], 4) \
                and not gap_in_essential(strands, pw_aln[0], 4):
            pass_threshold.append(hit_val)
        # else:
        #     print 'not accepted'
        #     print pw_aln[0]
        #     print pw_aln[1]
        #     print ''.join(a if i in strands else '_' for i, a in enumerate(query.fasta))
    print 'found %i seqs that passed the thresholds' % len(pass_threshold)
    # import sys
    # sys.exit()
    pass_thresh_fasta_name = args['path']+args['name']+'_passed_thresholds.fasta'
    sm.write_multi_seqs_to_file(pass_threshold, pass_thresh_fasta_name)
    os.system('muscle -in ' + pass_thresh_fasta_name + ' -out ' + args['path']+args['name']+'_passed_thresholds.aln')
    aln_fastas = sm.read_multi_fastas(args['path']+args['name']+'_passed_thresholds.aln')
    os.system('weblogo -f %s -o %s --annotate %s' % (args['path']+args['name']+'_passed_thresholds.aln',
                                                      args['path']+args['name']+'_passed_thresholds.eps',
                                                      ','.join(a for a in aln_fastas[args['name']]['seq'].replace('-', '-'))))
    sm.write_seq_from_seq(args['name'], args['path'], query.fasta)
    sm.run_psi_blast_for_pssm(args['path']+args['name']+'.fasta', args['path']+args['name']+'_passed_thresholds.aln',
                              args['path']+args['name']+'.pssm')


def get_strands_seq(name, seq):
    """
    :param name: name
    :param seq: whole sequence
    :return:a list of all positions in seq that are within the known strands, found manually in the PDBs
    """
    known_strands = {'1ohz': ['ANCDFVFR', 'KSFDTAIYPDRKIIVFLFAEDSG', 'ITFDEVGGFADNDLVE'],
                     '2b59': ['AGFQVNIV', 'GPIQIADNDPEKGILNFALAYSYI', 'QDTLSMPGAILGTQLFDWDGEVI'],
                     '3ul4': ['VSSDFVIEY', 'ESFSYNVVEKDEIIAVLYLEETGL', 'GISPIKFESFGATADNDMNEMT'],
                     '4dh2': ['QSFNLSLYYDSK', 'NNFDYNIVYKDSEIVFLFDDDKQ', 'TKKYSLITFGESNFCDFDL'],
                     '4fl4': ['ASGDFVVSY', 'KSFDTAVYPDRKMIVFLFAEDSG', 'GFSAIEISEFGAFADNDLVE']}
    result = []
    result.extend(range(seq.find(known_strands[name][0]), seq.find(known_strands[name][0])+len(known_strands[name][0])))
    result.extend(range(seq.find(known_strands[name][1]), seq.find(known_strands[name][1])+len(known_strands[name][1])))
    result.extend(range(seq.find(known_strands[name][2]), seq.find(known_strands[name][2])+len(known_strands[name][2])))
    return result


def get_helices_seq(name, seq):
    """
    :param name: name
    :param seq: whole sequnce as found in the script
    :return: a list of all positions is seq that are part of the known helices
    looking at the PDBs for each entry, I designated the helices manually for each one. seeing that many end in a
    positive charge I added these.
    >>> name = '1ohz'
    >>> seq = 'ESSSVLLGDVNGDGTINSTDLTMLKRSVLRAITLTDDAKARADVDKNGSINSTDVLLLSRYLLRVIDKFPVAENP'
    >>> test_helices_seq(name, seq)
    [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62]
    """
    known_helices = {'1ohz': ['DVNGDGTINSTDLTMLKRSVLR', 'DVDKNGSINSTDVLLLSRYLLR'],
                     '2b59': ['DNSINLLDVAEVIRCF', 'DINRNGAINMQDIMIVHKHF'],
                     '2ozn': ['DFNGNGKIDIGDLAMVSKNI', 'DLNKDGSIDEYEISFINHRILNLE'],
                     '2vn5': ['ALDFAGLKKYIM', 'STDLAILKKYLLGMV'],
                     '2vn6': ['DYNNDGNVDSTDFAGLKKYIMAAD', 'DVNLDNEVNAFDLAILKKYLLGM'],
                     '2y3n': ['DLNGDGVINMADVMILAQSFGK', 'DLNNDGVINSDDAIILAQYFGK'],
                     '3ul4': ['DLNRNGIVNDEDYILLKNYLLR', 'DVNKDGKVNSTDCLFLKKYIL'],
                     '4dh2': ['DVNADGVVNISDYVLMKRYILR', 'DVNGDNVINDIDCNYLKRYLL'],
                     '4fl4': ['DVNDDGKVNSTDLTLLKRYVLK', 'DVNRDGRVNSSDVTILSRYLIR'],
                     '4fl5': ['DIVKDNSINLLDVAEVIRCF', 'DINRNGAINMQDIMIVHK'],
                     '4iu2': ['DTDLNNIVDGRDATATLTYYAATST', 'DGRDASSILTFYTKSSV'],
                     '4uyp': ['DVDGNGSVRSIDAVLIRDYVLGK', 'DVDGNGSIKINDAVLVRDYVLGK'],
                     '4uyq': ['DVDGNGSVRSIDAVLIRDYVLGK', 'DVDGNGSIKINDAVLVRDYVLGK'],
                     '5new': ['DLDGDGEVDVFDLILMRKAVE', 'DLNCDGVIDSDDLTYHSEYLH']}
    result = []
    result.extend(range(seq.find(known_helices[name][0]), seq.find(known_helices[name][0])+len(known_helices[name][0])))
    result.extend(range(seq.find(known_helices[name][1]), seq.find(known_helices[name][1])+len(known_helices[name][1])))
    return result


def within_essentials(essentials, seq):
    return ''.join(a for i, a in enumerate(seq) if i in essentials)


def gap_in_essential(essentials, aln_seq, M=1):
    """
    :param helices: a list of all positions determined to be essential
    :param aln_seq: an laigned sequence
    :return: True iff all locations within essential are non gaps
    >>> essentials = [1, 2, 3]
    >>> aln_seq1 = 'abcdefg'
    >>> aln_seq2 = 'ab-defg'
    >>> gap_in_essential(essentials, aln_seq1)
    False
    >>> gap_in_essential(essentials, aln_seq2)
    True
    """
    res = 0
    for loc in essentials:
        if aln_seq[loc] == '-':
            res += 1
    if res >= M:
        return True
    else:
        return False



if __name__ == '__main__':
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument('-name', type=str)
    parser.add_argument('-path', default=os.getcwd()+'/', type=str)
    parser.add_argument('-mode', default=None, type=str)
    args = vars(parser.parse_args())
    if args['mode'] == 'analyse_doc':
        analyse_doc(args)
    elif args['mode'] == 'analyse_coh':
        analyse_coh(args)
    elif args['mode'] is None:
        print 'no mode was chosen'
    # args = {'name': '1ohz_B', 'path': '/Volumes/jonathaw/doc_analysis_16Aug/'}
