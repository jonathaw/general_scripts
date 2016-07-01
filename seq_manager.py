#!/usr/bin/env python
class WorkSeq():
    def __init__(self, name, dir, fasta=None, evalue=0.0001, blast_identity=0.4, blast_percent_query_len=0.9, cdhit_threshold=0.97):
        self.name = name
        self.dir = dir
        self.blast_evalue = evalue
        self.blast_identity = blast_identity
        self.blast_percent_query_len = blast_percent_query_len
        self.blast_name = name + '_blast.xml'
        self.fasta_name = name + '.fasta'
        self.blast_fastas = name + '_blast.fasta'
        self.cdhit_name = name + '.cdhit'
        self.msa_name = name + '.msa'
        self.cdhit_threshold = cdhit_threshold
        if fasta is None:
            self.fasta = read_fasta(self)
        else:
            self.fasta = fasta
        self.length = len(self.fasta)

    def __repr__(self):
        msg = ''
        for k, v in self.__dict__.items():
            msg += str(k) + ' ' + str(v) + '\n'
        return msg

def get_seq(name):
    from Bio import PDB
    parser = PDB.PDBParser()
    struct = parser.get_structure(name, name)
    ppd = PDB.PPBuilder()
    peptides = ppd.build_peptides(struct)
    seq = ''.join([str(pep.get_sequence()) for pep in peptides])
    return seq


def write_seq_form_file(name, in_path, out_path):
    seq = get_seq(in_path + name + '.pdb')
    out_path_name = out_path + name + '.fasta'
    with open(out_path_name, 'wr+') as o:
        o.write('>%s\n' % name.split('/')[-1].split('.')[0])
        o.write('%s\n' % seq)
    return out_path_name


def write_seq_from_seq(name, out_path, seq):
    out_path_name = out_path + name + '.fasta'
    with open(out_path_name, 'wr+') as o:
        o.write('>%s\n' % name.split('/')[-1].split('.')[0])
        o.write('%s\n' % seq)
    return out_path_name


def write_multi_seqs_to_file(seqs, file_i, query=None):
    """
    :param seqs: {name: {'name': 'a', 'seq': 'AAA'}} dict
    :param file_i: file with path and name
    :param query: a wuery to be first in the file
    :return: writes a file
    """
    print 'writer recived', seqs
    with open(file_i, 'wr+') as fout:
        if query is not None:
            fout.write('>%s\n' % query)
            fout.write('%s\n' % seqs[query]['seq'])
        for s in seqs.values():
            if query is not None:
                if query == s['name']:
                    continue
            fout.write('>%s\n' % s['name'])
            fout.write('%s\n' % s['seq'])


def pair_wise_aln_from_pdbs(name1, name2, matrix=None, gap_open=-10, gap_extend=-0.5):
    """
    :param name1: a pdb name
    :param name2: a pdb name
    :param matrix: a substitution matrix
    :param gap_open: gap_open penalty
    :param gap_extend: gap_extend penalty
    :return: aln for seq1, aln for seq2, score, beginning and end for the best alignment
    """
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
    if matrix is None:
        matrix = matlist.blosum62
    seq1 = get_seq(name1)
    seq2 = get_seq(name2)
    alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)[0]
    return alns


def pair_wise_aln_from_seqs(seq1, seq2, matrix=None, gap_open=-10, gap_extend=-0.5):
    """
    :param seq1: a fasta sequence
    :param seq2: a fasta sequence
    :param matrix: a substitution matrix
    :param gap_open: gap_open penalty
    :param gap_extend: gap_extend penalty
    :return: aln for seq1, aln for seq2, score, beginning and end for the best alignment
    """
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
    if matrix is None:
        matrix = matlist.blosum62
    alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)[0]
    return alns


def blastp(name_full, name, in_path, out_path, evalue, max_target_seqs):
    from os import system
    suffix = '_blast.xml'
    write_seq_from_file(name, in_path, out_path)
    print ('blastp -query %s '
           '-db /shareDB/nr/nr '
           '-evalue %f '
           '-max_target_seqs %i '
           '-outfmt 5 '
           '-out %s ' % (out_path+name+'.fasta', evalue, max_target_seqs, out_path+name+suffix))

    system('blastp -query %s '
           '-db /shareDB/nr/nr '
           '-evalue %f '
           '-max_target_seqs %i '
           '-outfmt 5 '
           '-out %s ' % (out_path+name+'.fasta', evalue, max_target_seqs, out_path+name+suffix))
    return out_path+name+suffix


def blastp2seqs(name, in_path, out_path):
    out_path_name = out_path+name+'_seqs.fasta'
    print "reading %s, writing sequences to %s" % (in_path+name, out_path_name)
    with open(in_path+name+'.xml', 'r') as f:
        xml = f.read().split('<Hit>')
    blast_reuslts = []
    for hit in xml[1:]:
        result = {}
        split = hit.split('\n')
        for line in split:
            split_tri = line.split('<')
            if len(split_tri) < 3: continue
            if split_tri[2] == '/Hit_id>': result['name'] = split_tri[1].split('>')[1]
            if split_tri[2] == '/Hsp_hseq>': result['seq'] = split_tri[1].split('>')[1]
        blast_reuslts.append(result)
    with open(out_path_name, 'wr+')as o:
        for res in blast_reuslts:
            o.write('>%s\n' % res['name'])
            o.write('%s\n' % res['seq'])
    return out_path_name


def parse_blast_xml(ws):
    """
    :param ws: a WorkSeq object
    :return: makes a file with WorkSeq's name with fastas from the blast. includes the query as the first entry
    """
    from Bio.Blast.NCBIXML import read
    xml_handle = open(ws.blast_name, 'r')
    xml = read(xml_handle)
    fout = open(ws.blast_fastas, 'wr+')
    fout.write(">%s\n" % ws.name)
    fout.write("%s\n" % ws.fasta)
    for aln in xml.alignments:
        for hsp in aln.hsps:
            # print vars(hsp)  # use to see what's in hsp
            coverage = float(hsp.align_length) / float(ws.length) #  Adi's coverage.
            identity = float(hsp.identities) / float(hsp.align_length)
            evalue = float(hsp.expect)
            percent_query_len = (float(hsp.query_end)-float(hsp.query_start)) / float(ws.length)
            seq = str(hsp.sbjct).translate(None, "-")
            if percent_query_len > ws.blast_percent_query_len and identity > ws.blast_identity \
                    and evalue < ws.blast_evalue and ws.length*0.9 <= len(seq) <= ws.length*1.1:
                fout.write(">%s_%f\n" % (aln.accession, identity))
                fout.write("%s\n" % seq)
                # print "%s %f" % (aln.accession, percent_query_len)
    xml_handle.close()
    fout.close()


def cd_hit_cluster(ws):
    from os import system
    system('cd-hit -i %s -o %s -c %f ' % (ws.dir+ws.blast_fastas, ws.dir+ws.cdhit_name, ws.cdhit_threshold))


def read_fasta(ws):
    with open(ws.fasta_name, 'r') as f:
        cont = f.read().split('>')
    for line in cont:
        split = line.split()
        if len(split) < 2:
            continue
        if split[0].replace('_', '.') == ws.name.replace('_', '.'):
            return split[1]


def read_multi_fastas(file):
    with open(file, 'r') as f:
        cont = f.read().split('>')
    result = {}
    for entry in cont:
        split_entry = entry.split('\n')
        if len(split_entry) < 2:
            continue
        name = '_'.join(split_entry[0].rstrip().split())
        seq = ''.join(a.rstrip() for a in split_entry[1:])
        result[name] = {'name': name, 'seq': seq}
    return result


def insure_query_in_cd_hit(ws):
    cdhit_fastas = read_multi_fastas(ws.dir+ws.cdhit_name)
    with open(ws.cdhit_name, 'w') as fout:
        fout.write(">%s\n" % ws.name)
        fout.write("%s\n" % ws.fasta)
        for hit in cdhit_fastas.values():
            if hit['name'] == ws.name:
                continue
            else:
                fout.write(">%s\n" % hit['name'])
                fout.write("%s\n" % hit['seq'])


def aln_identity(aln1, aln2):
    """
    :param aln1: alignment sequence (with gaps)
    :param aln2: alignment sequence (with gaps)
    :return: the identity calculated by: (# identities)/(aln1 length, no gaps)
    >>> a = '-ABC-D'
    >>> b = '-ABCED'
    >>> aln_identity(a, b)
    1.0
    >>> b = '-BBCED'
    >>> aln_identity(a, b)
    0.75
    """
    res = 0.
    for i, aa in enumerate(aln1):
        res += 1. if aln1[i] == aln2[i] != '-' else 0.
    length = float(len(aln1.replace('-', '')))
    return res/length


def run_muscle(ws):
    import os
    os.system('muscle -in ' + ws.dir+ws.cdhit_name + ' -out ' + ws.dir+ws.msa_name)


def run_psi_blast_for_pssm(query_file, msa_file, out_file):
    import os
    os.system('psiblast -subject %s -in_msa %s -out_ascii_pssm %s' % (query_file, msa_file, out_file))

def main():
    import argparse
    import os
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument('-name', type=str)  # nargs='+'
    parser.add_argument('-name2', default=None, type=str)
    parser.add_argument('-mode', default='get_seq', type=str)
    parser.add_argument('-chain', default='all')
    parser.add_argument('-out_path', default=os.getcwd()+'/')
    parser.add_argument('-blastp_evalue', default=0.0001, type=float)
    parser.add_argument('-blastp_max_target_seqs', default=1500, type=int)
    parser.add_argument('-cd_hit_threshold', default=0.97, type=float)
    args = vars(parser.parse_args())
    args['name_full'] = args['name']
    args['name'] = args['name_full'].split('/')[-1].split('.')[0]
    args['name_full2'] = args['name2']
    if args['name2'] is not None:
        args['name2'] = args['name_full2'].split('/')[-1].split('.')[0]
    if args['mode'] == 'get_seq':
        print get_seq(args['name'])
    args['in_path'] = '/'.join(args['name_full'].split('/')[:-1]) + '/'
    if args['mode'] == 'pairwise':
        aln1, aln2, score, begin, end = pair_wise_aln(args['name_full'], args['name_full2'])
        print "%-10s %s" % (args['name'], aln1)
        print "%-10s %s" % (args['name2'], aln2)
        print "score %f" % score
    if args['mode'] == 'blastp':
        blastp(args['name_full'], args['name'], args['in_path'], args['out_path'], args['blastp_evalue'],
               args['blastp_max_target_seqs'])
    if args['mode'] == 'blastp2seqs':
        blastp2seqs(args['name'], args['in_path'], args['out_path'])
    if args['mode'] == 'cd_hit':
        cd_hit_cluster(args['name']+'.fasta', args['in_path'], args['out_path'], args['cd_hit_threshold'])
    if args['mode'] == 'read_fasta':
        print read_fasta(args['name'], args['in_path'])
    if args['mode'] == 'insure_query_in_cd_hit':
        insure_query_in_cd_hit(args['name'], args['in_path'])


def new_main():
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument('-name', type=str)
    parser.add_argument('-path', default=os.getcwd()+'/')
    parser.add_argument('-mode', type=str)
    parser.add_argument('-cdhit_threshold', default=0.9, type=float)
    parser.add_argument('-blast_evalue', default=0.0001, type=float)
    parser.add_argument('-blast_identity', default=0.4, type=float)
    parser.add_argument('-blast_percent_query_len', default=0.9, type=float)

    args = vars(parser.parse_args())
    ws = WorkSeq(args['name'], args['path'], args['blast_evalue'], args['blast_identity'],
                 args['blast_percent_query_len'], args['cdhit_threshold'])
    if args['mode'] == 'blast2fasta':
        parse_blast_xml(ws)
    elif args['mode'] == 'cdhit':
        cd_hit_cluster(ws)
    elif args['mode'] == 'insure_query_in_cdhit':
        insure_query_in_cd_hit(ws)
    elif args['mode'] == 'muscle':
        run_muscle(ws)


if __name__ == '__main__':
    # ws = WorkSeq('q9kx51', './')# p0ab98
    new_main()
    # parse_blast_xml(ws)
    #print "ws %i" % ws.length