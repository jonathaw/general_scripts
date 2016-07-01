def seq_from_pdb(pdb_file, chain):
    aas_3_to_1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                  'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
                  'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}
    with open(pdb_file, 'r') as f:
        seq = ''
        for line in f:
            # print line
            line_split = line.split()
            if line_split[0] == 'ATOM' and int(line_split[5]) != len(seq) and line_split[4].upper() == chain.upper():
                seq += aas_3_to_1[line_split[3]]
        return seq


def chain_seq_from_pdb(pdb_file, chain):
    aas_3_to_1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                  'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
                  'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}
    with open(pdb_file, 'r') as f:
        seqs = {}
        up_to_now = 0
        for line in f:
            line_split = line.split()
            if line_split[0] == 'ATOM' and int(line_split[5]) > up_to_now:
                if line_split[4].upper() not in seqs.keys():
                    seqs[line_split[4].upper()] = ''
                seqs[line_split[4].upper()] += aas_3_to_1[line_split[3]]
                up_to_now += 1
    return seqs[chain.upper()]

if __name__ == '__main__':
    import sys
    seq = chain_seq_from_pdb(sys.argv[1], sys.argv[2])
    print seq