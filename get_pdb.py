def download_and_get_chains():
    from Bio.PDB import PDBParser, PDBIO
    failed = []
    pdbs_dict = read_rostdb_entries()
    io = PDBIO()
    pdbl = PDBList()
    for pdb_e, chains in pdbs_dict.items():
        for chain_e in chains:
            try:
                pdbl.retrieve_pdb_file(pdb_e, pdir='./')
                pdb = PDBParser().get_structure(pdb_e, 'pdb'+pdb_e.lower()+'.ent')
                for chain in pdb.get_chains():
                    if chain.get_id() == chain_e:
                        io.set_structure(chain)
                        io.save(pdb.get_id() + '_' + chain.get_id() + '.pdb')
            except:
                failed.append((pdb_e, chain_e))
    print "failures:", failed


def download_pdb():
    from Bio.PDB import PDBParser, PDBIO, PDBList
    io = PDBIO()
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(args['name'], pdir=args['path'])


if __name__ == '__main__':
    import argparse
    import sys
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', type=str, default='single')
    parser.add_argument('-name', type=str, default=sys.argv[2])
    parser.add_argument('-chain', type=str, default=sys.argv[3])
    parser.add_argument('-path', type=str, default='./')
    args = vars(parser.parse_args())
    if args['mode'] == 'single':
        download_pdb()