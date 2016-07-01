"""
this script is for truncating Rafael's gigantic PDB (multiple coh-doc pairs) into separate files, in separate folders,
with BB atoms only. this is for reading RB later.
"""
import fileinput, os
folder_num = 1
pdb_num = 1
bb_atoms = ['N', 'CA', 'C', 'O']
current_pdb = []
for line in fileinput.input('1ohz.dynamics.pdb'):
    split = line.split()
    if split[0] == 'END':
        if not os.path.exists('./dir_' + str(pdb_num/1000) + '/'):
            os.makedirs('./dir_' + str(pdb_num/1000) + '/')
        with open('./dir_%i/1ohz_%i.pdb' % (pdb_num/1000, pdb_num), 'wr+') as o:
            for w_line in current_pdb:
                o.writelines(w_line)
        current_pdb = []
        pdb_num += 1
        # break
    elif split[2] in bb_atoms:
        current_pdb.append(line)


