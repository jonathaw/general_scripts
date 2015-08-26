coh_names = ['ct11', 'ct12', 'ct13', 'ct15', 'ct16', 'ct17', 'ct29', 'ct31', 'ct33', 'ct36', 'ct38', 'ct41',
                   'ct44', 'ct45', 'ct46', 'ct47', 'ct48', 'ct49', 'ct50', 'ct51', 'ct52', 'ct53', 'ct54', 'ct55',
                   'ct59']
doc_names = ['ct11', 'ct12', 'ct13', 'ct15', 'ct16', 'ct17', 'ct29', 'ct31', 'ct33', 'ct36', 'ct38', 'ct41',
                   'ct44', 'ct45', 'ct46', 'ct47', 'ct48', 'ct49', 'ct50', 'ct51', 'ct52', 'ct53', 'ct54', 'ct55',
                   'ct59']

def main():
    import numpy as np
    from rosetta_score_files import how_many_purples_in_file
    import os
    import re
    global coh_names, doc_names
    design_list = ['ct11', 'ct12', 'ct13', 'ct15', 'ct16', 'ct17', 'ct29', 'ct31', 'ct33', 'ct36', 'ct38', 'ct41',
                   'ct44', 'ct45', 'ct46', 'ct47', 'ct48', 'ct49', 'ct50', 'ct51', 'ct52', 'ct53', 'ct54', 'ct55',
                   'ct59']
    # df = DataFrame({name: Series([-1 * len(design_list)], index=design_list) for name in design_list})
    mtrx = np.zeros([len(design_list), len(design_list)], dtype=int)
    score_file_list = [x for x in os.listdir(
        '/Users/jonathan/eden/no_backup/designs/Ct_8parts_10.2/prediction/results/')
                       if re.match('.*\.score', x)]
    for score_file in score_file_list:
        coh_name = score_file.split('_')[1]
        doc_name = score_file.split('_')[3]
        purple_num = int(how_many_purples_in_file(
            '/Users/jonathan/eden/no_backup/designs/Ct_8parts_10.2/prediction/results/' + score_file))
        mtrx[design_list.index(coh_name)][design_list.index(doc_name)] = 1 if purple_num >= 12 else 0
    dof_vec = find_degree_vector(mtrx)
    print dof_vec
    print mtrx
    mtrx = clean_all_zeros(mtrx)
    while not are_all_ones(dof_vec):
        new_all_ones(dof_vec)
        to_remove = dof_vec[-1].values()[0]
        mtrx = remove_from_matrix(mtrx, to_remove)
        mtrx = clean_all_zeros(mtrx)
        dof_vec = find_degree_vector(mtrx)
        # break
        print mtrx
        print dof_vec
        print coh_names
        print doc_names


def find_degree_vector(matrix):
    cohs = [{val: (i, 'coh')} for i, val in enumerate(matrix.sum(axis=1))]
    docs = [{val: (i, 'doc')} for i, val in enumerate(matrix.sum(axis=0))]
    cohs.extend(docs)
    return sorted(cohs)


def clean_all_zeros(matrix):
    row_sum = matrix.sum(axis=1)
    col_sum = matrix.sum(axis=0)
    coh_removed = 0
    doc_removed = 0
    for i, row in enumerate(row_sum):
        if row == 0:
            matrix = remove_from_matrix(matrix, (i-coh_removed, 'coh'))
            coh_removed += 1
    for i, col in enumerate(col_sum):
        if col == 0:
            matrix = remove_from_matrix(matrix, (i-doc_removed, 'doc'))
            doc_removed += 1
    return matrix


def remove_from_matrix(matrix, to_remove):
    import numpy as np
    if to_remove[1] == 'coh':
        coh_names.remove(coh_names[to_remove[0]])
    else:
        doc_names.remove(doc_names[to_remove[0]])
    return np.delete(matrix, to_remove[0], 1 if to_remove[1] == 'doc' else 0)


def are_all_ones(dof_vec):
    return True if sum([a.keys()[0] for a in dof_vec]) == len(dof_vec) else False


def new_all_ones(dfv):
    print dfv
    print all([a.keys()[0] == 1 for a in dfv])


if __name__ == '__main__':
    main()