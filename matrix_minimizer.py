def main():
    from pandas import DataFrame, Series
    from rosetta_score_files import how_many_purples_in_file
    import os
    import re
    # design_list = ['ct11', 'ct12', 'ct13', 'ct15', 'ct16', 'ct17', 'ct29', 'ct31', 'ct33', 'ct36', 'ct38', 'ct41',
    #                'ct44', 'ct45', 'ct46', 'ct47', 'ct48', 'ct49', 'ct50', 'ct51', 'ct52', 'ct53', 'ct54', 'ct55',
    #                'ct59']
    # df = DataFrame({name: Series([-1], index=design_list) for name in design_list})
#    score_file_list = [x for x in os.listdir('/Users/jonathan/eden/no_backup/designs/Ct_8parts_10.2/prediction/results/')
#                       if re.match('.*\.score', x)]
    score_file_list = [x for x in os.listdir('./')
                       if re.match('.*\.score', x)]
    coh_name_list = sorted(list(set([a.split('_')[1] for a in score_file_list])))
    doc_name_list = sorted(list(set([a.split('_')[3] for a in score_file_list])))
    df = DataFrame({coh_name: Series([-1 * len(doc_name_list)], index=doc_name_list) for coh_name in coh_name_list})
    for score_file in score_file_list:
        coh_name = score_file.split('_')[1]
        doc_name = score_file.split('_')[3]
#        purple_num = int(how_many_purples_in_file('/Users/jonathan/eden/no_backup/designs/Ct_8parts_10.2/prediction/results/'+score_file))
#        df[coh_name][doc_name] = 1 if purple_num >= 10 else 0
        purple_num = int(how_many_purples_in_file('./'+score_file))
        df[coh_name][doc_name] = 1 if purple_num >= 10 else 0

    i = 1
    while not all_dof_ones(df):
        dof_vec = find_degree_vector(df)
        df = remove_from_df(df, dof_vec[-1].values()[0])

        df = clean_zeroes(df)
        print 'printing dof for %i time' % i
        print dof_vec[-1]
        print find_degree_vector(df)
        if i > -1:
            show_prediction_heat_map(df)
        i += 1


def find_degree_vector(dfi):
    from pandas import DataFrame as DF
    results = [{dof: (dfi.columns[i], 'coh')} for i, dof in enumerate(DF.sum(dfi))]
    temp = [{dof: (dfi.index[i], 'doc')} for i, dof in enumerate(DF.sum(dfi, axis=1))]
    results.extend(temp)
    # results = [{i[0]: (i[1], 'coh')} for i in DF.sum(dfi)]
    return sorted(results)


def remove_from_df(dfi, to_remove):
    return dfi.drop(to_remove[0], axis=1 if to_remove[1] == 'coh' else 0)


def clean_zeroes(dfi):
    dof_vec_i = find_degree_vector(dfi)
    print '\n cleaners dof', dof_vec_i
    for dif in dof_vec_i:
        if dif.keys()[0] == 0:
            print 'cleaner is removing', dif
            dfi = remove_from_df(dfi, dif.values()[0])
    return dfi


def show_prediction_heat_map(df):
    import matplotlib.pyplot as plt
    from matplotlib import colors
    from numpy import array, arange
    cmap = colors.ListedColormap(['white', 'red', 'green', 'purple'])
    bounds = [-100, 0, 1, 15, 100]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    heatmap = plt.pcolor(array(df), cmap=cmap, norm=norm, edgecolors='k')
    for y in range(array(df.shape)[0]):
        for x in range(array(df.shape)[1]):
            if array(df)[y, x] >= 0:
                plt.text(x+0.5, y+0.5, array(df)[y, x], horizontalalignment='center', verticalalignment='center')
    plt.yticks(arange(0.5, len(df.index), 1), df.index)
    plt.xticks(arange(0.5, len(df.columns), 1), df.columns, rotation='vertical')
    plt.xlabel('cohesin name')
    plt.ylabel('dockerin name')
    plt.title('cohesin dockerin cross binding')
    plt.show()


def all_dof_ones(dfi):
    dof_i = find_degree_vector(dfi)
    for i in dof_i:
        print 'all deg is one found %i' % i.keys()[0]
        if i.keys()[0] > 1:
            return False
    return True


if __name__ == '__main__':
    main()