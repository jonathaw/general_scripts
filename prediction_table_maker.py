def main():
    from pandas import DataFrame, Series
    from rosetta_score_files import how_many_purples_in_file
    import os
    import re
    design_list = ['ct11', 'ct12', 'ct13', 'ct15', 'ct16', 'ct17', 'ct29', 'ct31', 'ct33', 'ct36', 'ct38', 'ct41',
                   'ct44', 'ct45', 'ct46', 'ct47', 'ct48', 'ct49', 'ct50', 'ct51', 'ct52', 'ct53', 'ct54', 'ct55',
                   'ct59']
    # df = DataFrame({name: Series([-1 * len(design_list)], index=design_list) for name in design_list})
    score_file_list = [x for x in os.listdir('.') if re.match('.*\.score', x)]

    coh_name_list = sorted(list(set(['_'.join(a.split('_VS_')[0].split('_')[1:]) for a in score_file_list])))
    doc_name_list = sorted(list(set(['_'.join(a.split('_VS_')[1].split('_')[:-1]) for a in score_file_list])))
    # print coh_name_list
    # print doc_name_list
    df = DataFrame({coh_name: Series([-1 * len(doc_name_list)], index=doc_name_list) for coh_name in coh_name_list})

    for score_file in score_file_list:
        # coh_name = score_file.split('_')[1]
        # doc_name = score_file.split('_')[3]
        coh_name = '_'.join(score_file.split('_VS_')[0].split('_')[1:])
        doc_name = '_'.join(score_file.split('_VS_')[1].split('_')[:-1])
        # print coh_name, doc_name
        purple_num = int(how_many_purples_in_file(score_file))
        df[coh_name][doc_name] = purple_num
    # pandas.set_option('display.max_columns', None)
    # print df
    show_prediction_heat_map(df.copy())


def show_prediction_heat_map(df):
    import matplotlib.pyplot as plt
    from matplotlib import colors
    from numpy import array, arange
    cmap = colors.ListedColormap(['white', 'red', 'green', 'purple'])
    bounds = [-100, 0, 10, 15, 100]
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


if __name__ == '__main__':
    main()