#!/usr/bin/env python3
import os
import glob
import time
import argparse
import pickle
import matplotlib as mpl
# mpl.rc_context(fname="/home/labs/fleishman/jonathaw/.matplotlib/" +
               # "publishable_matplotlibrc")
# mpl.use('Agg')
import matplotlib.pyplot as plt
from sklearn import linear_model
# from sklearn.preprocessing import Imputer
from sklearn.metrics import r2_score
from matplotlib.widgets import Slider
import RosettaFilter as Rf
from Logger import Logger
import pandas as pd
import numpy as np
# from AnnotateFinder import AnnoteFinder
# from FollowDotCursor import FollowDotCursor
filters = ['total', 'a_sasa']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', default='q')
    parser.add_argument('-pc_old')
    parser.add_argument('-score_old')
    parser.add_argument('-obj_old')
    parser.add_argument('-pc_new')
    parser.add_argument('-score_new')
    parser.add_argument('-obj_new')
    parser.add_argument('-pc_wj')
    parser.add_argument('-score_wj')
    parser.add_argument('-obj_wj')
    parser.add_argument('-sc')
    parser.add_argument('-mp_sc')
    parser.add_argument('-rs_sc')
    parser.add_argument('-pc')
    parser.add_argument('-mp_pc')
    parser.add_argument('-rs_pc')
    parser.add_argument('-span_threshold', default=0.3, type=float)
    parser.add_argument('-names', default=None, help='file list names to show')
    parser.add_argument('-log_path', default='./', help='path to for log file')
    parser.add_argument('-percent', type=int, default=100)
    parser.add_argument('-best', type=bool, default=False)
    parser.add_argument('-terms', nargs='+',
                        default=['score', 'a_shape', 'a_pack', 'a_ddg',
                                 'res_solv'])
    parser.add_argument('-threshold', type=int, default=5)
    parser.add_argument('-show', default='show')
    parser.add_argument('-all4', default=False)
    parser.add_argument('-dir', help='which folder to use for fnd results')
    parser.add_argument('-original_sc')
    parser.add_argument('-name')

    args = vars(parser.parse_args())
    args['logger'] = Logger('logeer_%s.log' % time.strftime("%d.%0-m"),
                            args['log_path'])

    if args['mode'] == 'old':
        analyse_old(args)

    elif args['mode'] == 'new':
        analyse_new(args)

    elif args['mode'] == 'wj':
        analyse_wj(args)

    elif args['mode'] == 'q':
        quick_rmsd_total(args)

    elif args['mode'] == 'slider':
        slide_ddg(args)

    elif args['mode'] == 's_by_s':
        side_by_side(args)

    elif args['mode'] == 'mutant_table':
        mutant_table(args)

    elif args['mode'] == 'erbb2_mutants':
        erbb2_mutants(args)

    elif args['mode'] == 'design_fnd':
        design_fnd_scatter(args)

    elif args['mode'] == 'all_fnds':
        draw_all_fnds(args)

    else:
        print('no mode')
    args['logger'].close()


def erbb2_mutants(args):
    """
    draw what rosetta thinks about assaf's mutants of ErbB2
    """
    score_dir = "/home/labs/fleishman/jonathaw/elazaridis/fold_and_dock/" + \
        "erbb2/mutations/all_results/"
    exp_table = "/home/labs/fleishman/jonathaw/elazaridis/fold_and_dock/" + \
        "erbb2/mutations/general_data/mut_table.txt"
    exp_df = parse_erbb2_exp_table(exp_table)
    wt_score_file = score_dir + "all_erbb2v4_wt_28Feb.score"
    wt_df = Rf.score_file2df(wt_score_file)
    wt_ddg = Rf.get_term_by_threshold(wt_df, 'score', 5, 'a_ddg', 'mean')

    exp_df['rosetta'] = np.nan
    # exp_df['rosetta_score'] = np.nan
    for sc_file in [a for a in os.listdir(score_dir)
                    if '.score' in a and 'wt' not in a]:
        df = Rf.score_file2df(score_dir+sc_file)
        ddg = Rf.get_term_by_threshold(df, 'score', 5, 'a_ddg', 'mean')
        # scr = Rf.get_term_by_threshold(df, 'score', 5, 'score', 'mean')
        name = sc_file.split('_')[2]
        # print(sc_file, name)
        wt = name[0]
        pos = int(name[1:-1])
        mut = name[-1]
        exp_df.set_value((exp_df['pos'] == pos) & (exp_df['wt'] == wt) &
                         (exp_df['mut'] == mut), 'rosetta', ddg-wt_ddg)
    print(exp_df)
    exp_df = exp_df.dropna()
    print(exp_df.to_string())
    plt.scatter(exp_df['rosetta'], exp_df['exp'])
    plt.ylabel('experimental ∆∆G')
    plt.xlabel('rosetta ∆∆G')
    plt.axhline(0)
    plt.axvline(0)
    for i, row in exp_df.iterrows():
        plt.annotate('%s%i%s' % (row['wt'], row['pos'], row['mut']),
                     (row['rosetta'], row['exp']))
    plt.show()


def parse_erbb2_exp_table(file_name: str) -> pd.DataFrame:
    result = {}
    with open(file_name, 'r') as fin:
        first = True
        for l in fin:
            s = l.split()
            if first:
                for i, a in enumerate(s):
                    result[i+1] = {'wt': a}
                first = False
            else:
                for i, a in enumerate(s[1:]):
                    result[i+1][s[0]] = float(a)
    # df = pd.DataFrame.from_dict(result, "index")
    new_df = pd.DataFrame(columns=["pos", "wt", "mut", "exp"])
    for pos, dct in result.items():
        for k, v in dct.items():
            if k != 'wt':
                t = {'pos': pos, 'wt': dct['wt'], 'mut': k, 'exp': v}
                new_df = new_df.append(t, ignore_index=True)
    return new_df


def mutant_table( args: dict ):
    """
    a function to find and display the correlation between ResSolv and
    MPFrameWork and experimental results from both Doung 2006 and Assaf
    """
    scores_dir = '/home/labs/fleishman/jonathaw/elazaridis/fold_and_dock/gpa/mutant_results/%s' % args['dir']
    mp_dir = '/home/labs/fleishman/jonathaw/elazaridis/fold_and_dock/gpa/mutant_results/mpframework_18Dec/'
    main_df = pd.read_csv("/home/labs/fleishman/jonathaw/elazaridis/" +
                          "fold_and_dock/gpa/mutant_results/" +
                          "experimental_results.tsv", sep='\s+')
    wt_beta_score_file = [a for a in os.listdir(scores_dir)
                          if 'wt' in a and '.score' in a][0]
    wt_beta_df = Rf.score_file2df(scores_dir + '/' + wt_beta_score_file)
    wt_beta_ddg = Rf.get_term_by_threshold(wt_beta_df, 'score', 5, 'a_ddg',
                                           'mean')
    wt_mp_df = Rf.score_file2df('%sall_gpav1_wt_mpframework_25Oct.score' % mp_dir)
    wt_mp_ddg = Rf.get_term_by_threshold(wt_mp_df, 'score', 5, 'a_ddg', 'mean')
    results = {'rs': {}, 'mp': {}}

    for sc_file in [a for a in os.listdir(scores_dir)+os.listdir(mp_dir)
                    if '.score' in a]:
        if 'mpframework' in sc_file:
            df = Rf.score_file2df('%s/%s' % (mp_dir, sc_file))
        else:
            df = Rf.score_file2df('%s/%s' % (scores_dir, sc_file))
        name = sc_file.split('_')[2]
        if '16Mar' in sc_file:
            name = '%s%i%s' % (name[0], int(name[1:-1])+72, name[-1])
        # if name[-1] == 'M': continue
        # threshold = np.percentile(df['score'].values, 5)
        min_ddg = Rf.get_term_by_threshold(df, 'score', 5, 'a_ddg', 'mean')
        if 'mpframework' in sc_file:
            results['mp'][name] = min_ddg
            main_df.set_value(main_df['name'] == name, 'mp', min_ddg-wt_mp_ddg)
        else:
            results['rs'][name] = min_ddg
            main_df.set_value(main_df['name'] == name, 'rs',
                              min_ddg-wt_beta_ddg)

    print(main_df)
    # main_df = main_df.dropna( how='any' )
    args['logger'].log(main_df)

    if args['all4']:
        fig = plt.figure(figsize=(10, 10), facecolor='w')
        i = 1
        for scfxn in ['rs', 'mp']:
            for exp in ['dstbl', 'Doung']:
                ax = plt.subplot(2, 2, i)
                model = linear_model.LinearRegression()
                model.fit(main_df[scfxn].to_frame(), main_df[exp].to_frame())
                line_x = np.linspace(main_df[scfxn].min(), main_df[scfxn].max())
                line_y = model.predict(line_x[:, np.newaxis])
                r2 = r2_score(main_df[exp].values,
                              model.predict(main_df[scfxn].to_frame()))
                plt.scatter(main_df[scfxn], main_df[exp])
                plt.plot(line_x, line_y)
                scfxn_name = 'ResSolv' if scfxn == 'rs' else 'MPFrameWork'
                exp_name = 'Doung 2006' if exp == 'Doung' else r'dsT$\beta$L'
                plt.title('%s Vs. %s' % (scfxn_name, exp_name))
                plt.text(0.8, 0.1, r'$R^2=%.2f$' % r2, fontsize=15,
                         horizontalalignment='center',
                         verticalalignment='center', transform=ax.transAxes)
                plt.axhline(0, color='k')
                plt.axvline(0, color='k')
                if i == 3:
                    plt.xlabel('Rosetta ∆∆G', fontsize=18)
                    plt.ylabel('Experimental ∆∆G', fontsize=18)
                i += 1
        plt.show()

    else:
        fig = plt.figure(facecolor='w')
        ax1 = plt.subplot(1, 2, 1)
        model = linear_model.LinearRegression()
        rs_df = main_df[['name', 'dstbl', 'rs']].dropna(how='any')

        model.fit(rs_df['rs'].to_frame(), rs_df['dstbl'].to_frame())
        line_x = np.linspace(rs_df['rs'].min(), rs_df['rs'].max())
        line_y = model.predict(line_x[:, np.newaxis])
        r2 = r2_score(rs_df['dstbl'].values,
                      model.predict(rs_df['rs'].to_frame()))
        plt.scatter(rs_df['rs'], rs_df['dstbl'])
        plt.plot(line_x, line_y)
        plt.title('%s Vs. %s' % ('ResSolv', r'dsT$\beta$L'))
        plt.text(0.8, 0.1, r'$R^2=%.2f$' % r2, fontsize=15,
                 horizontalalignment='center', verticalalignment='center',
                 transform=ax1.transAxes)
        plt.axhline(0, color='k')
        plt.axvline(0, color='k')
        plt.xlabel('Rosetta ∆∆G', fontsize=18)
        plt.ylabel(r'dsT$\beta$L experimental results', fontsize=18)
        for x, y, n in zip(rs_df['rs'], rs_df['dstbl'], rs_df['name']):
            ax1.annotate(n, (x, y))

        ax2 = plt.subplot(1, 2, 2)
        model = linear_model.LinearRegression()
        mp_df = main_df[['name', 'dstbl', 'mp']].dropna(how='any')
        model.fit(mp_df['mp'].to_frame(), mp_df['dstbl'].to_frame())
        line_x = np.linspace(mp_df['mp'].min(), mp_df['mp'].max())
        line_y = model.predict(line_x[:, np.newaxis])
        r2 = r2_score(mp_df['dstbl'].values,
                      model.predict(mp_df['mp'].to_frame()))
        plt.scatter(mp_df['mp'], mp_df['dstbl'])
        plt.plot(line_x, line_y)
        plt.title('%s Vs. %s' % ('MPFrameWork', r'dsT$\beta$L'))
        plt.text(0.8, 0.1, r'$R^2=%.2f$' % r2, fontsize=15,
                 horizontalalignment='center', verticalalignment='center',
                 transform=ax2.transAxes)
        plt.axhline(0, color='k')
        plt.axvline(0, color='k')
        # plt.xlabel( 'Rosetta ∆∆G', fonctsize=18 )
        # plt.ylabel( r'dsT$\beta$L experimental results', fonctsize=18 )
        plt.show()
        plt.savefig('%s/dsTbL_alone.pdf' % scores_dir)


def analyse_wj(args):
    if not os.path.isfile(args['obj_wj']):
        print('creating wj df')
        df = get_fold_and_dock_df(args['pc_wj'], args['score_wj'])
        pickle.dump(df, open(args['obj_wj'], 'wb'))
    else:
        print('reading wj df')
        df = pickle.load(open(args['obj_wj'], 'rb'))
    print('finished getting df')
    print(df)
    draw_scatter(df, title=args['obj_wj'])  # [df['angle'] <= 90.0][df['a_sasa'] > 800], title='wt_jiggle')
    return df


def analyse_old(args) -> pd.DataFrame:
    if not os.path.isfile(args['obj_old']):
        print('creating old df')
        df = get_fold_and_dock_df(args['pc_old'], args['score_old'])
        pickle.dump(df, open(args['obj_old'], 'wb'))
    else:
        print('reading old df')
        df = pickle.load(open(args['obj_old'], 'rb'))
    print('finished getting df')

    draw_scatter(df[df['angle'] <= 90.0][df['a_sasa'] > 800], title='Full')
    return df


def analyse_new(args) -> pd.DataFrame:
    if not os.path.isfile(args['obj_new']):
        print('creating new df')
        df = get_fold_and_dock_df(args['pc_new'], args['score_new'])
        pickle.dump(df, open(args['obj_new'], 'wb'))
    else:
        print('reading new df')
        df = pickle.load(open(args['obj_new'], 'rb'))
    print('finished getting df')
    draw_scatter(df[df['angle'] <= 90.0][df['a_sasa'] > 800], 'ResSolv')
    return df


def get_fold_and_dock_df(pc_file: str, score_file: str) -> pd.DataFrame:
    df = pd.DataFrame(columns=['name', 'rmsd', 'angle', 'filter'])
    pc = parse_pc_all(pc_file)
    missed = []
    for name, flt_dict in generate_score(score_file, filters):
        if name not in pc.keys():
            missed.append(name)
            continue
        res = {'name': name, 'rmsd': pc[name]['rmsd'], 'angle': pc[name]['angle']}
        for flt, r in flt_dict.items():
            res[flt] = r
        df = df.append(res, ignore_index=True)
    print('found %i rows, and missed %i' % (len(df), len(missed)))
    return df


def draw_scatter(df, title=None) -> None:
    print(df)
    plt.scatter(df['rmsd'].values, df['total'].values, label=df['name'])
    for name, rmsd, filter in zip(df['name'], df['rmsd'], df['total']):
        plt.annotate(name, xy=(rmsd, filter), xytext=(rmsd, filter))
    # plt.xlim([0, 20])
    # plt.ylim([-300, 100])
    plt.title(title)
    # print_best_scores(df)
    plt.show()


def print_best_scores(df_, filter='total', percentile=0.1):
    perc_threshold = np.nanpercentile(df_[filter].values, percentile)
    print('threshold for %.2f percentile is %.2f' % (percentile, perc_threshold))
    for name, rmsd, filter in zip(df_['description'], df_['rmsd_calc'], df_[filter]):
        if filter <= perc_threshold:
            print('%s %.2f %.2f' % (name, rmsd, filter))


def parse_pc_all(file_name) -> dict:
    """
    parse a file with all pymol calc data with order "pc_NAME.txt ANGLE RMSD"
    """
    result = {}
    for l in open(file_name, 'r'):
        s = l.split()
        try:
            result[s[0].split('pc_')[1].split('.txt')[0].split('.pdb')[0]] = {'angle': float(s[1]), 'rmsd': float(s[2])}
        except:
            result[s[0]] = {'angle': float(s[1]), 'rmsd': float(s[2])}
            pass
    return result


def generate_score(file_name, filters: list):
    result = {}
    for i, l in enumerate(open(file_name, 'r')):
        s = l.split()
        if i == 0:
            fields = {a: i for i, a in enumerate(s)}
            # if 'score' not in fields.keys():
            #     fields['total'] = fields['total_score']
        else:
            # try:
            # print({filter: float(s[fields[filter]]) for filter in filters})
            yield s[fields['description']], {filter: float(s[fields[filter]]) for filter in filters}


def get_rmsds_from_table(pymol_calc_file:str) -> pd.DataFrame:
    t = pd.read_csv(pymol_calc_file, header=None, names=['description', 'pc_rmsd'], sep='\s', engine='python')
    return t


def slide_ddg(args):
    global scat, ax, ddg_slider, sc_df, fig
    sc_df = Rf.score_file2df(args['sc'], args['names'])
    pc_df = get_rmsds_from_table(args['pc'])
    a = sc_df.merge(pc_df, on='description')
    sc_df = a.copy()

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.25, bottom=0.25)

    scat = ax.scatter(sc_df['pc_rmsd'].values, sc_df['score'].values)
    # sc_df.plot(kind='scatter', x='pc_rmsd', y='score')

    slider_ax = plt.axes([0.25, 0.1, 0.65, 0.03])
    ddg_slider = Slider(slider_ax, 'ddG', np.min(sc_df['a_ddg'].values), +5, 0) #np.max(sc_df['a_ddg'].values), valinit=-8)

    ddg_slider.on_changed(update)
    plt.show()


def update(val):
    global scat#, ax
    threshold = ddg_slider.val
    new_df = sc_df[ sc_df['a_ddg'] < threshold ]
    scat.remove()
    scat = ax.scatter(new_df['pc_rmsd'].values, new_df['score'].values)
    ax.set_xlim(0, np.max(new_df['pc_rmsd'].values)+1)
    min_ddg, max_ddg = np.min(new_df['score'].values), np.max(new_df['score'].values)
    ax.set_ylim(min_ddg-1, max_ddg+1)
    fig.canvas.draw_idle()


def side_by_side(args):
    mp_sc = Rf.score_file2df(args['mp_sc'])
    mp_pc = get_rmsds_from_table(args['mp_pc'])
    a = mp_sc.merge(mp_pc, on='description')
    mp_sc = a.copy()
    rs_sc = Rf.score_file2df(args['rs_sc'])
    rs_pc = get_rmsds_from_table(args['rs_pc'])
    b = rs_sc.merge(rs_pc, on='description')
    rs_sc = b.copy()

    names_dict = {}
    for i, d in enumerate(mp_sc['description'].values):
        names_dict[d] = i if 'MP' in d else i + 100
    for i, d in enumerate(rs_sc['description'].values):
        if d not in names_dict.keys():
            names_dict[d] = i if 'MP' in d else i + 100
    for k, v in names_dict.items():
        print(v, k)

    axmp = plt.subplot(121)
    axmp.scatter(mp_sc['pc_rmsd'].values, mp_sc['tot_mp_fa'].values, label=mp_sc['description'].values)
    for x, y, d in zip(mp_sc['pc_rmsd'].values, mp_sc['tot_mp_fa'], mp_sc['description']):
        axmp.annotate(names_dict[d], (x, y))
    # axmp.title('MP')

    axrs = plt.subplot(122)
    axrs.scatter(rs_sc['pc_rmsd'].values, rs_sc['tot_rs_fa'].values, label=rs_sc['description'].values)
    for x, y, d in zip(rs_sc['pc_rmsd'].values, rs_sc['tot_rs_fa'], rs_sc['description']):
        axrs.annotate(names_dict[d], (x, y))
    # axrs.title('RS')

    plt.show()


def draw_all_fnds(args):
    dir = '/home/labs/fleishman/jonathaw/elazaridis/design/polyA_13Nov/' + \
        'chosen_designs_fnd/seq_diver_13Mar'
    current = [l.split()[0] for l in open('%s/current.txt' % dir, 'r')]
    # fnd_dirs = [a for a in os.listdir(dir) if 'poly' in a]
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                        wspace=0.15, hspace=0.45)
    for i, fnd_dir in enumerate(current):
        print(fnd_dir)
        os.chdir('%s/fnd_13.3_beta' % fnd_dir)
        score_file = glob.glob('*.score')[0]
        temp_args = {'sc': score_file, 'pc': glob.glob('all_pc_*.txt')[0],
                     'percent': 5,
                     'original_sc': '../../all_chosen_scores/score.sc',
                     'name': fnd_dir, 'mode': 'all_fnds', 'names': None,
                     'logger': args['logger']}
        sc_df = design_fnd_scatter(temp_args)

        plt.subplot(6, 4, i+1)
        plt.scatter(sc_df['pc_rmsd'].values, sc_df['score'].values, marker='o',
                    c=sc_df['a_span_topo'].values, picker=True,
                    cmap=plt.cm.coolwarm)

        min_energy = np.nanmin(list(sc_df['score'].values))
        max_energy = np.nanmax(list(sc_df['score'].values))
        plt.ylim([min_energy - 1, max_energy + 1])
        plt.xlim([0, 15])
        plt.title(fnd_dir.split('poly')[1])

        os.chdir(dir)

    plt.show()


def design_fnd_scatter(args):
    sc_df = Rf.score_file2df(args['sc'], args['names'])
    args['logger'].log('found %i structs in sc_df' % len(sc_df))
    pc_df = get_rmsds_from_table(args['pc'])
    args['logger'].log('found %i structs in pc' % len(pc_df))
    a = sc_df.merge(pc_df, on='description')
    sc_df = a.copy()
    sc_df = sc_df[sc_df['a_tms_span'] > 0.5]

    threshold = np.percentile(sc_df['score'], args['percent'])
    sc_df = sc_df[sc_df['score'] < threshold]

    original_df = Rf.score_file2df(args['original_sc'])
    for d in original_df['description']:
        if args['name'] in d:
            row_name = d
    original_row = original_df[original_df['description'] == row_name]
    term_dict = {'total_score': {'ou': 'under'},
                 'a_sasa': {'ou': 'over'},
                 'a_pack': {'ou': 'over'},
                 # 'a_shape': {'ou': 'over'},
                 'a_res_solv': {'ou': 'under'},
                 'a_ddg': {'ou': 'under'},
                 'a_span_topo': {'ou': 'over'}}
    for term in term_dict.keys():
        if term == 'a_res_solv':
            term_dict[term]['threshold'] = 0.5 * original_row[term].values[0]
        else:
            term_dict[term]['threshold'] = 0.8 * original_row[term].values[0]

    sc_df, fail_msg = Rf.remove_failed_dict(sc_df, term_dict)
    for k, v in fail_msg.items():
        print(v)

    if args['mode'] == 'all_fnds':
        return sc_df

    fig, ax = plt.subplots()
    ax.scatter(sc_df['pc_rmsd'].values, sc_df['score'].values, marker='o',
               c=sc_df['a_span_topo'].values, picker=True, cmap=plt.cm.coolwarm)

    min_energy = np.nanmin(list(sc_df['score'].values))
    max_energy = np.nanmax(list(sc_df['score'].values))
    plt.ylim([min_energy - 1, max_energy + 1])
    plt.xlim([0, 15])
    plt.title(args['name'])

    z_score, rmsd_threshold = rf.get_z_score_by_rmsd_percent(sc_df)
    plt.text(0.75, 0.2, "zscore=%.2f" % z_score, transform=ax.transaxes)
    plt.axvline(rmsd_threshold)

    point_label_cols = list(set(args['terms'] + ['description', 'a_sasa',
                                                 'a_res_solv', 'a_pack',
                                                 'a_span_topo', 'a_ddg',
                                                 'fa_elec']))
    pl = PointLabel(sc_df, ax, fig, 'pc_rmsd', 'score', point_label_cols,
                    args['logger']) # a_shape ???
    fig.canvas.mpl_connect('pick_event', pl.onpick)

    if args['show'] == 'show':
        plt.show()
    else:
        plt.savefig('%s.png' % args['name'])


def quick_rmsd_total(args):

    y_axis_term = 'score'

    sc_df = Rf.score_file2df(args['sc'], args['names'])
    args['logger'].log('found %i structs in sc_df' % len(sc_df))
    pc_df = get_rmsds_from_table(args['pc'])
    args['logger'].log('found %i structs in pc' % len(pc_df))
    a = sc_df.merge(pc_df, on='description')
    sc_df = a.copy()

    # if 'a_hha' in sc_df.columns:
        # sc_df['angle'] = sc_df['a_hha'] > 0

    args['logger'].log('left with %i in merged df' % len(sc_df))

    args['logger'].log('examining %s with span_topo threshold %f' % (args['sc'], args['span_threshold']))
    fig, ax = plt.subplots()

    if args['best']:
        sc_df = sc_df[sc_df['a_tms_span_fa'] > 0.5 ]
        threshold = np.percentile(sc_df[y_axis_term], args['percent'])
        sc_df = sc_df[ sc_df[y_axis_term] < threshold ]
        sc_df = sc_df[ sc_df['a_span_topo'] >= 0.99 ]
        sc_df_pass = Rf.get_best_of_best(sc_df, args['terms'], args['threshold'])
        sc_df_fail = sc_df[ ~sc_df['description'].isin( sc_df_pass['description'] ) ]
        args['logger'].log('%i models returned from BEST' % len(sc_df_pass))
    else:
        args['logger'].log('total of %i models in score' % len(sc_df))
        sc_df = sc_df[sc_df['a_tms_span_fa'] > 0.5]
        args['logger'].log('%i models pass tms_span' % len(sc_df))
        threshold = np.percentile(sc_df[y_axis_term], args['percent'])
        sc_df = sc_df[ sc_df[y_axis_term] < threshold ]
        args['logger'].log('for percent %.2f found threshold to be %.2f and %i strucutres pass it' % (args['percent'], threshold, len(sc_df)))
        sc_df = sc_df[sc_df['a_shape'] >= 0.6]
        sc_df = sc_df[sc_df['a_sasa'] > 700]
        args['logger'].log('%i passed sasa 600' % len(sc_df))
        sc_df = sc_df[sc_df['a_ddg'] < -5]
        args['logger'].log('%i passed ddg' % len(sc_df))
        # sc_df = sc_df[sc_df['a_pack'] > 0.6]
        sc_df = sc_df[sc_df['a_unsat'] < 1]
        args['logger'].log('%i passed unsat' % len(sc_df))
        sc_df['pass'] = sc_df['a_span_topo'] > args['span_threshold']
        sc_df = sc_df[sc_df['a_res_solv'] < -10]
        args['logger'].log('%i passed res_solv -10' % len(sc_df))

        sc_df_pass = sc_df[sc_df['a_span_topo'] > args['span_threshold']]
        args['logger'].log('%i models passed span_topo threshold' % len(sc_df_pass))
        sc_df_fail = sc_df[sc_df['a_span_topo'] <= args['span_threshold']]
        args['logger'].log('%i models failed span_topo threshold' % len(sc_df_fail))

    # ax.scatter(sc_df_fail['rmsd_calc'].values, sc_df_fail['score'].values, color='r', marker='.')

    x_array = np.ndarray(buffer=sc_df_pass['pc_rmsd'].values, shape=(len(sc_df),))
    y_array = np.ndarray(buffer=sc_df_pass[y_axis_term].values, shape=(len(sc_df)))
    if 'a_hha' in sc_df.columns:
        ax.scatter(sc_df_pass['pc_rmsd'].values, sc_df_pass[y_axis_term].values, marker='o',
                c=sc_df_pass['a_hha'].values, picker=True, cmap=plt.cm.coolwarm)
    else:
        ax.scatter(sc_df_pass['pc_rmsd'].values, sc_df_pass[y_axis_term].values, marker='o',
                c=sc_df_pass['a_span_topo'].values, picker=True, cmap=plt.cm.coolwarm)

    # min_energy = np.nanmin(list(sc_df_pass['score'].values)+list(sc_df_fail['score'].values))
    min_energy = np.nanmin(list(sc_df_pass[y_axis_term].values))
    max_energy = np.nanmax(list(sc_df_pass[y_axis_term].values))
    plt.ylim([min_energy - 1, max_energy + 1])
    plt.xlim([0, 15])
    plt.title(args['sc']+'_pass')

    z_score, rmsd_threshold = Rf.get_z_score_by_rmsd_percent(sc_df_pass)
    plt.text(0.75, 0.2, "Zscore=%.2f" % z_score, transform=ax.transAxes)
    plt.axvline(rmsd_threshold)
    # if 'a_hha' in sc_df.columns:
        # ax.scatter(sc_df_fail['pc_rmsd'].values, sc_df_fail[y_axis_term].values, marker='x',
                # c=sc_df_fail['a_hha'].values, picker=True, cmap=plt.cm.coolwarm, s=5, alpha=90)#, markersize=200)
    # else:
        # ax.scatter(sc_df_fail['pc_rmsd'].values, sc_df_fail[y_axis_term].values, marker='x',
                # c=sc_df_fail['a_span_topo'].values, picker=True, cmap=plt.cm.coolwarm, s=5, alpha=90)#, markersize=200)

    # af = PrintLabel(sc_df_pass, 'rmsd_calc', 'score', ['description', 'pass'])
    # fig.canvas.mpl_connect('button_press_event', af)
    point_label_cols = list(set(args['terms'] + ['description', 'a_sasa', 'a_res_solv', 'a_pack', 'a_span_topo', 'a_ddg', 'fa_elec']))
    pl = PointLabel(sc_df_pass, ax, fig, 'pc_rmsd', y_axis_term, point_label_cols,
                    args['logger']) # a_shape ???
    fig.canvas.mpl_connect('pick_event', pl.onpick)
    plt.xlabel('RMSD')
    plt.ylabel(y_axis_term)
    if args['show'] == 'show':
        # fig.canvas.mpl_connect('pick_event', on_pick3)
        # cursor = FollowDotCursor(ax, sc_df_pass['pc_rmsd'], sc_df_pass[y_axis_term])
        plt.show()
    else:
        plt.savefig('%s.png' % args['sc'].split('.score')[0])


def on_pick3(event):
    ind = event.ind
    print('pick', ind, np.take(x_array, ind), np.take(y_array, ind))


class PointLabel:
    def __init__(self, df: pd.DataFrame, ax, fig, x_axis: str, y_axis: str,
                 labels: list, file_=None):
        self.df = df.copy()
        self.axis = ax
        self.fig = fig
        self.x_axis = x_axis
        self.y_axis = y_axis
        self.labels = labels
        if 'score' in self.labels:
            self.labels.remove('score')
        self.has_written_title = False
        if file_ is not None:
            self.file_handler = file_
        else:
            self.file_handler = Logger('point_label.log')
        self.label_order = [ self.x_axis, self.y_axis ] + [label for label in self.labels if label != 'description']

    def onpick(self, event):
        """
        inspired by http://matplotlib.org/examples/event_handling/pick_event_demo.html
        """
        ind = event.ind[0]
        row = self.df.iloc[ind]
        if not self.has_written_title:
            # self.file_handler.log('%s %s %s description' %
                                  # (self.x_axis, self.y_axis, ''.join('{:>10}'.format( label ) for label in self.labels if label != 'description')), skip_stamp=True)
            self.file_handler.log(''.join('{:<10}'.format( label ) for label in self.label_order+['description']))
            self.has_written_title = True
        ord_row = [ row[self.x_axis], row[self.y_axis] ] + [row[label] for label in self.labels if label != 'description']
        # self.file_handler.log('picker %.2f\t%.2f\t%s %s' %
              # (row[self.x_axis], row[self.y_axis],
               # '\t'.join("%.2f" % row[label] for label in self.labels if label != 'description'),
               # row['description']), skip_stamp=True)
        # print( ord_row )
        self.file_handler.log('%s %s' % ( ''.join(['{:<10.2f}'.format( r ) for r in ord_row]), row['description'] ), skip_stamp=1)


class PrintLabel(object):
    def __init__(self, df: pd.DataFrame, x_axis: str, y_axis: str, labels: list):
        self.df = df
        self.x_axis = x_axis
        self.y_axis = y_axis
        self.label = labels

    def __call__(self, event):
        clk_x = event.xdata
        clk_y = event.ydata
        self.closest_row(clk_x, clk_y)

    def closest_row_(self, x, y):
        x_df = self.df.copy()
        x_diff = x_df[self.x_axis].apply(lambda z: abs(x - z))
        x_diff.sort()
        print('hhh', x_diff[:100])
        x_inds = x_diff.index[:10]
        # print('for x', x, 'found these', x_df.ix[x_diff.index[:5]]['rmsd_calc'])

        y_df = self.df.copy()
        y_diff = y_df[self.y_axis].apply(lambda z: abs(y - z))
        y_diff.sort()
        y_inds = y_diff.index[:10]
        # print('for y', y, 'found these', y_diff.index[:5])

        both = set(list(x_inds)) & set(list(y_inds))
        print('closest to %.2f, %.2f found these %i points' % (x, y, len(both)))
        for a in both:
            print(self.df.ix[int(a)][[self.x_axis, self.y_axis] + self.label])

    def closest_row(self, x, y):
        print('self.df', self.df)
        x_df = self.df.copy()
        x_df = x_df.iloc([self.x_axis] - x).abs().argsort()
        y_df = self.df.copy()
        y_df = y_df.iloc[(y_df[self.y_axis] - y).abs().argsort()]
        print('x', x_df.index)
        print('y', y_df.index)


if __name__ == '__main__':
    main()
