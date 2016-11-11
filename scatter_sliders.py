#!/usr/bin/env python3.5
import os
import time
import argparse
import pickle
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import RosettaFilter as Rf
from Logger import Logger
import pandas as pd
import numpy as np
from AnnotateFinder import AnnoteFinder

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
    parser.add_argument('-pc')
    parser.add_argument('-span_threshold', default=0.3, type=float)
    parser.add_argument('-names', default=None, help='file list of names to show')
    parser.add_argument('-log_path', default='./', help='path to place log file')

    args = vars(parser.parse_args())
    args['logger'] = Logger('logeer_%s.log' % time.strftime("%d.%0-m"), args['log_path'])

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

    else:
        print('no mode')
    args['logger'].close()


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
            # except:
            #     pass


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


def quick_rmsd_total(args):

    y_axis_term = 'score'

    sc_df = Rf.score_file2df(args['sc'], args['names'])
    pc_df = get_rmsds_from_table(args['pc'])
    a = sc_df.merge(pc_df, on='description')
    sc_df = a.copy()

    args['logger'].log('examining %s with span_topo threshold %f' % (args['sc'], args['span_threshold']))
    fig, ax = plt.subplots()

    args['logger'].log('total of %i models in score' % len(sc_df))
    sc_df = sc_df[sc_df['a_tms_span_fa'] > 0.5]
    args['logger'].log('%i models pass tms_span' % len(sc_df))
    # sc_df = sc_df[sc_df['a_shape'] >= 0.6]
    # sc_df = sc_df[sc_df['a_sasa'] > 900]
    sc_df = sc_df[sc_df['a_ddg'] < -6]
    # sc_df = sc_df[sc_df['a_pack'] > 0.6]
    sc_df['pass'] = sc_df['a_span_topo'] > args['span_threshold']

    sc_df_pass = sc_df[sc_df['a_span_topo'] > args['span_threshold']]
    args['logger'].log('%i models passed span_topo threshold' % len(sc_df_pass))
    sc_df_fail = sc_df[sc_df['a_span_topo'] <= args['span_threshold']]
    args['logger'].log('%i models failed span_topo threshold' % len(sc_df_fail))

    # ax.scatter(sc_df_fail['rmsd_calc'].values, sc_df_fail['score'].values, color='r', marker='.')
    ax.scatter(sc_df_pass['pc_rmsd'].values, sc_df_pass[y_axis_term].values, marker='o',
               c=sc_df_pass['a_span_topo'].values, picker=True, cmap=plt.cm.coolwarm)

    # min_energy = np.nanmin(list(sc_df_pass['score'].values)+list(sc_df_fail['score'].values))
    min_energy = np.nanmin(list(sc_df_pass[y_axis_term].values))
    max_energy = np.nanmax(list(sc_df_pass[y_axis_term].values))
    plt.ylim([min_energy - 1, max_energy + 1])
    plt.xlim([0, 30])
    plt.title(args['sc']+'_pass')

    ax.scatter(sc_df_fail['pc_rmsd'].values, sc_df_fail[y_axis_term].values, marker='x',
               c=sc_df_fail['a_span_topo'].values, picker=True, cmap=plt.cm.coolwarm)#, markersize=200)

    # af = PrintLabel(sc_df_pass, 'rmsd_calc', 'score', ['description', 'pass'])
    # fig.canvas.mpl_connect('button_press_event', af)
    pl = PointLabel(sc_df_pass, ax, fig, 'pc_rmsd', y_axis_term,
                    ['description', 'a_sasa', 'a_res_solv',
                     'a_pack', 'a_span_topo', 'a_ddg'], args['logger']) # a_shape ???
    fig.canvas.mpl_connect('pick_event', pl.onpick)
    # print('for pass')
    # print_best_scores(sc_df_pass, 'score', percentile=0.05)
    # print('for. fail')
    # print_best_scores(sc_df_fail, 'score', percentile=0.05)
    plt.xlabel('RMSD')
    plt.ylabel(y_axis_term)
    plt.show()

    # Rf.df2boxplots(sc_df_pass)


class PointLabel:
    def __init__(self, df: pd.DataFrame, ax, fig, x_axis: str, y_axis: str,
                 labels: list, file_=None):
        self.df = df.copy()
        self.axis = ax
        self.fig = fig
        self.x_axis = x_axis
        self.y_axis = y_axis
        self.labels = labels
        self.has_written_title = False
        if file_ is not None:
            self.file_handler = file_
        else:
            self.file_handler = Logger('point_label.log')

    def onpick(self, event):
        """
        inspired by http://matplotlib.org/examples/event_handling/pick_event_demo.html
        """
        ind = event.ind[0]
        row = self.df.iloc[ind]
        if not self.has_written_title:
            self.file_handler.log('picker %s %s %s description' %
                  (self.x_axis, self.y_axis, '\t'.join(label for label in self.labels if label != 'description')), skip_stamp=True)
            self.has_written_title = True
        self.file_handler.log('picker %.2f\t%.2f\t%s %s' %
              (row[self.x_axis], row[self.y_axis],
               '\t'.join("%.2f" % row[label] for label in self.labels if label != 'description'),
               row['description']), skip_stamp=True)
        # self.axis.text(row[self.x_axis], row[self.y_axis], row['description'], transform=self.axis.transAxes,
        #                bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10})
        # self.fig.canvas.draw()


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
