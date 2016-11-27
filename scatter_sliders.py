#!/usr/bin/env python3.5
import os
import time
import argparse
import pickle
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import RosettaFilter as Rf
from Logger import Logger
from ELazaridis_fold_and_dock_analysis import PointLabel
import pandas as pd
from prettypandas import PrettyPandas
from tabulate import tabulate
import numpy as np
from AnnotateFinder import AnnoteFinder

filters = ['total', 'a_sasa']


def main():
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', default='q')
    parser.add_argument('-sc')
    parser.add_argument('-pc', default=None)
    parser.add_argument('-names', default=None, help='file list of names to show')
    parser.add_argument('-log_path', default='./', help='path to place log file')
    parser.add_argument('-terms', nargs='+', default=['score', 'a_ddg', 'a_sasa'])
    parser.add_argument('-percent', default=100)
    parser.add_argument('-x', default='pc_rmsd')
    parser.add_argument('-y', default='score')

    args = vars(parser.parse_args())
    args['logger'] = Logger('logeer_%s.log' % time.strftime("%d.%0-m"), args['log_path'])

    slide_ddg(args)

    args['logger'].close()


def get_rmsds_from_table(pymol_calc_file:str) -> pd.DataFrame:
    t = pd.read_csv(pymol_calc_file, header=None, names=['description', 'pc_rmsd'], sep='\s', engine='python')
    return t


def slide_ddg(args):
    global new_df, radio, color_by, picked
    global scat, ax, sliders, sc_df, fig, cm, cbar
    sc_df = Rf.score_file2df(args['sc'], args['names'])
    args['logger'].log('score file has %i entries' % len(sc_df))
    if args['pc'] is not None:
        pc_df = get_rmsds_from_table(args['pc'])
        args['logger'].log('pc file had %i entries' % len(pc_df))
        a = sc_df.merge(pc_df, on='description')
        args['logger'].log('combined there are %i entries' % len(a))
        sc_df = a.copy()

    if args['percent'] != 100:
        threshold = np.percentile(sc_df[args['y']], args['percent'])
        sc_df = sc_df[ sc_df[args['y']] < threshold ]

    color_by = args['y']
    picked = False

    new_df = sc_df.copy()
    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.25, bottom=0.25)

    cm = plt.cm.get_cmap('RdYlBu')

    scat = ax.scatter(sc_df[args['x']].values, sc_df[args['y']].values, s=40, cmap=cm, c=sc_df[color_by], picker=True)
    cbar = plt.colorbar(scat)
    sliders = {}
    for i, term in enumerate(args['terms']):
        slider_ax = plt.axes([0.25, 0.01+i*0.035, 0.65, 0.03])
        sliders[term] = Slider(slider_ax, term, np.min(sc_df[term].values), np.max(sc_df[term].values), 0)
        sliders[term].on_changed(update)

    ax.set_xlim(np.min(new_df[args['x']].values)-1, np.max(new_df[args['x']].values)+1)
    ax.set_ylim(np.min(new_df[args['y']].values)-1, np.max(new_df[args['y']].values)+1)

    ax.set_xlabel(args['x'])
    ax.set_ylabel(args['y'])

    resetax = plt.axes([0.025, 0.7, 0.15, 0.15]) #[0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', color='lightgoldenrodyellow', hovercolor='0.975')
    button.on_clicked(reset)

    printax = plt.axes([0.025, 0.3, 0.15, 0.15])
    printbutton = Button(printax, 'Print', color='green', hovercolor='red')
    printbutton.on_clicked(print_table)

    logax = plt.axes([0.025, 0.1, 0.15, 0.15])
    logbutton = Button(logax, 'log table', color='blue', hovercolor='red')
    logbutton.on_clicked(log_table)

    rax = plt.axes([0.025, 0.5, 0.15, 0.15], axisbg='white')
    radio = RadioButtons(rax, args['terms'], active=0)
    radio.on_clicked(colorfunc)

    # cbar = plt.colorbar(scat)
    pl = PointLabel(new_df, ax, fig, args['x'], args['y'], ['description', 'a_sasa', 'a_res_solv', 'a_pack', 'a_span_topo', 'a_ddg', 'fa_elec'], args['logger'])
    fig.canvas.mpl_connect('pick_event', pl.onpick)

    plt.show()


def log_table(event):
    global new_df, picked
    if picked:
        picked = False
        return
    with open('table_log_%s.score' % time.strftime("%d.%0-m"), 'w+') as fout:
        fout.write(new_df.to_string())


def print_table(event):
    global new_df, picked
    if picked:
        picked = False
        return
    print(new_df.to_string())
    print('total rows: %i' % len(new_df))


def colorfunc(term):
    global new_df, color_by, scat, cbar, picked
    if picked:
        picked = False
        return
    color_by = term
    scat = ax.scatter(new_df[args['x']].values, new_df[args['y']].values, c=new_df[term], s=40, cmap=cm, picker=True)
    cbar.set_clim([ np.min(new_df[term]), np.max(new_df[term]) ])
    scat.set_clim([ np.min(new_df[term]), np.max(new_df[term]) ])
    # plt.colorbar(scat)
    fig.canvas.draw_idle()


def update(val):
    global scat, sliders, args, new_df, ax, color_by, picked
    if picked:
        picked = False
        return
    thresholds = {}
    new_df = sc_df.copy()
    for term in args['terms']:
        thresholds[term] = sliders[term].val
        if 'sasa' in term or 'pack' in term or 'span' in term:
            new_df = new_df[ new_df[term] > thresholds[term]]
        else:
            new_df = new_df[ new_df[term] < thresholds[term] ]
    ax.cla()
    ax.scatter(sc_df[args['x']].values, sc_df[args['y']].values, color='r', s=10, alpha=0.5, picker=True)
    ax.scatter(new_df[args['x']].values, new_df[args['y']].values, s=40, cmap=cm, c=new_df[color_by], picker=True)
    ax.set_xlim(np.min(new_df[args['x']].values)-1, np.max(new_df[args['x']].values)+1)
    ax.set_ylim(np.min(new_df[args['y']].values)-1, np.max(new_df[args['y']].values)+1)
    ax.set_xlabel(args['x'])
    ax.set_ylabel(args['y'])
    fig.canvas.draw_idle()


def reset(event):
    global sliders, scat, new_df, sc_df, ax, picked
    if picked:
        picked = False
        return
    for slider in sliders:
        sliders[slider].reset()
    new_df = sc_df.copy()
    # scat.remove()
    ax.cla()
    # ax.scatter(sc_df[args['x']].values, sc_df[args['y']].values, color='r', s=10, alpha=0.5)
    ax.scatter(sc_df[args['x']].values, sc_df[args['y']].values, s=40, cmap=cm, c=sc_df[args['y']], picker=True)
    # ax.set_xlim(np.min(sc_df[args['x']].values)-1, np.max(sc_df[args['x']].values)+1)
    # min_ddg, max_ddg = np.min(sc_df[args['y']].values), np.max(sc_df[args['y']].values)
    # ax.set_ylim(min_ddg-1, max_ddg+1)
    ax.set_xlim(np.min(sc_df[args['x']].values)-1, np.max(sc_df[args['x']].values)+1)
    ax.set_ylim(np.min(sc_df[args['y']].values)-1, np.max(sc_df[args['y']].values)+1)
    ax.set_xlabel(args['x'])
    ax.set_ylabel(args['y'])
    fig.canvas.draw_idle()


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
        global picked
        picked = True
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
        picked = False


if __name__ == '__main__':
    main()
