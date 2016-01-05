#!/usr/bin/env python3.5
"""
a script to analyse a DoCohModeller run results
"""
from __future__ import unicode_literals
from RosettaFilter import Filter, RunFilters, score2dict
import numpy as np
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
import matplotlib
from pandas import DataFrame, Series
from collections import OrderedDict
import os
import re
import random

matplotlib.rc_context(fname="/home/labs/fleishman/jonathaw/.matplotlib/publishable_matplotlibrc")

cohesin_names = ['ac_a5', 'ac_b1', 'ac_c3', 'ac_d1', 'ac_d3', 'af_75', 'af_76', 'bc_a11', 'bc_a5', 'bc_b3', 'ca_a1',
                 'ca_a2', 'ca_a5', 'cc_a1', 'cc_a8', 'cc_ox', 'ct_2p2', 'ct_a1', 'ct_a2', 'ct_a3', 'ct_oa', 'ct_ob1',
                 'ct_ob4', 'ct_s1', 'rf_a3', 'rf_b1', 'rf_b6', 'rf_c', 'rf_e', '4ums', '1ANU', '1AOH', '1G1K', '1QZN',
                 '1TYJ', '2BM3', '2JH2', '2VO8', '2W1N', '2W5F', '2XBT', '2XDH', '2ZF9', '3BWZ', '3FNK', '3L8Q', '4JO5',
                 '1OHZ', '2B59', '2CCL', '2OZN', '2VN5', '2VN6', '2Y3N', '3KCP', '3UL4', '4DH2', '4FL4', '4IU2', '4UYP',
                 '4UYQ', '4FL5', '5NEW']

dockerin_names = ['ct_48b', 'ct_11b', 'ct_cipa', 'ctcipadoc', 'cc_5a', 'rf_44b', 'rf_scaa', 'rf_scab', 'ca_9e',
                  'ac_scaa', 'ac_scab', 'ac_9b', 'bc_48a', 'bc_scaa', 'af_doc', '4ums', '1OHZ', '2B59', '2CCL', '2OZN',
                  '2VN5', '2VN6', '2Y3N', '3KCP', '3UL4', '4DH2', '4FL4', '4IU2', '4UYP', '4UYQ', '4FL5', '5NEW']


def monomer_data():
    """
    :return: a dictionary of thresholds for the different monomers
    """
    cohs = {'1ohz': {'score': -276.802000, 'packstat': 0.676000},
            '2b59': {'score': -314.548000, 'packstat': 0.656000},
            '2ccl': {'score': -281.684000, 'packstat': 0.686000},
            '2ozn': {'score': -259.101000, 'packstat': 0.698000},
            '2vn5': {'score': -280.649000, 'packstat': 0.740000},
            '2vn6': {'score': -279.445000, 'packstat': 0.698000},
            '2y3n': {'score': -319.608000, 'packstat': 0.649000},
            '3kcp': {'score': -311.021000, 'packstat': 0.713000},
            '3ul4': {'score': -274.482000, 'packstat': 0.719000},
            '4dh2': {'score': -289.640000, 'packstat': 0.675000},
            '4fl4': {'score': -279.734000, 'packstat': 0.672000},
            '4fl5': {'score': -319.303000, 'packstat': 0.646000},
            '4iu2': {'score': -366.587000, 'packstat': 0.674000},
            '4uyp': {'score': -276.661000, 'packstat': 0.739000},
            '4uyq': {'score': -277.328000, 'packstat': 0.698000},
            '5new': {'score': -297.519000, 'packstat': 0.679000},

            '1anu': {'score': -286.493000, 'packstat': 0.706000},
            '1aoh': {'score': -309.990000, 'packstat': 0.763000},
            '1g1k': {'score': -273.151000, 'packstat': 0.719000},
            '1qzn': {'score': -312.724000, 'packstat': 0.648000},
            '1tyj': {'score': -323.999000, 'packstat': 0.711000},
            '2bm3': {'score': -330.232000, 'packstat': 0.706000},
            '2jh2': {'score': -251.190000, 'packstat': 0.667000},
            '2vo8': {'score': -254.734000, 'packstat': 0.686000},
            '2w1n': {'score': -247.025000, 'packstat': 0.679000},
            '2w5f': {'score': -258.077000, 'packstat': 0.677000},
            '2xbt': {'score': -290.186000, 'packstat': 0.603000},
            '2xdh': {'score': -237.370000, 'packstat': 0.605000},
            '2zf9': {'score': -317.812000, 'packstat': 0.652000},
            '3bwz': {'score': -323.190000, 'packstat': 0.700000},
            '3fnk': {'score': -319.249000, 'packstat': 0.645000},
            '3l8q': {'score': -284.827000, 'packstat': 0.662000},
            '4jo5': {'score': -328.457000, 'packstat': 0.752000},
            '4ums': {'score': -258.205000, 'packstat': 0.684000},
            }

    docs = {'1ohz': {'score': -84.981000, 'packstat': 0.632000},
            '2b59': {'score': -301.66300, 'packstat': 0.653000},
            '2ccl': {'score': -93.799000, 'packstat': 0.606000},
            '2ozn': {'score': -258.16000, 'packstat': 0.639000},
            '2vn5': {'score': -21.212000, 'packstat': 0.483000},
            '2vn6': {'score': -92.409000, 'packstat': 0.550000},
            '2y3n': {'score': -21.242000, 'packstat': 0.479000},
            '3kcp': {'score': -132.19300, 'packstat': 0.663000},
            '3ul4': {'score': -90.540000, 'packstat': 0.498000},
            '4dh2': {'score': -90.261000, 'packstat': 0.648000},
            '4iu2': {'score': -105.470000, 'packstat': 0.642000},
            '4fl4': {'score': -13.360000, 'packstat': 0.632000},
            '4fl5': {'score': 187.356000, 'packstat': 0.648000},
            '4uyp': {'score': 95.3980000, 'packstat': 0.702000},
            '4uyq': {'score': -126.87700, 'packstat': 0.670000},
            '5new': {'score': -89.463000, 'packstat': 0.641000},
    }

    return cohs, docs


def dimer_data():
    """
    :return: a dictionary of thresholds measured on dimers
    """
    return {'a_score': -394.147000,
            'a_sasa': 1464, 'a_sasa_4iu2': 919.987,
            'a_shape': 0.623000,
            'a_ddg': -22.643, 'a_ddg_4iu2': -17.442,
            'a_packstat': 0.641000,
            'a_buried_2': 3} # 'total_score': 69640.368000,


def generate_run_filters(args=None, filter_thresholds=None) -> RunFilters:
    """
    the original thresholds were: ddg -22.5, score -390, sasa 1460, shape 0.623, packstat 0.641, buried_2 3
    :rtype : RunFilters
    """
    if args is None:
        args = {'ddg': -18.0, 'sasa': 1200, 'shape': 0.6, 'packstat': 0.6}
    run_filters = RunFilters()
    if filter_thresholds is None:
        run_filters.append_filter(Filter(name='a_ddg', typ='ddg', threshold=-abs(args['ddg']), limits=[-10000, 10000],
                                         under_over='under', g_name='$\Delta$$\Delta$G'))
        run_filters.append_filter(Filter(name='a_score', typ='score', threshold=-390.0, limits=[-10000, 10000],
                                         under_over='under', g_name='Score'))
        run_filters.append_filter(Filter(name='a_sasa', typ='sasa', threshold=args['sasa'], limits=[0, 100000],
                                         under_over='over', g_name='SASA'))
        run_filters.append_filter(Filter(name='a_shape', typ='shape', threshold=args['shape'], limits=[0.0, 1.0],
                                         under_over='over', g_name='Shape Complementarity'))
        run_filters.append_filter(Filter(name='a_packstat', typ='packstat', threshold=args['packstat'], limits=[0.0, 1.0],
                                         under_over='over', g_name='PackStat'))
        run_filters.append_filter(Filter(name='a_buried_2', typ='buried_2', threshold=args['buried_2'], limits=[0, 100],
                                         under_over='under', g_name='UnsatisfiedHBonds'))
        run_filters.append_filter(Filter(name='a_rms', typ='rmsd', threshold=1000, limits=[0, 1000],
                                         under_over='under', g_name='RMSD'))
        run_filters.append_filter(Filter(name='a_hbonds', typ='hbonds', threshold=-abs(args['hbonds']),
                                         limits=[-1000, 1000], under_over='under', g_name='H. bonds'))
        # run_filters.append_filter(Filter(name='coh_packstat', typ='packstat', threshold=))
    return run_filters


def score_passes(score, args=None):
    """
    :param score: a score dict {'filter_name': float(result)}
    :param args: run arguments
    :return: True iff score passes all thresholds
    """
    if args['cohs_data'] is not None:
        cohs_data = args['cohs_data']
        docs_data = args['docs_data']
        dimers_data = args['dimer_data']
    else:
        cohs_data, docs_data = monomer_data()
        dimers_data = dimer_data()
        args['ignore_monomers'] = False
    try:
        desc_s = score['description'].split('_')
        on_is_at = desc_s.index('on')
        coh_bb = desc_s[on_is_at + 1]
        doc_bb = desc_s[on_is_at + 3]
    except:
        desc_s = score['description'].split('_')
        A_is_at = desc_s.index('A')
        coh_bb = desc_s[A_is_at - 1]
        doc_bb = desc_s[A_is_at + 1]
    if not args['ignore_monomers']:
        if score['coh_score'] >= cohs_data[coh_bb]['score'] + 100 or \
                        score['coh_packstat'] <= cohs_data[coh_bb]['packstat'] - 0.1 or \
                        score['doc_score'] >= docs_data[doc_bb]['score'] + 100 or \
                        score['doc_packstat'] <= docs_data[doc_bb]['packstat'] - 0.1:
    #     print('fail stability'
    #     print('coh_score', score['coh_score'], cohs_data[coh_bb]['score'], score['coh_score'] >= cohs_data[coh_bb]['score'] + 50
    #     print('coh_packstat', score['coh_packstat'], cohs_data[coh_bb]['packstat'], score['coh_packstat'] <= cohs_data[coh_bb]['packstat'] - 0.1
    #     print('doc_score', score['doc_score'], docs_data[doc_bb]['score'], score['doc_score'] >= docs_data[doc_bb]['score'] + 50
    #     print('doc_packstat', score['doc_packstat'], docs_data[doc_bb]['packstat'], score['doc_packstat'] <= docs_data[doc_bb]['packstat'] - 0.1
    #     # raw_input()
            return False
    # and score['total_score'] >= dimers_data['total_score'] - 100\
    # and score['total_score'] <= dimers_data['total_score'] - 0\
    if score['a_score'] <= dimers_data['a_score'] + 0\
            and score['a_sasa'] >= dimers_data['a_sasa'] - 0\
            and score['a_shape'] >= dimers_data['a_shape'] - 0.\
            and score['a_ddg'] <= dimers_data['a_ddg'] + 0\
            and score['a_packstat'] >= dimers_data['a_packstat'] - 0\
            and score['a_buried_2'] <= dimers_data['a_buried_2']:
        return True
    # and score['total_score'] <= dimers_data['total_score'] - 0\
    elif coh_bb == '4iu2' \
            and score['a_score'] <= dimers_data['a_score'] + 0\
            and score['a_sasa_4iu2'] >= dimers_data['a_sasa'] - 0\
            and score['a_shape'] >= dimers_data['a_shape'] - 0.\
            and score['a_ddg_4iu2'] <= dimers_data['a_ddg'] + 0\
            and score['a_packstat'] >= dimers_data['a_packstat'] - 0.\
            and score['a_buried_2'] <= dimers_data['a_buried_2']:
        return True
    else:
        # print('FAIL'
        # print('a_score', score['a_score'], dimers_data['a_score'], score['a_score'] <= dimers_data['a_score'] + 150
        # print('a_sasa', score['a_sasa'], dimers_data['a_sasa'], score['a_sasa'] >= dimers_data['a_sasa'] - 100
        # print('a_shape', score['a_shape'], dimers_data['a_shape'], score['a_shape'] >= dimers_data['a_shape'] - 0.1
        # print('total_score', score['total_score'], dimers_data['total_score'], score['total_score'] >= dimers_data['total_score'] - 100
        # print('a_ddg', score['a_ddg'], dimers_data['a_ddg'], score['a_ddg'] <= dimers_data['a_ddg'] + 5
        # print('a_packstat', score['a_packstat'], dimers_data['a_packstat'], score['a_packstat'] >= dimers_data['a_packstat'] - 0.1
        # raw_input()
        return False


def all_who_pass(args, score_dict=None):
    """
    :param args: run arguments
    :return: a list of scores that passed the thresholds
    """
    if score_dict is None:
        score_dict = score2dict_new(args['score_file'])
    passed = []
    for name, entry in score_dict.items():
        if score_passes(entry, args):
            passed.append(entry)
    return passed


def all_who_pass_run_filters(args: dict, score_dict: dict, runfilters: RunFilters) -> (dict, dict):
    """
    :param args: run arguments
    :param score_dict: score dictionary
    :param runfilters: run filters
    :return: a list of dict of scores
    """
    passed, failed = {}, {}
    msg = False
    for k, v in score_dict.items():
        res, msg = runfilters.test_all(v)
        if res:
            passed[k] = v
        else:
            failed[k] = v
    if msg:
        print(msg)
    return passed, failed


def best_n_structures(args, score_dict=None):
    """
    :param args: run arguments. requires 'filter' and 'n'
    :param score_dict: if None, will make score dict.
    :return: names of n best structrues by filter
    """
    sorted_score = sorted(score_dict.values(), key=lambda x: x[args['filter']])
    return sorted_score[:args['n']]


def average_filter(args):
    """
    :param args: run arguments
    :return: prints statistics of args['filter'], and plots it as a boxplot with scatter
    """
    import matplotlib.pyplot as plt
    import numpy as np
    dimers_data = dimer_data()
    score_dict = score2dict_new(args['score_file'])
    resutls = []
    over, under = 0, 0
    for entry in score_dict.values():
        resutls.append(entry[args['filter']])
        under += 1 if entry[args['filter']] <= dimers_data[args['filter']] else 0
        under += 1 if entry[args['filter']] >= dimers_data[args['filter']] else 0
        print(entry[args['filter']])
    print('average of %s is %f' % (args['filter'], np.average(resutls)))
    print('stdev   of %s is %r' % (args['filter'], np.std(resutls)))
    print('max     of %s is %r' % (args['filter'], max(resutls)))
    print('min     of %s is %r' % (args['filter'], min(resutls)))
    print('for filter threshold %f' % dimers_data[args['filter']])
    print('over     %s    %i' % (args['filter'], over))
    print('under    %s    %i' % (args['filter'], under))
    plt.boxplot(resutls)
    plt.plot(np.random.normal(1, 0.04, size=len(resutls)), resutls, 'r.', alpha=0.05)
    plt.axhline(y=dimers_data[args['filter']])
    plt.show()


def normal_name(old_name):
    """
    :param old_name: a filter name
    :return: a presentable name
    """
    rename = {'ddg': '$\Delta$$\Delta$G', 'rms': 'RMS'}
    new_name = old_name
    if old_name[:2] == 'a_':
        new_name = old_name[2:]
    if new_name in rename.keys():
        return rename[new_name]
    return new_name.upper()


def this_vs_that(args: dict, runfilters: RunFilters, passed: dict, failed: dict, score_dict: dict, plot_num: int=None,
                 how_many_plots: int=1):
    """
    :param args: run rguments
    :return: plots a scatter of args['x'] Vs. args['y']
    """
    if plot_num is not None:
        if args['mode'] == 'jack_matrix':
            plt.subplot(440+plot_num)
        else:
            plt.subplot(how_many_plots*100+10+plot_num)
    x_passed, y_passed, x_failed, y_failed = [], [], [], []
    for k, v in passed.items():
        x_passed.append(v[args['x']])
        y_passed.append(v[args['y']])
    for k, v in failed.items():
        x_failed.append(v[args['x']])
        y_failed.append(v[args['y']])

    plt.scatter(x_failed, y_failed, marker='.', c='b', alpha=0.1, s=50)
    plt.scatter(x_passed, y_passed, marker='o', c='purple', alpha=0.3, s=80)

    # plt.xlim([runfilters[args['x']].min_tested, runfilters[args['x']].max_tested])
    # plt.ylim([runfilters[args['y']].min_tested, runfilters[args['y']].max_tested])
    all_together_x = x_passed + x_failed
    all_together_y = y_passed + y_failed
    plt.xlim([np.percentile(all_together_x, 0.01), np.percentile(all_together_x, 99.9)])
    plt.ylim([np.percentile(all_together_y, 0.01), np.percentile(all_together_y, 99.9)])
    # plt.xlim([0, 30])
    # plt.ylim([-30, 0])
    plt.xlabel(runfilters[args['x']].g_name)
    plt.ylabel(runfilters[args['y']].g_name)

    plt.title('%s Vs. %s for %s over %i structures' % (normal_name(args['x']), normal_name(args['y']),
                                                       args['score_file'], len(list(score_dict.keys()))))
    try:
        plt.suptitle('%i passed, %f percent' % (len(x_passed), 100.*float(len(x_passed))/float(len(x_failed))))
    except:
        pass
    if args['show_fig']:
        plt.show()
    if args['save_fig']:
        plt.savefig('%s.png' % args['score_file'].split('.')[0])


def multiple_plots(args: dict, runfilters: RunFilters, passed: dict, failed: dict, score_dict: dict):
    args['show_fig'] = False
    plt.figure(1)
    for x, y, i in zip(['sasa', 'shape', 'sasa'], ['ddg', 'packstat', 'score'], [1, 2, 3]):
        args['x'] = x
        args['y'] = y
        this_vs_that(args, runfilters, passed, failed, score_dict, plot_num=i, how_many_plots=3)
    plt.show()


def multiple_plots_rmsd(args: dict, runfilters: RunFilters, passed: dict, failed: dict, score_dict: dict):
    args['show_fig'] = False
    plt.figure(1)
    args['x'] = 'rmsd'
    filters = ['ddg', 'packstat', 'score', 'packstat', 'buried_2', 'sasa']
    for i, y in enumerate(filters):
        args['y'] = y
        this_vs_that(args, runfilters, passed, failed, score_dict, plot_num=i+1, how_many_plots=len(filters))
    plt.show()


def find_thresholds_by_rmsd(args):
    """
    :param args: run args
    :return: minimal/maximal (depending on threshold type) of the different filters that will pass all structures with
    RMSD under args['rmsd_threshold'].
    """
    score_dict = score2dict(args['score_file'])
    filter_thresholds = dict(a_score=-100000., a_sasa=100000., a_shape=1.0, total_score=1000000., a_ddg=-100.,
                             a_packstat=1.0, a_buried_2=0)
    passed_rmsd = []
    for name, sc in score_dict.items():
        if sc['rmsd'] <= args['rmsd_threshold']:
            filter_thresholds['a_score'] = max([filter_thresholds['a_score'], sc['score']])
            filter_thresholds['a_sasa'] = min([filter_thresholds['a_sasa'], sc['sasa']])
            filter_thresholds['a_shape'] = min([filter_thresholds['a_shape'], sc['shape']])
            filter_thresholds['total_score'] = min([filter_thresholds['total_score'], sc['total_score']])
            filter_thresholds['a_ddg'] = max([filter_thresholds['a_ddg'], sc['ddg']])
            filter_thresholds['a_packstat'] = min([filter_thresholds['a_packstat'], sc['packstat']])
            filter_thresholds['a_buried_2'] = max([filter_thresholds['a_buried_2'], sc['buried_2']])
            passed_rmsd.append(sc)
    print('found %i scores with rmsd <= %f' % (len(passed_rmsd), args['rmsd_threshold']))
    # print('the old thresholds were:\n%s' % '\n'.join(['%s %f' % (k, v) for k, v in dimer_data().items()]))
    print('the old thresholds were:', generate_run_filters().report())
    print('defined these new filters:\n%s' % '\n'.join(['%s %f' % (k, v) for k, v in filter_thresholds.items()]))
    args['dimer_data'] = filter_thresholds
    run_filters_updated = RunFilters()
    run_filters_updated.append_filter(Filter(name='a_ddg', typ='ddg', threshold=filter_thresholds['a_ddg'],
                                             limits=[-10000, 10000], under_over='under', g_name='$\Delta$$\Delta$G'))
    run_filters_updated.append_filter(Filter(name='a_score', typ='score', threshold=filter_thresholds['a_score'],
                                             limits=[-10000, 10000], under_over='under', g_name='Score'))
    run_filters_updated.append_filter(Filter(name='a_sasa', typ='sasa', threshold=filter_thresholds['a_sasa'],
                                             limits=[0, 100000], under_over='over', g_name='SASA'))
    run_filters_updated.append_filter(Filter(name='a_shape', typ='shape', threshold=filter_thresholds['a_shape'],
                                             limits=[0.0, 1.0], under_over='over', g_name='Shape Complementarity'))
    run_filters_updated.append_filter(Filter(name='a_packstat', typ='packstat',
                                             threshold=filter_thresholds['a_packstat'], limits=[0.0, 1.0],
                                             under_over='over', g_name='PackStat'))
    run_filters_updated.append_filter(Filter(name='a_buried_2', typ='buried_2',
                                             threshold=filter_thresholds['a_buried_2'], limits=[0, 100],
                                             under_over='under', g_name='UnsatisfiedHBonds'))
    run_filters_updated.append_filter(Filter(name='a_rms', typ='rmsd', threshold=1000, limits=[0, 1000],
                                             under_over='under', g_name='RMSD'))
    passed, failed = all_who_pass_run_filters(args, score_dict, run_filters_updated)
    # this_vs_that(args, run_filters_updated, passed, failed, score_dict)
    multiple_plots(args, run_filters_updated, passed, failed, score_dict)
    args['x'], args['y'] = 'rmsd', 'ddg'
    this_vs_that(args, run_filters_updated, passed, failed, score_dict)
    plt.show()


def spread_filters(args: dict, run_filters: RunFilters, score_dict: dict):
    plt.figure(1)
    col_num = len(list(run_filters.keys()))
    lims = {'sasa': [0, 2000], 'ddg': [-40, 5], 'packstat': [0.0, 1.0], 'shape': [0.0, 1.0], 'score': [-500, 10],
            'rmsd': [0.0, 30.0]}
    for i, flt in enumerate(run_filters.filters.values()):
        plt.subplot(1, col_num, i+1)
        plt.title(flt.filter_type)
        if flt.filter_type == 'buried_2':
            plt.hist(flt.all_seen)
        else:
            plt.boxplot(flt.all_seen)
            plt.ylim(lims[flt.filter_type])
    plt.show()


def score_diagonal(args: dict):
    results = {}
    run_filters = generate_run_filters(args)
    sc_files = [a for a in os.listdir(args['score_dir']) if a[-6:] == '.score']
    for sc_file in sc_files:
        score_dict = score2dict(sc_file)
        passed, failed = all_who_pass_run_filters(args, score_dict, run_filters)
        results[sc_file] = len(list(passed.keys()))
        print(sc_file, len(list(passed.keys())))


def what_coh(description: str, args=None) -> str:
    """
    :param description: score description
    :return: coh_name
    """
    if not args:
        args= dict()
    if args['naming'] == 'coh_on_doc':
        return description.split(('_on_'))[0].split('_')[1]
    ynum = re.search('y[0-9]{3}', description)

    if '_ON_' in description:
        return description.split('_ON_')[0]
    if ynum:
        return ynum.group(0)
    if '_on_' not in description:
        return description.split('_')[1]
    s = description.split('_')
    coh = '_'.join(s[1:3])
    return coh


def what_doc(description: str, args=None) -> str:
    """
    :param args: run arguments
    :param description: score description
    :return: doc_name
    """
    if not args:
        args = dict()
    if args['naming'] == 'coh_on_doc':
        return description.split(('_on_'))[1].split('_')[0]
    desc = re.compile('dz_[0-9a-z]{4}_[Aa]_[0-9a-z]{4}_')
    ynum = re.search('y[0-9]{3}', description)

    if 'ctcipadoc' in description:
        return 'ctcipadoc'
    if '_ON_' in description:
        return description.split('_ON_')[1]
    if ynum:
        return ynum.group(0)
    if '_on_' not in description:
        return description.split('_')[1]
    s = description.split('_')
    doc = '_'.join(s[4:6])
    return doc


def analyse_matrix(args: dict):
    # run_filters = generate_run_filters(args={'ddg': 22.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6})
    # run_filters = generate_run_filters(args={'ddg': 24.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6})
    # run_filters = generate_run_filters(args={'ddg': 25.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6})
    run_filters = generate_run_filters(args={'ddg': 24.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6, 'buried_2': 3,
                                             'hbonds': -10.})
    # run_filters = generate_run_filters(args={'ddg': 25.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.65, 'buried_2': 3})
    # print(run_filters.report())
    score_files = [a for a in os.listdir('./') if a[-6:] == '.score']
    coh_name_list = sorted(list(set([what_coh(a, args) for a in score_files])))
    doc_name_list = sorted(list(set([what_doc(a, args) for a in score_files])))
    # print(coh_name_list)
    # print(doc_name_list)
    # df = DataFrame([{coh: Series([-1] * len(doc_name_list))} for coh in coh_name_list], index=doc_name_list)
    df = DataFrame([{coh: -1 for coh in coh_name_list}] * len(doc_name_list), index=doc_name_list)
    for sc in score_files:
        score_dict = score2dict(sc)
        passed, failed = all_who_pass_run_filters(args, score_dict, run_filters)
        coh = what_coh(sc, args)
        doc = what_doc(sc, args)
        df[coh][doc] = len(passed)
    # print(df)
    show_prediction_heat_map(df)
    
    
def show_prediction_heat_map(df: DataFrame):
    from matplotlib import colors
    from numpy import array, arange
    axis = plt.gca()
    cmap = colors.ListedColormap(['white', 'cornflowerblue', 'darkturquoise', 'darkorange'])
    bounds = [-100, 0, 10, 20, 100]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    heatmap = plt.pcolor(array(df), cmap=cmap, norm=norm, edgecolors='k', linewidth=2)
    print(df)

    for y in range(array(df.shape)[0]):
        for x in range(array(df.shape)[1]):
            if array(df)[y, x] >= 0:
                plt.text(x+0.5, y+0.5, array(df)[y, x], horizontalalignment='center', verticalalignment='center')

    plt.yticks(arange(0.5, len(df.index), 1), df.index)
    plt.xticks(arange(0.5, len(df.columns), 1), df.columns, rotation=70)

    # plt.yticks(arange(0.5, len(df.index), 1), [official_names(n, {'naming': 'coh_on_doc'}, 'coh') for n in df.index])
    # plt.xticks(arange(0.5, len(df.columns), 1), [official_names(n, {'naming': 'coh_on_doc'}, 'doc') for n in df.columns],
    #            rotation=70)
    plt.xlabel('Cohesin name', style='oblique')
    plt.ylabel('Dockerin name', style='oblique')
    axis.set_aspect('equal')
    plt.title('Cohesin dockerin cross binding')
    plt.show()


def mock_matrix_prediction():
    names = ['A', 'B', 'C', 'D', 'E', 'F']#  , 'G', 'H', 'I', 'J']
    df = DataFrame([{coh: -1 for coh in names}] * len(names), index=names)
    for n1 in names:
        for n2 in names:
            df[n1][n2] = random.randrange(0, 30, 5)
            # if n1 == n2:
            #     df[n1][n2] = random.randrange(20, 25, 5)
    show_prediction_heat_map(df)


def official_names(name):
    if 'j' == name[0]:
        return name
    if name == 'ctcipadoc':
        return 'Ct CipA-Doc'
    else:
        s = name.split('_')
        if 'sca' in s[1]:
            return s[0][0].upper()+s[0][1]+' '+'Sca'+s[1][-1].upper()
        if 'cip' in s[1]:
            return s[0][0].upper()+s[0][1]+' '+'Cip'+s[1][-1].upper()
        if 'doc' == s[1]:
            return s[0][0].upper()+s[0][1]+' '+'Doc'
        return s[0][0].upper()+s[0][1]+' '+s[1].upper()


def multi_filters_plot(args):
    """
    :param args: run arguments, not used
    :return: draws multiple plots for the different run filter configurations (for my thesis)
    """
    from matplotlib import colors
    from numpy import array, arange

    cmap = colors.ListedColormap(['white', 'cornflowerblue', 'darkturquoise', 'darkorange'])
    bounds = [-100, 0, 10, 20, 100]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    fig = plt.figure(figsize=(6.02, 6.38))
    plt.subplots_adjust(left=0.25, bottom=0.15, right=None, top=None, wspace=0.1, hspace=0.2)

    various_filters = OrderedDict()
    various_filters['ddG'] = generate_run_filters(
        args={'ddg': 24.0, 'sasa': 0000, 'shape': 0., 'packstat': 0.0, 'buried_2': 30, 'hbonds': 0.})
    various_filters['ddG_SASA'] = generate_run_filters(
        args={'ddg': 24.0, 'sasa': 1400, 'shape': 0., 'packstat': 0.0, 'buried_2': 30, 'hbonds': 0.})
    various_filters['ddG_SASA_pack'] = generate_run_filters(
        args={'ddg': 24.0, 'sasa': 1400, 'shape': 0., 'packstat': 0.6, 'buried_2': 30, 'hbonds': 0.})
    various_filters['ddG_SASA_pack_shape'] = generate_run_filters(
        args={'ddg': 24.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6, 'buried_2': 30, 'hbonds': 0.})
    various_filters['ddG_SASA_pack_shape_buried'] = generate_run_filters(args={
        'ddg': 24.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6, 'buried_2': 3, 'hbonds': 0.})
    various_filters['ddG_SASA_pack_shape_buried_hbonds'] = generate_run_filters(
        args={'ddg': 24.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6, 'buried_2': 3, 'hbonds': -10.})
    num_letter = {1: 'A', 2: 'B', 3: 'C', 4: 'D', 5: 'E', 6: 'F'}

    score_files = [a for a in os.listdir('./') if a[-6:] == '.score']
    coh_name_list = sorted(list(set([what_coh(a) for a in score_files])))
    doc_name_list = sorted(list(set([what_doc(a) for a in score_files])))
    scores_dict = {}
    for sc in score_files:
        scores_dict[sc] = score2dict(sc)
    i = 1
    for name, filters in various_filters.items():
        df = DataFrame([{coh: -1 for coh in coh_name_list}] * len(doc_name_list), index=doc_name_list)
        for sc, sc_dict in scores_dict.items():
            passed, failed = all_who_pass_run_filters({}, sc_dict, filters)
            coh = what_coh(sc)
            doc = what_doc(sc)
            df[coh][doc] = len(passed)
        plt.subplot(3, 2, i)
        axis = plt.gca()
        heatmap = plt.pcolor(array(df), cmap=cmap, norm=norm, edgecolors='k', linewidth=2)
        for y in range(array(df.shape)[0]):
            for x in range(array(df.shape)[1]):
                if array(df)[y, x] >= 0:
                    plt.text(x+0.5, y+0.5, array(df)[y, x], horizontalalignment='center', verticalalignment='center', fontsize=6)
        # make sure labels are only on outer subplots
        if i in [1, 3, 5]:
            plt.yticks(arange(0.5, len(df.index), 1), [official_names(n) for n in df.index], fontsize=10)
        else:
            plt.yticks([])
        if i in [5, 6]:
            plt.xticks(arange(0.5, len(df.columns), 1), [official_names(n) for n in df.columns], rotation=70, fontsize=10)
        else:
            plt.xticks([])
        plt.title(num_letter[i])
        # plt.title(official_title(name))
        i += 1
        # axis.set_aspect('equal')
    fig.text(0.5, 0.04, 'Cohesin name', ha='center', va='center', fontsize=24)
    fig.text(0.06, 0.5, 'Dockerin name', ha='center', va='center', rotation='vertical', fontsize=24)
    plt.savefig('mini_postdiction.png', dpi=600)
    plt.show()


def official_title(title: str) -> str:
    official = {'ddG': '∆∆G', 'SASA': 'SASA', 'pack': 'PACKSTAT', 'shape': 'Shape Complementrity',
                'buried': 'Buried Unsatisfied H.bonds', 'hbonds': 'Hydrogen bonds'}
    s = title.split('_')
    result = []
    for n in s:
        result.append(official[n])
    return ' '.join(result)


def filter_result_for_passed(args):
    """
    analyses a specific filter behaviour for all passed structures
    :param args: run arguments
    :return:
    """
    run_filters = generate_run_filters(args={'ddg': 25.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6})
    score_files = sorted(args['score_files'])
    results = {}
    for i in score_files:
        results[i] = {}
        score_dict = score2dict(i)
        passed, failed = all_who_pass_run_filters(args, score_dict, run_filters)
        results[i]['filter_passed'] = []
        for name, psd in passed.items():
            results[i]['filter_passed'].append(psd[args['filter']])
    plt.boxplot([a['filter_passed'] for a in results.values()], labels=[a.split('.score')[0] for a in results.keys()])
    plt.show()


def jack_matrix(args):
    """
    displays a matrix of ddG Vs. RMSD
    :param args: run arguments
    :return:
    """
    args['show_fig'] = False
    fig = plt.figure()
    # fig = plt.figure(figsize=(8.27, 11.69))
    # run_filters = generate_run_filters(args={'ddg': 25.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6})
    run_filters = generate_run_filters(args={'ddg': -16, 'sasa': 1200, 'shape': 0.5, 'packstat': 0.5, 'buried_2': 30,
                                             'hbonds': 10})
    jk_scores = [a for a in os.listdir('./') if a[-6:] == '.score' and 'SD' not in a and 'no_docking' not in a]
    sd_scores = {}
    for sc_f in [a for a in os.listdir('./') if a[-6:] == '.score' and 'SD' in a]:
        if sc_f[4:8] not in sd_scores.keys():
            sd_scores[sc_f[4:8]] = {}
        sd_scores[sc_f[4:8]][sc_f[12:16]] = sc_f
    no_dock_files = {}
    for sc_f in [a for a in os.listdir('./') if 'no_docking' in a]:
        if sc_f[:4] not in no_dock_files.keys():
            no_dock_files[sc_f[:4]] = {}
        no_dock_files[sc_f[:4]][sc_f[11:15]] = sc_f
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.15, hspace=0.25)
    matplotlib.rcParams['axes.linewidth'] = 0.8
    font0 = FontProperties()
    font = font0.copy()
    font.set_weight('semibold')
    sorted_scores = sorted(jk_scores)
    z = sorted_scores[1]
    sorted_scores[1] = sorted_scores[-2]
    sorted_scores[-2] = z
    sorted_scores = [sorted_scores[0], sorted_scores[1]] + sorted(sorted_scores[2:])
    for num, sc in enumerate(sorted_scores):
        ax = plt.subplot(3, 3, 1+num)
        sc_dict = score2dict(sc)

        passed, failed = all_who_pass_run_filters(args, sc_dict, run_filters)
        x_passed, y_passed, x_failed, y_failed = [], [], [], []
        for k, v in passed.items():
            x_passed.append(v[args['x']])
            y_passed.append(v[args['y']])
        for k, v in failed.items():
            x_failed.append(v[args['x']])
            y_failed.append(v[args['y']])

        # draw simple docking results
        if sc[4:8] in sd_scores.keys():
            crystals = list(sd_scores[sc[4:8]].keys())
            sd_sc_dict = score2dict(sd_scores[sc[4:8]][crystals[0]])
            plt.scatter([a['rmsd'] for a in sd_sc_dict.values()], [a['ddg'] for a in sd_sc_dict.values()], marker='.',
                        c='lightgrey', s=30, linewidths=0)
            if len(crystals) > 1:
                sd_sc_dict = score2dict(sd_scores[sc[4:8]][crystals[1]])
                plt.scatter([a['rmsd'] for a in sd_sc_dict.values()], [a['ddg'] for a in sd_sc_dict.values()],
                            marker='.', c='grey', s=30, linewidths=0)
        # draw DoCohModeller results
        plt.scatter(x_failed, y_failed, marker='.', c='blue', alpha=0.6, s=50, linewidth=0.2)
        plt.scatter(x_passed, y_passed, marker='.', c='red', alpha=0.6, s=100, linewidths=0.3)

        # draw no docking results
        no_dock_rmsd_1, no_dock_ddg_1 = get_no_docking_results(sc, list(no_dock_files[sc[4:8]].values())[0])
        print('marking %s with %f, %f' % (sc[4:8], no_dock_rmsd_1, no_dock_ddg_1))
        plt.scatter(no_dock_rmsd_1, no_dock_ddg_1, c='lightgrey', marker='^', s=60)
        if len(no_dock_files[sc[4:8]].keys()) > 1:
            no_dock_rmsd_2, no_dock_ddg_2 = get_no_docking_results(sc, list(no_dock_files[sc[4:8]].values())[1])
            print('marking %s with %f, %f' % (sc[4:8], no_dock_rmsd_2, no_dock_ddg_2))
            plt.scatter(no_dock_rmsd_2, no_dock_ddg_2, c='grey', marker='^', s=60)

        plt.xlim([0, 20])
        plt.ylim([-32.5, 0])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        plt.yticks([0, -10, -20, -30], fontsize=16)
        plt.xticks(fontsize=16)
        if num not in [0, 3, 6]:
            plt.setp(ax.get_yticklabels(), visible=False)
        if num not in [6, 7, 8]:
            plt.setp(ax.get_xticklabels(), visible=False)
        plt.title(sc.split('_')[1].upper(), fontsize=16, fontproperties=font)
    fig.text(0.5, 0.04, r'RMSD ($\AA$)', ha='center', va='center', fontsize=24)
    fig.text(0.06, 0.5, '∆∆G (R.E.U.)', ha='center', va='center', rotation='vertical', fontsize=24)
    plt.savefig('jack.png', dpi=100)
    plt.show()


def get_no_docking_results(sc_name: str, no_dock_file: str) -> (float, float):
    """
    :param sc_name: score name
    :param no_dock_file: score file
    :return: rmsd, ddg
    """
    sc_dict = score2dict(no_dock_file)
    min_rmsd, min_ddg = 1000., 1000.
    for sc in sc_dict.values():
        if sc['rmsd'] < min_rmsd:
            min_rmsd, min_ddg = sc['rmsd'], sc['ddg']
    return min_rmsd, min_ddg


def parse_names() -> dict:
    result = {}
    for l in open('/home/labs/fleishman/jonathaw/no_backup/designs/multi_docs_15Oct/recliques_4Nov/clique_6_pdbs/stabilisations_8Nov/all_stabilised_8Nov/names_list.txt', 'r'):
        s = l.split()
        if len(s) >= 2:
            result[s[0].split('.pdb.gz')[0]] = s[1]
    return result


def analyse_minidiagonal(args):
    with open('../minidiagonal.txt', 'w') as fout:
        run_filters = generate_run_filters(args={'ddg': 24.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6, 'buried_2': 3,
                                                 'hbonds': 10})
        counter = 0
        score_files = [a for a in os.listdir('./') if a[-3:] == '.sc']
        for sc in score_files:
            score_dict = score2dict(sc)
            passed, failed = all_who_pass_run_filters(args, score_dict, run_filters)
            if len(passed) > 5:
                fout.write('%s\t%i\n' % (sc, len(passed)))
                counter += 1
    print('%i passed minidiagonal' % counter)


def prism(args):
    run_filters = generate_run_filters(args={'ddg': 16.0, 'sasa': 1200, 'shape': 0.6, 'packstat': 0.6, 'buried_2': 3,
                                             'hbonds': 10})
    score_dict = score2dict(args['score_file'])
    passed, failed = all_who_pass_run_filters(args, score_dict, run_filters)
    print('RMSD\tPassed\tFailed')
    for v in passed.values():
        print('%3.3f\t%3.3f' % (v['rmsd'], v['ddg']))
    for v in failed.values():
        print('%3.3f\t\t%3.3f' % (v['rmsd'], v['ddg']))


def prediction_results(args):
    sc_files = [a for a in os.listdir() if a[-6:] == '.score']
    run_filters = generate_run_filters(args={'ddg': 18.0, 'sasa': 1300, 'shape': 0.6, 'packstat': 0.6,
                                                 'buried_2': 3, 'hbonds': 0})
    for sc_file in sc_files:
        sc_dict = score2dict(sc_file)
        passed, failed = all_who_pass_run_filters(args, sc_dict, run_filters)
        names = sc_file.split('_on_')
        coh = names[0].split('all_')[1]
        doc = names[1].split('_0')[0]
        print("%-10s %-10s %-4i %-3f %i" % (coh, doc, len(passed), 100*len(passed)/(len(failed)+len(passed)), len(sc_dict.keys())))


if __name__ == '__main__':
    import argparse
    import sys
    first = sys.argv[1]
    parser = argparse.ArgumentParser()
    parser.add_argument('-score_file', default=first)
    parser.add_argument('-mode', default='analyse')
    parser.add_argument('-filter', default='ddg')
    parser.add_argument('-x', default='rmsd')
    parser.add_argument('-y', default='ddg')
    parser.add_argument('-rmsd_threshold', default=5.0, type=float)
    parser.add_argument('-cohs_data', default=None)
    parser.add_argument('-docs_data', default=None)
    parser.add_argument('-dimer_data', default=None)
    parser.add_argument('-n', type=int, default=5)
    parser.add_argument('-ignore_monomers', default=False)
    parser.add_argument('-show_fig', default=True)
    parser.add_argument('-save_fig', default=True)
    parser.add_argument('-show_all', default=False)
    parser.add_argument('-score_dir')
    parser.add_argument('-sasa', default=1200, type=int)
    parser.add_argument('-ddg', default=-18.0, type=float)
    parser.add_argument('-packstat', default=0.6, type=float)
    parser.add_argument('-shape', default=0.6, type=float)
    parser.add_argument('-buried_2', default=3)
    parser.add_argument('-hbonds', type=float, default=0)
    parser.add_argument('-score_files', nargs='*')
    parser.add_argument('-naming', default='coh_on_doc')

    args = vars(parser.parse_args())

    if args['mode'] == 'how_many_pass':
        score_dict = score2dict_new(args['score_file'])
        all_who = all_who_pass(args, score_dict)
        percentage = 100.0*float(len(all_who))/float(len(score_dict.values()))
        print('there were %i purples out of %i, which is %f' % (len(all_who),
                                                                len(score_dict.values()), percentage))
        if args['show_all']:
            for k in all_who:
                print(k['description']+'.pdb')

    elif args['mode'] == 'analyse':
        # args['run_filters'] = generate_run_filters()
        score_dict = score2dict(args['score_file'])
        passed = all_who_pass(args, score_dict)
        print('found %i purples' % len(passed))
        if len(passed) == 0:
            print('non passed, so the lowest ddg are:')
            best_n_structs = best_n_structures(args, score_dict)
        else:
            print('had passed, so best ddg out of those:')
            passed_score_dict = {a['description']: a for a in passed}
            best_n_structs = best_n_structures(args, passed_score_dict)
        print('the best %i strucutres by filter %s are \n%s.pdb' % (args['n'], args['filter'],
                                                                '.pdb '.join([a['description'] for a in best_n_structs])))

        if args['show_all']:
            print('showing ALL structures that passed the thresholds:')
            print('.pdb '.join([a['description'] for a in passed]) + '.pdb')
        # else:
        #     print('the best %i strucutres by filter %s are \n%s.pdb' % (args['n'], args['filter'],
        #                                                                 '.pdb '.join([a['description'] for a in
        #                                                                               best_n_structs]))
        this_vs_that(args, score_dict)

    elif args['mode'] == 'average':
        print(average_filter(args))

    elif args['mode'] == 'thresholds_by_rmsd':
        find_thresholds_by_rmsd(args)

    elif args['mode'] == 'test':
        score_dict = score2dict(args['score_file'])
        # run_filters = generate_run_filters(args)  # for use with low filters
        # for use with tested filters:
        run_filters = generate_run_filters(args={'ddg': 25.0, 'sasa': 1400, 'shape': 0.6, 'packstat': 0.6,
                                                 'buried_2': 3, 'hbonds': -10.})
        passed, failed = all_who_pass_run_filters(args, score_dict, run_filters)
        print('Found %i scores that passed out of %i scores' % (len(passed), len(list(score_dict.keys()))))
        if len(passed) == 0:
            print('non passed, so %i best by %s are' % (args['n'], args['filter']))
            best_n_structs = best_n_structures(args, score_dict)
        else:
            print('Found passed, so these are the %i best by %s from there' % (args['n'], args['filter']))
            best_n_structs = best_n_structures(args, passed)
        print('.pdb '.join([a['description'] for a in best_n_structs])+'.pdb')

        if args['show_all']:
            print('these are all that passed')
            for k in passed.keys():
                print('%s.pdb\n' % k)
        for flt in run_filters.values():
            print('a', flt)

        multiple_plots(args, run_filters, passed, failed, score_dict)
        run_filters.report()
        # args['x'] = 'rmsd'
        # args['y'] = 'ddg'
        # this_vs_that(args, run_filters, passed, failed, score_dict)
        # multiple_plots_rmsd(args, run_filters, passed, failed, score_dict)

    elif args['mode'] == 'spread_filters':
        run_filters = generate_run_filters(args)
        score_dict = score2dict(args['score_file'])
        passed, failed = all_who_pass_run_filters(args, score_dict, run_filters)
        spread_filters(args, run_filters, score_dict)

    elif args['mode'] == 'diagonal':
        score_diagonal(args)

    elif args['mode'] == 'analyse_matrix':
        analyse_matrix(args)

    elif args['mode'] == 'compare_filter_result_for_passed':
        filter_result_for_passed(args)

    elif args['mode'] == 'jack_matrix':
        jack_matrix(args)

    elif args['mode'] == 'a':
        parse_names()

    elif args['mode'] == 'analyse_minidiagonal':
        analyse_minidiagonal(args)

    elif args['mode'] == 'prism':
        prism(args)

    elif args['mode'] == 'multi_filters':
        multi_filters_plot(args)

    elif args['mode'] == 'mock':
        mock_matrix_prediction()

    elif args['mode'] == 'prediction_results':
        prediction_results(args)

    else:
        print('no mode found')
