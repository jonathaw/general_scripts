#!/usr/bin/env python2.7
"""
script and classes for the analysis of Yeast Surface Display results, provided as csv
"""
import os
import sys
# import wx
import argparse
import numpy as np
import pandas as pd
import FlowCytometryTools as fct
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from collections import OrderedDict

channels = OrderedDict({'alexa': 'Alexa 488-A',
                        'apc': 'APC-A',
                        'fsca': 'FSC-A',
                        'ssca': 'SSC-A'})


class Gate:
    def __init__(self, type_, thresholds, region, param):
        self.type = type_
        self.thresholds = thresholds
        self.region = region
        self.parameter = param

    def df_pass(self, df, passed):
        if self.type == 'line':
            if self.region == 'above' and passed:
                return df[df[self.parameter] > self.thresholds[0]].copy()
            elif self.region == 'above' and not passed:
                return df[df[self.parameter] < self.thresholds[0]].copy()
            elif self.region == 'below' and  passed:
                return df[df[self.parameter] < self.thresholds[0]].copy()
            elif self.region == 'below' and  not passed:
                return df[df[self.parameter] > self.thresholds[0]].copy()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', type=str)
    parser.add_argument('-folder')
    parser.add_argument('-name')
    parser.add_argument('-ligand_conc_nM', nargs='+', type=float)
    parser.add_argument('-well')
    parser.add_argument('-desc_file')

    args = vars(parser.parse_args())

    if args['mode'] == 'Kd':
        determine_Kd(args)

    elif args['mode'] == 'plate':
        plate_analysis(args)

    elif args['mode'] == 'interactive':
        interactive_view(args)

    else:
        print('no mode found')


def interactive_view(args):
    # read in entire plate, drop NAs
    file_position_dict = get_file_name_to_well_name_dict(args)
    plate = fct.FCPlate.from_dir(ID=args['name'], path=args['folder'], parser=file_position_dict,
                                 position_mapper='name')
    plate.dropna()

    plate[args['well']].view_interactively()


def plate_analysis(args):
    # get plate description
    plate_desc, pos_dict = parse_plate_description_file(args['desc_file'])
    plate_dict = get_plate_dfs(args)

    # normal alexa and apc to positives
    min_alexa = min([min(df[channels['alexa']].values) for df in plate_dict.values()])
    for well in plate_dict.keys():
        plate_dict[well][channels['alexa']] = plate_dict[well][channels['alexa']] - min_alexa

    min_apc = min([min(df[channels['apc']].values) for df in plate_dict.values()])
    for well in plate_dict.keys():
        plate_dict[well][channels['apc']] = plate_dict[well][channels['apc']] - min_apc

    # set expression gate
    expression_gate = Gate(type_='line', thresholds=[10000.0], region='above', param=channels['alexa'])

    for well, df in plate_dict.items():
        # well = 'G1'
        # df = plate_dict[well]

        df_passed = expression_gate.df_pass(df, passed=True)
        df_not_passed = expression_gate.df_pass(df, passed=False)

        perc_passed = float(len(df_passed)) / float(len(df))
        passed_median_apc = np.median(df_passed[channels['apc']])
        print '%s: passed: %.2f, percentage: %.2f, median APC: %.2f' % (well, len(df_passed), 100.*perc_passed,
                                                                        passed_median_apc)
        draw_scats_hists(df_passed, df_not_passed, well, expression_gate)
        # break


def determine_Kd(args):
    pwd = '/Volumes/labs/fleishman/jonathaw/YSD/13Apr_5711.A_5517.A/j5711.A_829-4518/96_Well-U_bottom_002/'
    pwd = './'
    file_name = pwd+'Specimen_002_A1_A01_001.fcs'
    sample = fct.FCMeasurement(ID='test', datafile=file_name)
    print(sample.channel_names)
    tsample = sample.transform('hlog', channels=[alexa, apc, fsca, ssca])
    ax = plt.gca()
    ax.scatter(tsample.data[alexa], tsample.data[apc], c='blue', alpha=0.05, edgecolors='none')

    plt.show()


def get_plate_dfs(args ,verbose=False):
    """
    analyse a plate of fcs files, return Kds etc.
    """
    # read in entire plate, drop NAs
    file_position_dict = get_file_name_to_well_name_dict(args)
    plate = fct.FCPlate.from_dir(ID=args['name'], path=args['folder'], parser=file_position_dict, position_mapper='name')
    plate.dropna()

    # gate by interactive gate
    gate1 = fct.PolyGate([(8.815e+03, 7.037e+03), (3.385e+04, 1.520e+04), (6.450e+04, 4.241e+04), (6.705e+04, 7.030e+04),
                          (2.976e+04, 4.649e+04), (4.218e+03, 1.112e+04), (7.794e+03, 7.717e+03)], ('FSC-A', 'SSC-A'),
                         region='in', name='gate1')
    plate = plate.gate(gate1)

    processed_dfs = OrderedDict()
    for well in get_viable_wells(plate.data.keys()):
        df = plate.data[well].data.copy()
        if verbose:
            print 'looking at %s' % well
            print '\toriginally have %i rows' % len(df)

        df.dropna()

        if verbose:
            for name, channel in channels.items():
                print '\tchannel %s, mean %.2e, median %.2e, std %.2e, over %i point' % \
                      (name, np.mean(df[channel]), np.median(df[channel]), np.std(df[channel]), len(df))

        processed_dfs[well] = df
    return processed_dfs


def analyse_plate(args):
    """
    analyse a plate of fcs files, return Kds etc.
    """
    # get plate description
    plate_desc, pos_dict = parse_plate_description_file(args['desc_file'])
    meds_df = pd.DataFrame(columns=plate_desc['column_strains'], index=plate_desc['conc'])

    # read in entire plate, drop NAs
    file_position_dict = get_file_name_to_well_name_dict(args)
    plate = fct.FCPlate.from_dir(ID=args['name'], path=args['folder'], parser=file_position_dict, position_mapper='name')
    plate.dropna()

    gate1 = fct.PolyGate([(2.250e+04, 1.135e+04), (5.858e+04, 2.024e+04), (8.131e+04, 5.304e+04),
                          (8.329e+04, 8.584e+04), (4.623e+04, 6.466e+04), (2.695e+04, 2.844e+04),
                          (2.250e+04, 1.204e+04)], ('FSC-A', 'SSC-A'), region='in', name='gate1')
    plate = plate.gate(gate1)

    processed_dfs = OrderedDict()
    for well in get_viable_wells(plate.data.keys()):
    # for well in ['%s2' % a for a in list('ABCDEFGH')]:
        print 'looking at %s' % well
        df = plate.data[well].data.copy()
        print '\toriginally have %i rows' % len(df)

        df = df.drop(df[df[channels['alexa']] < 0].index)
        df = df.drop(df[df[channels['apc']] < 0].index)
        df.dropna()

        for name, channel in channels.items():
            print '\tchannel %s, mean %.2e, median %.2e, std %.2e, over %i point' % (name, np.mean(df[channel]), np.median(df[channel]), np.std(df[channel]), len(df))
        expression_gate = 100.0
        # calc median by rejecting all points further that
        passed_df = df[df[channels['alexa']] > expression_gate]
        passed_df_apc = passed_df[channels['apc']]
        print 'number of points that passed threshold %i, %.2f' % (len(passed_df_apc), 100*len(passed_df_apc)/len(df))
        apc_median = np.median(reject_outliers(passed_df_apc.values))

        print('found the median of APC to be %.2e' % apc_median)
        # draw_scats_hists(df, well)
        meds_df.set_value(pos_dict[well]['conc'], pos_dict[well]['column_strain'], apc_median)

        processed_dfs[well] = passed_df
        # break
    print meds_df
    draw_median_plot(meds_df)
    # draw_well_set(processed_dfs)


def draw_well_set(dfs_dict):
    i = 1
    for well, df in dfs_dict.items():
        # print(well)
        # print(df)
        ax = plt.subplot(420+i)
        # plt.scatter(df[channels['alexa']], df[channels['apc']], marker='.')
        plt.hist(df[df[channels['apc']] < 2000][channels['apc']].values, bins=40)
        ax.set_title(well)
        # ax.set_xscale('log')
        # ax.set_yscale('log')
        ax.set_xlim([0, 2000])
        i += 1
    plt.show()


def reject_outliers(data, m = 2.):
    """
    rejecting outliers by distance from median
    """
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return data[s<m]


def draw_median_plot(meds_df):
    conc = meds_df.index.values
    strains = meds_df.columns.values
    curves = {strain: meds_df[strain].values for strain in strains}
    for name, crv in curves.items():
        ax1 = plt.subplot(111)
        plt.plot(conc, crv, label=name)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
    plt.legend()
    plt.show()


def draw_scats_hists(df_passed, df_not_passed, well, gate):
    fig = plt.figure()
    fig.suptitle(well)


    # draw FSC-A Vs. SSC-A scatter
    ax1 = plt.subplot(221)
    ax1.scatter(df_passed[channels['fsca']], df_passed[channels['ssca']], alpha=0.8, marker='.', color='r')
    ax1.scatter(df_not_passed[channels['fsca']], df_not_passed[channels['ssca']], alpha=0.8, marker='.', color='b')
    ax1.set_xlabel('FSC-A')
    ax1.set_ylabel('SSC-A')

    ax2 = plt.subplot(222)
    # ax2.hist(df_passed[channels['alexa']].values, bins=10**np.linspace(np.log10(0.01), np.log10(1000000), 50), color='r', edgecolor='r', log=True)
    # ax2.hist(df_not_passed[channels['alexa']].values, bins=10**np.linspace(np.log10(0.01), np.log10(1000000), 50), color='b', edgecolor='b', log=True)
    ax2.hist(df_passed[channels['alexa']].values, color='r', edgecolor='r')
    ax2.hist(df_not_passed[channels['alexa']].values, color='b', edgecolor='b')
    ax2.set_xlabel('Alexa 488')
    ax2.set_ylabel('Counts')
    # ax2.set_xscale('log')

    ax3 = plt.subplot(223)
    ax3.loglog(df_passed[channels['alexa']], df_passed[channels['apc']], '.', alpha=0.5, color='r')
    ax3.loglog(df_not_passed[channels['alexa']], df_not_passed[channels['apc']], '.', alpha=0.5, color='b')
    ax3.set_xlabel('Alexa 488')
    ax3.set_ylabel('APC')

    ax4 = plt.subplot(224)
    ax4.hist(df_passed[channels['apc']].values, bins=10**np.linspace(np.log10(0.01), np.log10(1000000), 50), color='r', edgecolor='r', log=True)
    ax4.hist(df_not_passed[channels['apc']].values, bins=10**np.linspace(np.log10(0.01), np.log10(1000000), 50), color='b', edgecolor='b', log=True)
    ax4.set_xlabel('APC')
    ax4.set_ylabel('Counts')
    ax4.set_xscale('log')


    plt.show()


def get_viable_wells(data_keys):
    """
    return ordered wells with data
    """
    result = []
    for n in range(1, 13):
        for a in list('ABCDEFGH'):
            if '%s%i' % (a, n) in data_keys:
                result.append('%s%i' % (a, n))
    return result


def temp():
    # transform and plot FSC-A and SSC-A
    # plate = plate.transform('tlog', channels=[channels['fsca'], channels['ssca']], th=0.1)
    # plate.plot([channels['fsca'], channels['ssca']], bins=100)
    # plt.show()

    # transform and plot alexa, gate for expression
    # plate = plate.transform('hlog', channels=[channels['alexa'], channels['apc']], b=10)
    # plate = plate.transform('tlog', channels=[channels['alexa'], channels['apc']], th=0.1)

    t_plate = plate.apply(log_transform, output_format='collection')
    print(plate['A1'].data)
    print(t_plate['A1'].data)

    alexa_gate = fct.ThresholdGate(1000.0, channels['alexa'], region='above')
    g_plate = t_plate.gate(alexa_gate)

    print 'non gated\n', t_plate.apply(count_events)
    print 'gated\n', plate.gate(alexa_gate).apply(count_events)

    # t_plate.plot(channels['alexa'], bins=10, gates=[alexa_gate])
    # plt.show()

    # show APC Vs. alexa
    t_plate.plot([channels['alexa'], channels['apc']], gates=[alexa_gate])
    plt.show()

    # cacl median APC for gated pop
    print 'median APC', g_plate.apply(calc_channel_median)


def log_transform(original_sample):
    """ This function implements a log transformation on the data. """
    # Copy the original sample
    new_sample = original_sample.copy()
    new_data = new_sample.data

    # Our transformation goes here
    new_data[channels['apc']] = np.log10(new_data[channels['apc']])
    new_data[channels['alexa']] = np.log10(new_data[channels['alexa']])
    new_data = new_data.dropna() # Removes all NaN entries
    new_sample.data = new_data


    return new_sample


def count_events(well):
    """ Counts the number of events inside of a well. """
    data = well.get_data()
    count = data.shape[0]
    return count


def calc_channel_median(well, channel='apc'):
    data = well.get_data()
    return data[channels[channel]].median()


def get_file_name_to_well_name_dict(args):
    """
    return the well name in C9 format
    """
    result = {}
    for a in os.listdir(args['folder']):
        if '.fcs' in a:
            result['./'+a] = a.split('_')[2]

    if result == dict():
        print('no .fcs found !!!')
        sys.exit()
    return result


def parse_plate_description_file(file_name):
    """

    """
    result = {'missing_cols': [], 'missing_rows': []}
    for l in open(file_name, 'r'):
        s = l.split()
        if s[0] == 'ligand':
            result[s[0]] = s[1]
        elif s[0] == 'column_strains':
            result['column_strains'] = s[1:]
        elif s[0] == 'conc':
            result['conc'] = [int(a) for a in s[1:]]
        elif s[0] == 'missing_cols':
            result['missing_cols'] = [int(a) for a in s[1:]]
        elif s[0] == 'missing_rows':
            result['missing_rows'] = s[1:]
        else:
            print 'unrecognised row!!!', l
            sys.exit()

    pos_dict = {col: {} for col in result['column_strains']}
    for ni, n in enumerate(range(1, 13)):
        if n in result['missing_cols']:
            continue
        for ai, a in enumerate(list('ABCDEFGH')):
            if a in result['missing_rows']:
                continue
            # pos_dict[result['column_strains'][ni]][result['conc'][ai]] = '%s%i' % (a, n)
            pos_dict['%s%i' % (a, n)] = {'conc': result['conc'][ai],
                                         'column_strain': result['column_strains'][ni]}

    return result, pos_dict

if __name__ == '__main__':
    # print parse_plate_description_file('palte_description.txt')
    main()
