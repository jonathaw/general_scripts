#!/usr/bin/env python3.5
import os
import time
import argparse

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

from scipy.interpolate import griddata
from Logger import Logger
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-log')
    args = vars(parser.parse_args())
    args['logger'] = Logger('logeer_%s.log' % time.strftime("%d.%0-m"))

    z, d, e = parse_rosetta_log(args)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    z1 = np.linspace(min(z), max(z), len(z))
    d1 = np.linspace(min(d), max(d), len(d))

    z2, d2 = np.meshgrid(z1, d1)
    e2 = griddata((z, d), e, (z2, d2), method='cubic')

    surf = ax.plot_surface(z2, d2, e2, rstride=10, cstride=10, cmap=cm.coolwarm, linewidth=0, antialiased=False)

    ax.set_xlabel('Z')
    ax.set_ylabel('D')
    ax.set_zlabel('e')
    # Z, D = np.meshgrid(z, d)
    # args['logger'].log('finished calculating meshgrid...')
    # args['logger'].log('Z is sized %i, D is sized %i' % (len(Z), len(D)))
    # surf = ax.plot_surface(Z, D, e, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)

    plt.show()


def parse_rosetta_log(args):
    args['logger'].log('reading log file, requires lines to have format\nspline z d dielectric')
    args['logger'].log('MAKE SURE TO HAVE THE spline numbers only ONCE in the log file, rosetta prints them multiple times!!!!')
    z = []
    d = []
    e = []

    for l in open(args['log'], 'r'):
        s = l.split()
        if len(s) < 4: continue
        # if s[0] == 'spline' or s[0] == 'spl_':
        if s[0] == 'spl_':
            # if '.' in s[1] or '.' in s[2] and float(s[2]) > 1.4: continue # skip half A z, to quicken things
            z.append(float(s[1]))
            d.append(float(s[2]))
            e.append(float(s[3]))

    args['logger'].log('found %i lines' % len(z))
    args['logger'].log('z values had min %.2f, max %.2f, mean %.2f' % (np.max(z), np.min(z), np.mean(z)))
    args['logger'].log('d values had min %.2f, max %.2f, mean %.2f' % (np.max(d), np.min(d), np.mean(d)))
    args['logger'].log('e values had min %.2f, max %.2f, mean %.2f' % (np.max(e), np.min(e), np.mean(e)))

    return z, d, e


if __name__ == '__main__':
    main()
