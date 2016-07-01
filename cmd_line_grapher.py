#!/usr/bin/env python3.5
import argparse
import numpy as np
import matplotlib.pyplot as plt
from Equation import Expression


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-type', default='hist', help='graph type to use. suitable options are hist, plot and scatter '
                                                      'and box')
    parser.add_argument('-x', nargs='+', type=float)
    parser.add_argument('-y', nargs='+', type=float)
    parser.add_argument('-x_lim', nargs='+', type=float, default=None)
    parser.add_argument('-y_lim', nargs='+', type=float, default=None)
    parser.add_argument('-eq', type=str)

    args = vars(parser.parse_args())

    if args['type'] == 'hist':
        draw_hist(args)

    elif args['type'] == 'scatter':
        draw_scatter(args)

    elif args['type'] == 'plot':
        draw_plot(args)

    elif args['type'] == 'box':
        draw_boxplot(args)

    elif args['type'] == 'eq':
        draw_equation(args)

    else:
        print('no type given')


def draw_equation(args):
    eq = Expression(args['eq'])
    rng = np.arange(args['x_lim'][0], args['x_lim'][1], (args['x_lim'][1]-args['x_lim'][0])/1000.0)
    plt.plot(rng, [eq(a) for a in rng])

    if args['x_lim'] is not None:
        plt.xlim(args['x_lim'])
    if args['y_lim'] is not None:
        plt.xlim(args['y_lim'])

    plt.show()


def draw_boxplot(args) -> None:
    plt.boxplot(args['x'])
    plt.show()


def draw_scatter(args) -> None:
    plt.scatter(args['x'], args['y'])

    if args['x_lim'] is not None:
        plt.xlim(args['x_lim'])
    if args['y_lim'] is not None:
        plt.xlim(args['y_lim'])

    plt.show()


def draw_plot(args) -> None:
    plt.scatter(args['x'], args['y'])

    if args['x_lim'] is not None:
        plt.xlim(args['x_lim'])
    if args['y_lim'] is not None:
        plt.xlim(args['y_lim'])

    plt.show()


def draw_hist(args: dict) -> None:
    plt.hist(args['x'])

    if args['x_lim'] is not None:
        plt.xlim(args['x_lim'])
    if args['y_lim'] is not None:
        plt.xlim(args['y_lim'])

    plt.show()

if __name__ == '__main__':
    main()
