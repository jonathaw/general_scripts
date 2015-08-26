#!/usr/bin/env python
def main():
    import matplotlib.pyplot as plt
    import argparse
    import numpy as np
    parser = argparse.ArgumentParser()
    parser.add_argument('-eq')
    parser.add_argument('-rngx', nargs=2, default=[-10, 10])
    parser.add_argument('-rngy', nargs=2, default=[-10, 10])
    args = vars(parser.parse_args())
    for i, t in enumerate(args['rngx']):
        if str(t)[0] == 'm':
            args['rngx'][i] = - int(t[1:])
        else:
            args['rngx'][i] = int(t)
    for i, t in enumerate(args['rngy']):
        if str(t)[0] == 'm':
            args['rngy'][i] = - int(t[1:])
        else:
            args['rngy'][i] = int(t)
    print args
    x = np.array(np.arange(args['rngx'][0], args['rngx'][1], 0.01))
    print x
    res = eval(args['eq'])
    print 'res', res
    plt.xlim(args['rngx'])
    plt.ylim(args['rngy'])
    plt.plot(x, res)
    plt.show()


if __name__ == '__main__':
    main()