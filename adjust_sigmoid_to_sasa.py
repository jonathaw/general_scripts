#!/usr/bin/env python3.5
import scipy
from scipy import optimize
from sklearn.metrics import r2_score

import numpy as np
import matplotlib.pyplot as plt
work_path = '/home/labs/fleishman/jonathaw/elazaridis/sigmoid_sdjustments_7Apr/'
work_path = './'


def main():
    sasa_dict = parse_rsa(work_path+'3o7q.rsa')
    atom_dict = parse_atom_count(work_path+'atom_counts.lst')

    # only use points with SASA < 100 (it's percentages)
    good_by_sasa = [k for k, v in sasa_dict.items() if v < 100.0]
    sasas = np.array([sasa_dict[a] for a in good_by_sasa])
    atoms = np.array([0]+[atom_dict[a] for a in good_by_sasa])

    # resize sasa to 0-1, and add a (0, 1) point at the beginning
    sasas_resized = np.concatenate((np.array([0.95]), resize(sasas)))

    # change the sigma for the first point such that the sigmoid will go through it
    sigma = np.ones(len(atoms))
    sigma[0] = 0.01

    # optimize curve fit
    popt, pcov = scipy.optimize.curve_fit(sigmoid, atoms, sasas_resized, sigma=sigma)

    # calculate R^2
    r2 = r2_score(sasas_resized, sigmoid(atoms, popt[0], popt[1]))
    print('found these parameters:\nintercept: %f\nslope %f' % (popt[1], popt[0]))
    print('R^2 %f' % r2)

    # draw both points and sigmoid
    x = np.linspace(min(atoms), max(atoms), 1000)
    plt.plot(atoms, sasas_resized, '.', x, sigmoid(x, popt[0], popt[1]))
    plt.show()


def sigmoid(x, slope, intercept, a=1.0, b=1.0):
    return 1 / (1 + np.exp((x-intercept) * slope))


def parse_atom_count(file_name):
    result = {}
    for l in open(file_name, 'r'):
        s = l.split()
        if len(s) == 2:
            result[int(s[0])] = int(s[1])
    return result


def parse_rsa(file_name):
    result = {}
    for l in open(file_name, 'r'):
        s = l.split()
        if s[0] == 'RES':
            result[int(s[3])] = float(s[7])
    return result


def resize(arr, lower=0.0, upper=1.0):
    arr=arr.copy()
    if lower>upper: lower,upper=upper,lower
    arr -= arr.min()
    arr *= (upper-lower)/arr.max()
    arr += lower
    return arr


# def calc_sigXsig():


def by_two_numbers():
    from mpl_toolkits.mplot3d import axes3d
    from PDB_Bfactor_by_list import repalce_Bfactor
    from MyPDB import parse_PDB
    by_6A = parse_atom_count(work_path+'6A.lst')
    by_12A = parse_atom_count(work_path+'12A.lst')
    results = {}

    slope6, intercept6 = 0.05, 30
    slope12, intercept12 = 0.05, 350

    for k, v6 in by_6A.items():
        results[k] = sigmoid(v6, slope6, intercept6) * sigmoid(by_12A[k], intercept12, intercept12)
        print(k, v6, by_12A[k], results[k])
    pdb = parse_PDB('3o7q.pdb')
    repalce_Bfactor(pdb, results, {'pdb': '3o7q.pdb', 'list': '6x12sigs.lst'})

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x = np.arange(0, 400, 1)
    y = np.arange(0, 400, 1)
    X, Y = np.meshgrid(x, y)
    zs = np.array([sigmoid(x, slope=slope6, intercept=intercept6)*
                   sigmoid(y, slope=slope12, intercept=intercept12)
                   for x, y in zip(np.ravel(X), np.ravel(Y))])
    Z = zs.reshape(X.shape)

    ax.plot_surface(X, Y, Z)
    ax.set_xlabel('6A')
    ax.set_ylabel('12A')
    ax.set_zlabel('energy')

    plt.show()
    # X, Y, Z = axes3d.get_test_data(0.05)
    # print('x', X)
    # print('y', Y)
    # print('z', Z)


if __name__ == '__main__':
    # main()
    by_two_numbers()