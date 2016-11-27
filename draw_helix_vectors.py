#!/usr/bin/env python3.5
"""
"""
from MyPDB import MyPDB, parse_PDB
import sys
import argparse
import numpy as np
from MyPDB_funcs import find_helix_vector, find_points_of_closest_distance, calc_dihedral
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as m3d


def main():
    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d as m3d
    parser = argparse.ArgumentParser()
    parser.add_argument('-pdb')
    parser.add_argument('-s1', type=int)
    parser.add_argument('-e1', type=int)
    parser.add_argument('-s2', type=int)
    parser.add_argument('-e2', type=int)
    args = vars(parser.parse_args())

    pdb = parse_PDB(args['pdb'])
    l1, data1 = find_helix_vector(pdb, args['s1'], args['e1'])
    l2, data2 = find_helix_vector(pdb, args['s2'], args['e2'])

    p1, p2 = find_points_of_closest_distance(l1, l2)

    cross_angle = calc_dihedral(l1[0], p1, p2, l2[0])
    dist = np.linalg.norm(p1-p2)
    print('cross angle', cross_angle)
    print('dist', dist)
    # draw_helix_scatters_vectors(data1, data2, l1, l2, p1, p2)


def draw_helix_scatters_vectors(data1, data2, l1, l2, p1, p2):
    ax = m3d.Axes3D(plt.figure())
    ax.scatter3D(data1[:, 0], data1[:, 1], data1[:, 2])
    ax.scatter3D(data2[:, 0], data2[:, 1], data2[:, 2])

    ax.scatter3D(p1[0], p1[1], p1[2], s=40, c='r')
    ax.scatter3D(p2[0], p2[1], p2[2], s=20, c='k')

    ax.plot3D(*l1.T)
    ax.plot3D(*l2.T)
    plt.show()


if __name__ == '__main__':
    main()
