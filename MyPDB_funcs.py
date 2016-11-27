#!/usr/bin/env python3.5
"""
"""
from MyPDB import MyPDB, parse_PDB
import numpy as np


def find_helix_vector(pdb: MyPDB, start: int, end: int):
    xs, ys, zs = [], [], []
    for i in range(start, end+1):
       res_i = pdb.get_res(i)
       for bb_atom in res_i.iter_bb():
           xs.append(bb_atom.xyz.x)
           ys.append(bb_atom.xyz.y)
           zs.append(bb_atom.xyz.z)
    xs_ = np.array(xs)
    ys_ = np.array(ys)
    zs_ = np.array(zs)

    dist = np.sqrt((xs[-1]-xs[0])**2 + (ys[-1]-ys[0])**2 + (zs[-1]-zs[0])**2  )

    data = np.concatenate((xs_[:, np.newaxis], ys_[:, np.newaxis], zs_[:, np.newaxis]), axis=1)
    datamean = data.mean(axis=0)
    uu, dd, vv = np.linalg.svd(data - datamean)
    linepts = vv[0] * np.mgrid[-dist/2:dist/2:2j][:, np.newaxis]
    linepts += datamean
    return linepts, data


def find_points_of_closest_distance(l1, l2):
    u = l1[1] - l1[0]
    v = l2[1] - l2[0]
    w = l1[0] - l2[0]
    a = np.dot(u, u)
    b = np.dot(u, v)
    c = np.dot(v, v)
    d = np.dot(u, w)
    e = np.dot(v, w)
    D = (a*c) - (b*b)
    sc, tc = 0.0, 0.0

    if D < 0.000001:
        sc = 0.0
        tc = d / b if b > c else e / c
    else:
        sc = (b*e - c*d) / D
        tc = (a*e - b*d) / D

    # sc = l1[0] + l1*sc
    # tc = l2[0] + l2*tc
    sc = l1[0] + sc * (l1[1] - l1[0])
    tc = l2[0] + tc * (l2[1] - l2[0])
    return sc, tc


def calc_dihedral(p1, p2, p3, p4):
    p = np.array([p1, p2, p3, p4])
    b = p[:-1] - p[1:]
    b[0] *= -1
    v = np.array([ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]]])
    return np.degrees(np.arccos( v[0].dot(v[1])/(np.linalg.norm(v[0]) * np.linalg.norm(v[1])) ))



if __name__ == '__main__':
    main()
