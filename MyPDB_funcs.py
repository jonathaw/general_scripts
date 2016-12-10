#!/usr/bin/env python3.5
"""
"""
import copy
from MyPDB import MyPDB, parse_PDB
import numpy as np
import scipy.linalg as linalg
import MyPDB as mp

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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


def find_points_of_closest_distance(l1, l2) -> tuple:
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


def calc_dihedral(p1, p2, p3, p4) -> float:
    p = np.array([p1, p2, p3, p4])
    b = p[:-1] - p[1:]
    b[0] *= -1
    v = np.array([ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]]])
    return np.degrees(np.arccos( v[0].dot(v[1])/(np.linalg.norm(v[0]) * np.linalg.norm(v[1])) ))


def translate_and_rotate_res_to_xy_plane( res: mp.Residue, atom_list: list ) -> mp.Residue:
    xyz = copy.deepcopy(res[ atom_list[0] ].xyz)
    xyz = xyz.scalar_multi( -1 )
    res.translate_xyz( xyz )
    # rotate a1 to xy plane
    proj_a2_xy = copy.deepcopy( res[ atom_list[1] ].xyz )
    proj_a2_xy.z = 0
    a2_copy = copy.deepcopy( res[ atom_list[1] ].xyz )
    ang_a2_xy = np.arccos( a2_copy.unit().dot( proj_a2_xy.unit() ) )
    axis_a2_xy = a2_copy.unit().cross( proj_a2_xy.unit() ).unit()
    rotation_matrix = rotation_matrix_around_vec( axis_a2_xy, ang_a2_xy)
    res.dot_matrix_me( rotation_matrix )

    # rotate a2 to xy plane
    proj_a2_xy = copy.deepcopy( res[ atom_list[2] ].xyz )
    proj_a2_xy.z = 0
    a2_copy = copy.deepcopy( res[ atom_list[2] ].xyz )
    ang_a2_xy = np.arccos( a2_copy.unit().dot( proj_a2_xy.unit()) )
    closest_point = point_on_normed_vec_closest_to_point( proj_a2_xy.as_nparray(), res[atom_list[1]].xyz.unit().as_nparray() )
    ang_a2_xy = angle_between_3_XYZs( a2_copy, closest_point, proj_a2_xy )
    axis = copy.deepcopy( res[ atom_list[1] ].xyz )
    rotation_matrix = rotation_matrix_around_vec( axis, -ang_a2_xy )
    res.dot_matrix_me( rotation_matrix )

    # rotate so that a2 and 3 are on both sides of the Y axis
    a1_copy = res[ atom_list[0] ].xyz
    a2_copy = res[ atom_list[1] ].xyz
    a3_copy = res[ atom_list[2] ].xyz
    ang_312 = angle_between_3_XYZs(a3_copy, a1_copy, a2_copy)
    ang_y12 = angle_between_3_XYZs( mp.XYZ(0, 1, 0), a1_copy, a2_copy )
    axis = mp.XYZ(0, 0, 1)
    rotation_matrix = rotation_matrix_around_vec( axis, -( ang_y12 + 0.5 * ang_312 ) )
    res.dot_matrix_me( rotation_matrix )


def rotation_matrix_around_vec( axis: mp.XYZ, theta: float ):
    """
    calcualted the rotation matrix to roatate an boject around axis theta radian degrees
    based on:
    http://www.euclideanspace.com/maths/algebra/vectors/angleBetween/
    and
    http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    """
    return linalg.expm3(np.cross(np.eye(3), axis.as_list()/np.float64(linalg.norm(axis.as_list()))*theta))


def write_residues_to_file(residues: list, file_name: str) -> None:
    with open(file_name, 'w+') as fout:
        for res in residues:
            for a in res.atoms.values():
                fout.write('%s\n' % a)


def angle_between_3_XYZs( p1: mp.XYZ, p2: mp.XYZ, p3: mp.XYZ ) -> float:
    """
    return the angle between 3 points, in radians:w

    """
    ba = p1.as_nparray() - p2.as_nparray()
    bc = p3.as_nparray() - p2.as_nparray()
    cosine_angle = np.dot( ba, bc ) / ( np.linalg.norm( ba ) * np.linalg.norm( bc ) )
    return np.arccos( cosine_angle )


def point_on_normed_vec_closest_to_point( p: np.array, v: np.array  ) -> mp.XYZ:
    vec = v * np.dot( p, v )
    return mp.XYZ( vec[0], vec[1], vec[2] )


if __name__ == '__main__':
    main()
