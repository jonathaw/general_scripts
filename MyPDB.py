#!/usr/bin/env python3.5
from List1 import List1
from AASeq import AASeq
import numpy as np

three_2_one = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
               'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
               'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}


class XYZ:
    def __init__(self, x: float = None, y: float = None, z: float = None):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self) -> str:
        return "(%.3f %.3f %.3f)" % (self.x, self.y, self.z)

    def __repr__(self) -> str:
        return "(%f %f %f)" % (self.x, self.y, self.z)

    def __eq__(self, other) -> bool:
        return all([self.x == other.x, self.y == other.y, self.z == other.z])

    def __sub__(self, other):
        """
        :param other:another XYZ instance
        :return:a vector resulting in the subtraction self-other
        >>> a = XYZ(1, 1, 1)
        >>> b = XYZ(2, 2, 2)
        >>> Z = XYZ(-1, -1, -1)
        >>> a-b == Z
        True
        """
        return XYZ(x=self.x - other.x, y=self.y - other.y, z=self.z - other.z)

    def __abs__(self) -> float:
        """
        :return:the magnitude of XYZ
        >>> a = XYZ(1, 1, 1)
        >>> from math import sqrt
        >>> abs(a) == sqrt(3)
        True
        """
        from math import sqrt
        return sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def __add__(self, other):
        return XYZ(x=self.x + other.x, y=self.y + other.y, z=self.z + other.z)

    def cross(self, other):
        """
        :param other: another XYZ instance
        :return:the cross product with self in the left
        >>> x = XYZ(2, 3, 4)
        >>> y = XYZ(5, 6, 7)
        >>> x.cross(y) == XYZ(-3, 6, -3)
        True
        """
        x = self.y * other.z - self.z * other.y
        y = self.z * other.x - self.x * other.z
        z = self.x * other.y - self.y * other.x
        return XYZ(x, y, z)

    def dot(self, other) -> float:
        """
        :param other: another XYZ instance
        :return: the dot product
        >>> a = XYZ(-6, 8, 0)
        >>> b = XYZ(5, 12, 0)
        >>> a.dot(b) == 66
        True
        """
        return self.x * other.x + self.y * other.y + self.z * other.z

    def unit(self):
        """
        :return: the unit vector
        >>> a = XYZ(-2, 1, 0)
        >>> from math import sqrt
        >>> a.unit() == XYZ(-2/sqrt(5), 1/sqrt(5), 0)
        True
        """
        absi = abs(self)
        return XYZ(self.x / absi, self.y / absi, self.z / absi)

    def distance(self, other) -> float:
        """
        :param other:another atom instance
        :return:the euclidian distance
        >>> a = XYZ(x=1, y=1, z=1)
        >>> b = XYZ(x=1, y=1, z=1)
        >>> a.distance(b)
        0.0
        >>> c = XYZ(x=0, y=0, z=0)
        >>> a.distance(c)
        1.7320508075688772
        """
        from math import sqrt
        return sqrt((self.x - other.x) ** 2 + (self.y - other.y) ** 2 + (self.z - other.z) ** 2)


class MembraneResidue:
    def __init__(self):
        self.thkn = XYZ()
        self.cntr = XYZ()
        self.norm = XYZ()
        self.chain = None
        self.res_num = None

    def __repr__(self):
        msg = '\tthkn: %r\n' % self.thkn
        msg += '\tcntr: %r\n' % self.cntr
        msg += '\tnorm: %r\n' % self.norm
        return msg


class Atom:
    def __init__(self, header=None, serial_num=None, name=None, alternate=None, res_type_3=None, chain=None,
                 res_seq_num=None, x=None, y=None, z=None, achar=None, occupancy=None, temp=None, si=None, element=None,
                 charge=None):
        self.header = header
        self.serial_num = serial_num
        self.name = name
        self.alternate = alternate
        self.res_type_3 = res_type_3
        self.chain = chain
        self.res_seq_num = res_seq_num
        # self.x = x
        # self.y = y
        # self.z = z
        self.achar = achar
        self.occupancy = occupancy
        self.temp = temp
        self.si = si
        self.element = element
        self.charge = charge
        self.xyz = XYZ(x=x, y=y, z=z)

    def __repr__(self) -> str:
        return "%-6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n" % \
               (self.header, self.serial_num, self.name, self.achar, self.res_type_3, self.chain, self.res_seq_num,
                self.si, self.xyz.x, self.xyz.y, self.xyz.z, self.occupancy, self.temp, self.element, self.charge)

    def __str__(self) -> str:
        return "%-6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s" % \
               (self.header, self.serial_num, self.name, self.achar, self.res_type_3, self.chain, self.res_seq_num,
                self.si, self.xyz.x, self.xyz.y, self.xyz.z, self.occupancy, self.temp, self.element, self.charge)

    def __cmp__(self, other) -> bool:
        return self.serial_num > other.serial_num

    def __gt__(self, other) -> bool:
        return self.serial_num > other.serial_num

    def __eq__(self, other) -> bool:
        return self.serial_num == other.serial_num

    def __ge__(self, other) -> bool:
        return self.serial_num >= other.serial_num

    def set_temp(self, temp: float) -> None:
        """
        B factor
        """
        self.temp = temp

    def set_occupancy(self, occupancy) -> None:
        self.occupancy = occupancy

    def distance(self, other) -> float:
        """
        :param other:another atom instance
        :return:the euclidian distance
        >>> a = Atom(x=1, y=1, z=1)
        >>> b = Atom(x=1, y=1, z=1)
        >>> a.distance(b)
        0.0
        >>> c = Atom(x=0, y=0, z=0)
        >>> a.distance(c)
        1.7320508075688772
        """
        from math import sqrt
        return sqrt((self.xyz.x - other.xyz.x) ** 2 + (self.xyz.y - other.xyz.y) ** 2 + (self.xyz.z - other.xyz.z) ** 2)

    def change_chain_name(self, new: str) -> None:
        self.chain = new

    def translate_xyz(self, xyz: XYZ) -> None:
        """
        :param xyz: XYZ point
        :return: None. translate atom by x, y, z
        """
        self.xyz += xyz


class Residue:
    def __init__(self, res_type_3: str = None, res_num: int = None, chain: str = None, atoms: dict = None):
        self.res_type_3 = res_type_3
        if res_type_3 in three_2_one.keys():
            self.res_type = three_2_one[res_type_3]
        else:
            self.res_type = res_type_3
        self.res_num = res_num
        self.chain = chain
        if atoms is None:
            self.atoms = {}
        else:
            self.atoms = atoms
        self.memb_z = None

    def __repr__(self) -> str:
        if self.memb_z is None:
            msg = 'Chain %s Res #%i Type %s with %i atoms' % (
                self.chain, self.res_num, self.res_type_3, len(self.atoms))
        else:
            msg = 'Chain %s Res #%i Type %s with %i atoms, memb Z %.2f' % (
                self.chain, self.res_num, self.res_type_3, len(self.atoms), self.memb_z)
        return msg

    def __getitem__(self, item: str) -> Atom:
        return self.atoms[item]

    def __iter__(self):
        for k, v in self.atoms.items():
            yield k, v

    def remove_atom(self, atom: Atom) -> None:
        new_res = {}
        for aid, a in self:
            if a != atom:
                # print(a)
                new_res[aid] = a
        self.atoms = new_res

    def keys(self):
        return self.atoms.keys()

    def add_atom(self, atom: Atom) -> None:
        self.atoms[atom.name] = atom

    def min_distance_res(self, other) -> float:
        """
        :param other:another residue
        :return:the minimu, distance between all atoms in residues
        >>> a = Residue(atoms={1: Atom(x=0, y=0, z=0), 2: Atom(x=1, y=1, z=1)})
        >>> b = Residue(atoms={1: Atom(x=2, y=2, z=2), 2: Atom(x=3, y=3, z=3)})
        >>> a.min_distance_res(b)
        1.7320508075688772
        """
        dists = []
        for mid, m in self:
            for oid, o in other:
                dists.append(m.distance(o))
        return min(dists)

    def phi(self, prev_res) -> float:
        return dihedral(prev_res['C'].xyz, self['N'].xyz, self['CA'].xyz, self['C'].xyz)

    def psi(self, next_res) -> float:
        return dihedral(self['N'].xyz, self['CA'].xyz, self['C'].xyz, next_res['N'].xyz)

    def min_dist_xyz(self, xyz: XYZ):
        try:
            return distance_two_point_to_point(self['CA'].xyz, self['CB'].xyz, xyz.xyz)
        except:
            print('failed at', self, xyz)
            return None

    def change_chain_name(self, new: str) -> None:
        self.chain = new
        for aid, a in self:
            a.change_chain_name(new)

    def translate_xyz(self, xyz: XYZ) -> None:
        """
        :param xyz: xyz point
        :return: None. translates all residue atoms by x, y, z
        """
        for aid, a in self:
            a.translate_xyz(xyz)


class Chain:
    def __init__(self, chain_id: str = None, residues: dict = None, non_residues: dict = None):
        self.chain_id = chain_id
        self.residues = residues if residues is not None else {}
        self.seq = AASeq(''.join(a.res_type for a in residues.values())) if residues is not None else \
            AASeq('', name=chain_id)
        self.non_residues = non_residues if non_residues is not None else {}
        self.non_residues_seq = AASeq(''.join(a.res_type for a in residues.values()), name=chain_id) if \
            non_residues is not None else AASeq('', name=chain_id)

    def __repr__(self) -> str:
        return "chain %s has %i residues" % (self.chain_id, len(self.residues))

    def __getitem__(self, item: int) -> Residue:
        try:
            return self.residues[item]
        except:
            return self.non_residues[item]

    def __iter__(self):
        for k, v in self.residues.items():
            yield k, v

    def add_residue(self, residue: Residue) -> None:
        if residue.res_type_3 in three_2_one.keys():
            self.seq.add_aa(residue.res_type)
            self.residues[residue.res_num] = residue
        else:
            self.non_residues_seq.add_aa(residue.res_type)
            self.non_residues[residue.res_num] = residue

    def min_distance_chain(self, other: Residue) -> float:
        distances = []
        for mrid, mres in self:
            for orid, ores in other:
                distances.append(mres.min_distance_res(ores))
        return min(distances)

    def keys(self):
        return self.residues.keys()

    def values(self):
        return self.residues.values()

    def COM(self) -> XYZ:
        """
        :return:the Center Of Mass of the chain as calculated by the averages over Xs, Ys and Zs of all CAs
        """
        Xs = []
        Ys = []
        Zs = []
        for res in self.values():
            if 'CA' in res.keys():
                Xs.append(res['CA'].xyz.x)
                Ys.append(res['CA'].xyz.y)
                Zs.append(res['CA'].xyz.z)
        return XYZ(np.mean(Xs), np.mean(Ys), np.mean(Zs))

    def change_chain_name(self, new: str) -> None:
        self.chain_id = new
        for rid, r in self:
            r.change_chain_name(new)

    def translate_xyz(self, xyz: XYZ) -> None:
        """
        :param xyz: an xyz point
        :return: None. translate all chain atoms by xyz
        """
        for rid, r in self:
            r.translate_xyz(xyz)


class MyPDB:
    def __init__(self, name: str = None, chains: dict = None, seqs: dict = None):
        self.name = name
        self.chains = chains if chains is not None else {}
        self.seqs = seqs if seqs is not None else {}
        self.memb_res = None

    @property
    def __repr__(self) -> str:
        msg = u"PDB {0:s} has {1:d} chains".format(self.name, len(self.chains))
        for c in self.chains:
            msg += repr(self.chains[c])
        return msg

    def __getitem__(self, item: str) -> Chain:
        return self.chains[item]

    def __iter__(self):
        for k, v in self.chains.items():
            yield k, v

    def add_chain(self, chain: Chain) -> None:
        """
        :param chain: a Chain
        :return: appends chain to PDB
        """
        self.chains[chain.chain_id] = chain

    def add_atom(self, atom: Atom) -> None:
        if atom.chain not in self.chains.keys():
            self.add_chain(chain=Chain(chain_id=atom.chain))
        if atom.res_seq_num not in self.chains[atom.chain].residues.keys():
            self[atom.chain].add_residue(residue=Residue(res_type_3=atom.res_type_3, res_num=atom.res_seq_num,
                                                         chain=atom.chain))
        self[atom.chain][atom.res_seq_num].add_atom(atom=atom)
        if atom.res_type_3 in three_2_one.keys():
            self.seqs[atom.chain] = self[atom.chain].seq
        else:
            self.seqs[atom.chain + '_non_res'] = self[atom.chain].non_residues_seq

    def change_chain_name(self, old: str, new: str) -> None:
        self.chains[old].change_chain_name(new)

    def renumber(self) -> None:
        """
        :return: renumbers self
        """
        i = 1
        for cid, c in self:
            for rid, r in c:
                for aid, a in r:
                    self[cid][rid][aid].serial_num = i
                    i += 1

    def remove_hydrogens(self) -> None:
        """
        removes all Hydrogen atoms from instance
        """
        for cid, c in self:
            for rid, r in c:
                for aid, a in r:
                    if a.element == 'H':
                        print('removing H at %s' % aid)
                        r.remove_atom(a)

    def translate_xyz(self, xyz: XYZ) -> None:
        """
        :param xyz: a point
        :return: None. translates all pdb points by x, y and z
        """
        for cid, c in self:
            c.translate_xyz(xyz)

    def add_memb_res(self, memb_res: MembraneResidue) -> None:
        self.memb_res = memb_res

        # go over all residues, assign membrane Z value
        for cid in sorted(self.chains.keys()):
            for rid, res in sorted(self[cid].residues.items()):
                if -15. <= res['CA'].xyz.x <= 15:
                    res.memb_z = res['CA'].xyz.x
                else:
                    res.memb_z = None

    def summarize(self):
        print('MyPDB instance with:')
        print('\t%i chains' % len(self.chains))
        print('\tsequences %s' % '\n\t'.join('>%s\n%s' % (k, v) for k, v in self.seqs.items()))
        print('\tmembrnae residue \n%r' % self.memb_res)


def parse_membrane_residue(pdb_lines: list) -> MembraneResidue:
    """
    pdb_lines: a list of text lines from .pdb file
    returns a MembraneResidue instance
    """
    result = MembraneResidue()
    for l in pdb_lines:
        s = l.split()
        if s != 0:
            if s[2] == 'THKN':
                result.thkn = XYZ(x=float(s[6]), y=float(s[7]), z=float(s[8]))
            elif s[2] == 'CNTR':
                result.cntr = XYZ(x=float(s[6]), y=float(s[7]), z=float(s[8]))
            elif s[2] == 'NORM':
                result.norm = XYZ(x=float(s[6]), y=float(s[7]), z=float(s[8]))
                result.chain = s[4]
                result.res_num = int(s[5])
    return result


def distance_two_point_to_point(p1: XYZ, p2: XYZ, x: XYZ) -> float:
    """
    :param p1:point 1 (XYZ) on line
    :param p2:point 2 (XYZ) on line
    :param x:point x (XYZ)
    :return: the minimal distance between the line pw-p1 and point x
    >>> p1 = XYZ(0, 0, 0)
    >>> p2 = XYZ(1, 0, 0)
    >>> x = XYZ(1, 1, 0)
    >>> distance_two_point_to_point(p1, p2, x)
    1.0
    >>> x = XYZ(-1, -1, 0)
    >>> distance_two_point_to_point(p1, p2, x)
    1.0
    >>> p2 = XYZ(-1, 0, 0)
    >>> distance_two_point_to_point(p1, p2, x)
    1.0
    """
    Ul = (p2 - p1).unit()  # the unit vector between p2 to p1
    w = x - p1  # the vector from p1 to x
    return abs(Ul.cross(w))


def is_point_infront_points_vec(p1: XYZ, p2: XYZ, x: XYZ) -> float:
    """
    :param p1:point 1 in vec
    :param p2:point 2 in vec
    :param x: a point in space
    :return: True iff point x is "in front" of the direction of vector p2-p1
    >>> p1 = XYZ(0, 0, 0)
    >>> p2 = XYZ(1, 0, 0)
    >>> x = XYZ(1, 1, 0)
    >>> is_point_infront_points_vec(p1, p2, x)
    True
    >>> x = XYZ(-1, -1, 0)
    >>> is_point_infront_points_vec(p1, p2, x)
    False
    """
    oreintation_vec = p2 - p1
    inves_vec = x - p1
    return oreintation_vec.dot(inves_vec) > 0


def dihedral(p0: XYZ, p1: XYZ, p2: XYZ, p3: XYZ) -> float:
    """
    used http://www.cgl.ucsf.edu/Outreach/pc204/lecture_notes/phipsi/structured/phipsi.py
    :param p0: XYZ
    :param p1: XYZ
    :param p2: XYZ
    :param p3: XYZ
    :return: the dihedral angle between p0-3 in degrees
    >>> p0 =
    """
    from math import acos, degrees
    v01 = p0 - p1
    v32 = p3 - p2
    v12 = p1 - p2
    v0 = v12.cross(v01)
    v3 = v12.cross(v32)
    angle = degrees(acos(v0.dot(v3) / abs(v0) / abs(v3)))
    if v0.cross(v3).dot(v12) > 0:
        angle = -angle
    return angle


def extract_seq(pdb: MyPDB) -> dict:
    seqs = {}
    for cid, c in pdb:
        seqs[cid] = AASeq(name='%s.%s' % (pdb.name, cid))
        seq = ''
        for rid, r in c:
            seq += r.res_type
        seqs[cid].set_seq(seq)
    return seqs


def parse_PDB(file_in: str, name: str = None, with_non_residue: bool = True) -> MyPDB:
    if file_in[-3:] == '.gz':
        import gzip
        fin = gzip.open(file_in, 'rb')
        cont = fin.read().decode('utf-8').split('\n')
    else:
        fin = open(file_in, 'r')
        cont = fin.read().split('\n')
    pdb = MyPDB(name=name)
    memb_res = []
    for l in cont:
        s = l.split()
        if len(s) < 1:
            continue
        if s[0] in ['ATOM', 'HETATM']:
            if not with_non_residue:
                if l[17:20] not in three_2_one.keys():
                    continue
            if 'H' in s[2] and ('1' in s[2] or '2' in s[2] or '3' in s[2]):
                continue
            if s[3] == 'MEM':
                memb_res.append(l)
                continue
            atom = Atom(header=s[0], serial_num=int(l[6:11]), name=l[12:16].replace(' ', ''),
                        alternate=l[16] if l[16] != ' ' else None,
                        res_type_3=l[17:20], chain=l[21].upper(), res_seq_num=int(l[22:26]), achar=l[26],
                        x=float(l[30:38]),
                        y=float(l[38:46]), z=float(l[47:55]), occupancy=float(l[55:61]), temp=float(l[60:66]),
                        si=l[72:76].replace(' ', ''), element=l[76:78].replace(' ', ''),
                        charge=str(l[78:80].replace(' ', '')))
            if atom.alternate is 'B':
                continue
            pdb.add_atom(atom)
    pdb.renumber()
    if memb_res != []:
        mm_res = parse_membrane_residue(memb_res)
        pdb.add_memb_res(mm_res)
    return pdb


def write_PDB(file_out: str, pdb: MyPDB) -> None:
    atoms = []
    for cid in sorted(pdb.chains.keys()):
        for rid in sorted(pdb[cid].residues.keys()):
            for aid in sorted(pdb[cid][rid].keys()):
                atoms.append(pdb[cid][rid][aid])
    with open(file_out, 'w') as fout:
        for a in sorted(atoms):
            fout.write(str(a) + '\n')


def extract_chain(file_out: str, pdb: MyPDB, chain: str = 'A') -> None:
    with open(file_out, 'wr+') as fout:
        for rid, r in pdb[chain]:
            for aid, a in r:
                fout.write(str(a) + '\n')


def draw_ramachadran(pdb: MyPDB) -> None:
    import matplotlib.pyplot as plt
    phis = {}
    psis = {}
    for cid, c in pdb:
        for rid, r in c:
            try:
                prev_res = c[rid - 1]
                phis[rid] = r.phi(prev_res)
            except:
                pass
            try:
                next_res = c[rid + 1]
                psis[rid] = r.psi(next_res)
            except:
                pass
    plt.scatter(list(phis.values()), list(psis.values()), alpha=0.5)
    plt.xlim((-180., 180))
    plt.ylim((-180., 180))
    plt.xlabel('Phi')
    plt.ylabel('Psi')
    plt.show()


def interface_residues(ch1: Chain, ch2: Chain, dist: float = 10.0) -> list:
    """
    :type ch1: Chain
    :type ch2: Chain
    :param ch1: chain of interest
    :param ch2: other chain in possible interface
    :param dist: maxiaml dist for two residues to be declared as "close"
    :return: a list of Residue instances from ch1 that are close to ch2, point to it's direction,
    and do not point to ch1's COM
    """
    residues = []
    ch1_com = ch1.COM()
    # print('ch1COM', ch1_com)
    for r1 in ch1.values():
        for r2 in ch2.values():
            if 'CB' in r1.keys() and 'CA' in r1.keys() and 'CA' in r2.keys():
                if is_point_infront_points_vec(r1['CA'].xyz, r1['CB'].xyz, r2['CA'].xyz) and \
                        not is_point_infront_points_vec(r1['CA'].xyz, r1['CB'].xyz, ch1_com) and \
                                r1.min_dist_xyz(r2['CA']) < dist:
                    residues.append(r1)
    return list(set(residues))


def com_residues(chain: Chain, residues: list) -> XYZ:
    """
    :param residues: list of residue numbers
    :return: XYZ describing the COM
    """
    Xs, Ys, Zs = [], [], []
    for res in residues:
        resi = chain[res]
        if 'CA' in resi.keys():
            Xs.append(resi['CA'].xyz.x)
            Ys.append(resi['CA'].xyz.y)
            Zs.append(resi['CA'].xyz.z)
    return XYZ(np.mean(Xs), np.mean(Ys), np.mean(Zs))


def memb_residues(pdb: MyPDB) -> set():
    """
    collect a set of residues with memb_z within [-15, 15]
    """
    result = []
    for ch in pdb.chains.values():
        for res in ch.values():
            if res.memb_z is not None:
                result.append(res)
    return result

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-mode')
    parser.add_argument('-in_file')
    parser.add_argument('-out_file', default='tmp.pdb')
    parser.add_argument('-chains', nargs='+', default=['A'])
    parser.add_argument('-name', default=None)
    parser.add_argument('-no_non_residues', default=True, type=bool)
    parser.add_argument('-old_chain')
    parser.add_argument('-new_chain')
    parser.add_argument('-seq_positions', type=int, nargs='+', help='sequence positions to extract', default=False)
    args = vars(parser.parse_args())

    if args['name'] is None:
        args['name'] = args['in_file'].split('.')[0].split('_')[0]

    if args['mode'] == 'extract_seq':
        pdb = parse_PDB(args['in_file'], args['name'])
        if not args['seq_positions']:
            for k, v in pdb.seqs.items():
                if k in args['chains']:
                    print('>%s.%s' % (args['name'], k))
                    print(v.get_seq)
        else:
            for k, v in pdb.seqs.items():
                if k in args['chains']:
                    print('>%s.%s' % (args['name'], k))
                    print('%s' % ''.join(v.get_positions(args['seq_positions'])))

    elif args['mode'] == 'extract_chain':
        pdb = parse_PDB(args['in_file'], args['name'])
        t_pdb = MyPDB(name=args['name'])
        for chain in args['chains']:
            ch = pdb.chains[chain.upper()]
            t_pdb.add_chain(ch)
        if args['out_file'] is None:
            args['out_file'] = args['in_file'].replace('.pdb', '_%s.pdb' % ''.join(args['chains']))
        write_PDB(args['out_file'], t_pdb)

    elif args['mode'] == 'copy':
        pdb = parse_PDB(args['in_file'], args['name'], with_non_residue=False)

    elif args['mode'] == 'min_distances':
        pdb = parse_PDB(args['in_file'], args['name'])
        while args['chains']:
            a = args['chains'].pop()
            while args['chains']:
                b = args['chains'].pop()
                print("The minimal distance between chains %s and %s is %f" % (a, b, pdb[a].min_distance_chain(pdb[b])))
    elif args['mode'] == 'ramachadran':
        pdb = parse_PDB(args['in_file'], args['name'])
        draw_ramachadran(pdb)

    elif args['mode'] == 'interface':
        pdb = parse_PDB(args['in_file'], args['name'])
        inter_0 = interface_residues(pdb[args['chains'][0]], pdb[args['chains'][1]], with_dir=True)
        inter_1 = interface_residues(pdb[args['chains'][1]], pdb[args['chains'][0]], with_dir=True)
        ch1 = set([a.res_num for a in inter_0])
        ch2 = set([a.res_num for a in inter_1])
        print("select ch_%s_inter, %s and res %s" % (
        args['chains'][0], args['in_file'][:-4], '+'.join([str(a) for a in ch1])))
        print("select ch_%s_inter, %s and res %s" % (
        args['chains'][1], args['in_file'][:-4], '+'.join([str(a) for a in ch2])))

    elif args['mode'] == 'ramachandran':
        pdb = parse_PDB('/home/labs/fleishman/jonathaw/scripts/general_scripts/test_1ohz_AB.pdb', args['name'])
        print(pdb.seqs)
        draw_ramachadran(pdb)

    elif args['mode'] == 'change_chain':
        pdb = parse_PDB(args['in_file'], args['name'])
        pdb.change_chain_name(args['old_chain'], args['new_chain'])
        write_PDB(args['out_file'], pdb)

    elif args['mode'] == 'remove_h':
        pdb = parse_PDB(args['in_file'], args['name'])
        pdb.remove_hydrogens()
        write_PDB(args['out_file'], pdb)

    elif args['mode'] == 'test':
        pdb = parse_PDB(args['in_file'], args['name'])
        pdb.summarize()

    else:
        print("mode not found")
