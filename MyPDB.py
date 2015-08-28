three_2_one = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
               'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
               'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}


class MyPDB():
    def __init__(self, name=None, chains=None, seqs=None):
        self.name = name
        self.chains = chains if chains is not None else {}
        self.seqs = seqs if seqs is not None else {}

    def __repr__(self):
        msg = "PDB %s has %i chains" % (self.name, len(self.chains.keys()))
        for c in self.chains:
            msg += repr(self.chains[c])
        return msg

    def __getitem__(self, item):
        return self.chains[item]

    def __iter__(self):
        for k, v in self.chains.items():
            yield k, v

    def add_chain(self, chain):
        self.chains[chain.chain_id] = chain

    def add_atom(self, atom):
        if atom.chain not in self.chains.keys():
            self.add_chain(chain=Chain(chain_id=atom.chain))
        if atom.res_seq_num not in self.chains[atom.chain].residues.keys():
            self[atom.chain].add_residue(residue=Residue(res_type_3=atom.res_type_3, res_num=atom.res_seq_num,
                                                         chain=atom.chain))
        self[atom.chain][atom.res_seq_num].add_atom(atom=atom)
        self.seqs[atom.chain] = self[atom.chain].seq


class Chain():
    def __init__(self, chain_id=None, residues=None, seq=''):
        self.chain_id = chain_id
        self.residues = residues if residues is not None else {}
        self.seq = seq

    def __repr__(self):
        return "chain %s has %i residues" % (self.chain_id, len(self.residues.keys()))

    def __getitem__(self, item):
        return self.residues[item]

    def __iter__(self):
        for k, v in self.residues.items():
            yield k, v

    def add_residue(self, residue):
        self.seq += residue.res_type
        self.residues[residue.res_num] = residue

    def min_distance_chain(self, other):
        distances = []
        for mrid, mres in self:
            for orid, ores in other:
                distances.append(mres.min_distance_res(ores))
        return min(distances)

    def keys(self):
        return self.residues.keys()

    def values(self):
        return self.residues.values()

    def COM(self):
        """
        :return:the Center Of Mass of the chain as calculated by the averages over Xs, Ys and Zs of all CAs
        """
        from numpy import mean
        Xs = []
        Ys = []
        Zs = []
        for res in self.values():
            if 'CA' in res.keys():
                Xs.append(res['CA'].xyz.x)
                Ys.append(res['CA'].xyz.y)
                Zs.append(res['CA'].xyz.z)
        return XYZ(mean(Xs), mean(Ys), mean(Zs))


class Residue():
    def __init__(self, res_type_3=None, res_num=None, chain=None, atoms=None):
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

    def __repr__(self):
        msg = 'Chain %s Res #%i Type %s with %i atoms' % (self.chain, self.res_num, self.res_type_3, len(self.atoms.keys()))
        return msg

    def __getitem__(self, item):
        return self.atoms[item]

    def __iter__(self):
        for k, v in self.atoms.items():
            yield k, v

    def keys(self):
        return self.atoms.keys()

    def add_atom(self, atom):
        self.atoms[atom.name] = atom

    def min_distance_res(self, other):
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

    def phi(self, prev_res):
        return dihedral(prev_res['C'].xyz, self['N'].xyz, self['CA'].xyz, self['C'].xyz)

    def psi(self, next_res):
        return dihedral(self['N'].xyz, self['CA'].xyz, self['C'].xyz, next_res['N'].xyz)

    def min_dist_xyz(self, xyz):
        try:
            return distance_two_point_to_point(self['CA'].xyz, self['CB'].xyz, xyz.xyz)
        except:
            print 'failed at', self, xyz
            return None


class Atom():
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
        self.x = x
        self.y = y
        self.z = z
        self.achar = achar
        self.occupancy = occupancy
        self.temp = temp
        self.si = si
        self.element = element
        self.charge = charge
        self.xyz = XYZ(x=x, y=y, z=z)

    def __repr__(self):
        return "%-6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n" % \
               (self.header, self.serial_num, self.name, self.achar, self.res_type_3, self.chain, self.res_seq_num,
                self.si, self.x, self.y, self.z, self.occupancy, self.temp, self.element, self.charge)

    def __str__(self):
            return "%-6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s" % \
                   (self.header, self.serial_num, self.name, self.achar, self.res_type_3, self.chain, self.res_seq_num,
                    self.si, self.x, self.y, self.z, self.occupancy, self.temp, self.element, self.charge)

    def distance(self, other):
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
        return sqrt((self.x-other.x)**2+(self.y-other.y)**2+(self.z-other.z)**2)


class XYZ():
    def __init__(self, x=None, y=None, z=None):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return "(%.3f %.3f %.3f)" % (self.x, self.y, self.z)

    def __repr__(self):
        return "(%f %f %f)" % (self.x, self.y, self.z)

    def __eq__(self, other):
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
        return XYZ(x=self.x-other.x, y=self.y-other.y, z=self.z-other.z)

    def __abs__(self):
        """
        :return:the magnitude of XYZ
        >>> a = XYZ(1, 1, 1)
        >>> from math import sqrt
        >>> abs(a) == sqrt(3)
        True
        """
        from math import sqrt
        return sqrt(self.x**2 + self.y**2 + self.z**2)

    def __add__(self, other):
        return XYZ(x=self.x+other.x, y=self.y+other.y, z=self.z+other.z)

    def cross(self, other):
        """
        :param other: another XYZ instance
        :return:the cross product with self in the left
        >>> x = XYZ(2, 3, 4)
        >>> y = XYZ(5, 6, 7)
        >>> x.cross(y) == XYZ(-3, 6, -3)
        True
        """
        x = self.y*other.z - self.z*other.y
        y = self.z*other.x - self.x*other.z
        z = self.x*other.y - self.y*other.x
        return XYZ(x, y, z)

    def dot(self, other):
        """
        :param other: another XYZ instance
        :return: the dot product
        >>> a = XYZ(-6, 8, 0)
        >>> b = XYZ(5, 12, 0)
        >>> a.dot(b) == 66
        True
        """
        return self.x*other.x + self.y*other.y + self.z*other.z

    def unit(self):
        """
        :return: the unit vector
        >>> a = XYZ(-2, 1, 0)
        >>> from math import sqrt
        >>> a.unit() == XYZ(-2/sqrt(5), 1/sqrt(5), 0)
        True
        """
        absi = abs(self)
        return XYZ(self.x/absi, self.y/absi, self.z/absi)

    def distance(self, other):
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
        return sqrt((self.x-other.x)**2 + (self.y-other.y)**2 + (self.z-other.z)**2)


def distance_two_point_to_point(p1, p2, x):
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


def is_point_infront_points_vec(p1, p2, x):
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


def dihedral(p0, p1, p2, p3):
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
    angle = degrees(acos(v0.dot(v3)/abs(v0)/abs(v3)))
    if v0.cross(v3).dot(v12) > 0:
        angle = -angle
    return angle


def extract_seq(pdb):
    print pdb
    seqs = {}
    for cid, c in pdb:
        print cid
        seqs[cid] = ''
        for rid, r in c:
            seqs[cid] += r.res_type_3
    return seqs


def parse_PDB(file_in, name=None):
    with open(file_in, 'r') as f:
        cont = f.read().split('\n')
    pdb = MyPDB(name=name)
    for l in cont:
        s = l.split()
        if len(s) < 1:
            continue
        if s[0] in ['ATOM', 'HETATM']:
            atom = Atom(header=s[0], serial_num=int(l[6:11]), name=l[12:16].replace(' ', ''), alternate=l[16]if l[16] != ' ' else None,
                        res_type_3=l[17:20], chain=l[21], res_seq_num=int(l[22:26]), achar=l[26], x=float(l[30:38]),
                        y=float(l[38:46]), z=float(l[47:55]), occupancy=float(l[55:61]), temp=float(l[60:66]),
                        si=l[72:76].replace(' ', ''), element=l[76:78].replace(' ', ''), charge=str(l[78:80].replace(' ', '')))
            if atom.alternate is not None:
                continue
            pdb.add_atom(atom)
    return pdb


def write_PDB(file_out, pdb):
    with open(file_out, 'wr+') as fout:
        for ch_id, chain in pdb:
            for res_id, res in chain:
                for a_id, a in res:
                    fout.write(str(a) + '\n')


def extract_chain(file_out, pdb, chain='A'):
    with open(file_out, 'wr+') as fout:
        for rid, r in pdb[chain]:
            for aid, a in r:
                fout.write(str(a) + '\n')


def draw_ramachadran(pdb):
    import matplotlib.pyplot as plt
    phis = {}
    psis = {}
    for cid, c in pdb:
        for rid, r in c:
            try:
                prev_res = c[rid-1]
                phis[rid] = r.phi(prev_res)
            except:
                pass
            try:
                next_res = c[rid+1]
                psis[rid] = r.psi(next_res)
            except:
                pass
    plt.scatter(phis.values(), psis.values())
    plt.xlim((-180., 180))
    plt.ylim((-180., 180))
    plt.show()


def interface_residues(ch1, ch2, dist=10.0):
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
    print 'ch1COM', ch1_com
    for r1 in ch1.values():
        for r2 in ch2.values():
            if 'CB' in r1.keys() and 'CA' in r1.keys() and 'CA' in r2.keys():
                if is_point_infront_points_vec(r1['CA'].xyz, r1['CB'].xyz, r2['CA'].xyz) and \
                        not is_point_infront_points_vec(r1['CA'].xyz, r1['CB'].xyz, ch1_com) and \
                                r1.min_dist_xyz(r2['CA']) < dist:
                    residues.append(r1)
    return residues


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode')
    parser.add_argument('-in_file')
    parser.add_argument('-out_file', default=None)
    parser.add_argument('-chains', nargs='+', default=['A'])
    parser.add_argument('-name', default=None)
    args = vars(parser.parse_args())
    if args['name'] is None:
        args['name'] = args['in_file'].split('.')[0].split('_')[0]
    if args['mode'] == 'extract_seq':
        pdb = parse_PDB(args['in_file'], args['name'])
        for k, v in pdb.seqs.items():
            print '>%s.%s' % (args['name'], k)
            print v
    elif args['mode'] == 'extract_chain':
        pdb = parse_PDB(args['in_file'], args['name'])
        t_pdb = MyPDB(name=args['name'])
        for chain in args['chains']:
            ch = pdb.chains[chain.upper()]
            t_pdb.add_chain(ch)
        if args['out_file'] is None:
            args['out_file'] =args['in_file'].replace('.pdb', '_%s.pdb' % ''.join(args['chains']))
        write_PDB(args['out_file'], t_pdb)
    elif args['mode'] == 'min_distances':
        pdb = parse_PDB(args['in_file'], args['name'])
        while args['chains']:
            a = args['chains'].pop()
            while args['chains']:
                b = args['chains'].pop()
                print "The minimal distance between chains %s and %s is %f" % (a, b, pdb[a].min_distance_chain(pdb[b]))
    elif args['mode'] == 'ramachadran':
        pdb = parse_PDB(args['in_file'], args['name'])
        draw_ramachadran(pdb)
    elif args['mode'] == 'interface':
        pdb = parse_PDB(args['in_file'], args['name'])
        inter_0 = interface_residues(pdb[args['chains'][0]], pdb[args['chains'][1]], with_dir=True)
        inter_1 = interface_residues(pdb[args['chains'][1]], pdb[args['chains'][0]], with_dir=True)
        ch1 = set([a.res_num for a in inter_0])
        ch2 = set([a.res_num for a in inter_1])
        print "select ch_%s_inter, %s and res %s" % (args['chains'][0], args['in_file'][:-4], '+'.join([str(a) for a in ch1]))
        print "select ch_%s_inter, %s and res %s" % (args['chains'][1], args['in_file'][:-4], '+'.join([str(a) for a in ch2]))
    elif args['mode'] == 'test':
        pdb = parse_PDB(args['in_file'], args['name'])
        print pdb['A'][79]
        print pdb['B'][158]
        a = pdb['A'][79].dist_xyz(pdb['B'][158]['CA'])
        print a