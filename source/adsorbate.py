import numpy as np
from random import random
from scipy.optimize import minimize
from sklearn.cluster import MeanShift
from ase import Atom
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.data import atomic_numbers, covalent_radii
from .zeotype import Zeotype


class Adsorbate(Zeotype):
    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, silent: bool = False, zeolite_type: str = '',
                 t_site_to_atom_indices=None, atom_indices_to_t_site=None):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms, charges, scaled_positions, cell,
                         pbc, celldisp, constraint, calculator, info, velocities, silent, zeolite_type,
                         t_site_to_atom_indices, atom_indices_to_t_site)

        self.host_zeotype = None
        host = self.host_zeotype  # host zeotype object
        self.indices_of_ads_in_host = None
        self.new_positions = None

    @staticmethod  # should this be here?
    def min_dist(pos):
        """
        minimum distance from position to a host atom
        :param pos: vector x,y,z position
        :return: float, minimum distance
        """
        dummy_atom = Atom('H', position=pos)
        dummy_host = host + dummy_atom
        return min(dummy_host.get_distances(-1, [i for i in range(len(host))], mic=True))

    @staticmethod  # should this be here?
    def avg_dist(pos):
        """
        average distance from position to all host atoms
        :param pos: vector x,y,z position
        :return: float, average distance
        """
        dummy_atom = Atom('H', position=pos)
        dummy_host = host + dummy_atom
        return np.average(dummy_host.get_distances(-1, [i for i in range(len(host))], mic=True))

    @staticmethod  # should this be here?
    def find_void(pos):
        """
        finds nearest empty region in host
        :param pos: vector x,y,z position
        :return: vector position of center of empty void
        """
        guess = pos
        func = lambda pos: -1 * min_dist(pos)  # 1 param function for scipy.minimize
        ans = minimize(func, guess)
        return ans.x

    @staticmethod  # should this be here?
    def sphere_sample(radius, num_pts=None):
        """
        generates random positions on the surface of a sphere of certain radius
        :param radius: radius of sphere surface to sample
        :param num_pts: number of points to try
        :return: list of x,y,z positions on surface of sphere
        """
        if num_pts is None:
            num_pts = 500
        vect_list = []
        for i in range(num_pts):
            x, y, z = [2 * random() - 1 for i in range(3)]
            # select points inside of unit sphere
            if x ** 2 + y ** 2 + z ** 2 > 1.0:
                pass
            else:
                unit_vec = [x, y, z] / np.linalg.norm([x, y, z])  # normalize vector length to 1
                vect_list.append(radius * unit_vec)  # lengthen vector to desired radius
        return vect_list

    def generate_poss_place(host, index, radius, cutoff, num_pts=None):
        """
        finds positions near host atom far enough from other framework atoms.
        :param index: index of host atom at center
        :param radius: radius around host atom to test points
        :param cutoff: minimum distance from other host atoms allowed for test points
        :param num_pts: number of points to try
        :return: list. positions of points which meet cutoff criteria
        """
        assert (radius > cutoff)

        guess_pos = sphere_sample(radius, num_pts)
        host_pos = host.get_positions()[index]
        viable_pos = []

        for pos in guess_pos:
            dist = min_dist(pos + host_pos, host)
            if dist > cutoff:
                viable_pos.append(pos + host_pos)
        return viable_pos

    def get_place_clusters(host, index, radius, cutoff, num_pts=None):
        """
        finds positions near host atom far enough from other framework atoms,
        clusters these viable positions and returns the centers of these clusters.
        If number of points is too small will return error
        :param index: index of host atom at center
        :param radius: radius around host atom to test points
        :param cutoff: minimum distance from other host atoms allowed for test points
        :param num_pts: number of points to try
        :return: list. center positions of clusters of points which meet criteria
        """

        viable_pos = generate_poss_place(host, index, radius, cutoff, num_pts)
        ms = MeanShift(bin_seeding=True)
        ms.fit(np.array(viable_pos))
        cluster_centers = ms.cluster_centers_
        return cluster_centers

    def find_best_place(host, index, radius, cutoff, num_pts=None):
        """
        picks the best location to place an adsorbate around the host atom
        :param index: index of host atom at center
        :param radius: radius around host atom to test points
        :param cutoff: minimum distance from other host atoms allowed for test points
        :param num_pts: number of points to try
        :return: array. x,y,z position of best location
        """

        if num_pts is None:
            num_pts = 500
        viable_pos = get_place_clusters(host, index, radius, cutoff, num_pts)
        best_avg_dist = 0
        for pos in viable_pos:
            void = find_void(pos, host)
            void_avg_dist = avg_dist(void)  # average dist at void nearest to pos
            if void_avg_dist > best_avg_dist:
                best_avg_dist = void_avg_dist  # selects pos with largest nearby void
                best_pos = pos
        return best_pos

    def pick_donor(self):
        """
        finds atom in adsorbate most likely to bind to metal
        Heuristic: N > O > P > S > X > C > H
        :param ads: adsorbate atoms object
        :return: index of donor atom in adsorbate
        """

        donors = ['N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'C', 'H']
        ads_symbols = [k.symbol for k in self]
        for symbol in donors:
            if symbol in ads_symbols:
                ads_ind = [s.index for s in self if s.symbol == symbol][
                    0]  # picks first index of atom with specified symbol
                return ads_ind
        else:
            return ValueError("Cannot find donor atom in ads")

    def pick_host_atom(host):
        """
        picks element with largest atomic number higher than 14 (Si), else, Al
        :return: index of host atom
        """

        atom_nums = host.get_atomic_numbers()
        max_num = max(atom_nums)
        symbol = Atom(max_num).symbol  # symbol of atom with highest atomic number
        host_ind = [i.index for i in host if i.symbol == symbol][0]  # picks first index of atom with specified symbol
        if max_num == 14:
            if 13 in atom_nums:
                al_ind = [j.index for j in host if j.symbol == 'Al'][0]
                return al_ind  # return index of an Al atom if no heavy atoms in host
            else:
                return ValueError("No heavy atoms other than Si in host")
        else:
            return host_ind

    def pick_pos(self, host, donor_ind, host_ind, radius=None, cutoff=None):
        """
        Finds a good position to add adsorbate to host
        :param donor_ind: index of donor atom on adsorbate
        :param host_ind: index of host binding site
        :param radius: distance between host atom and donor atom
        :param cutoff: minimum distance donor atom must be from host atoms
        :return: vector of best position to place adsorbate
        """
        donor_atom_symbol, host_atom_symbol = ads[donor_ind].symbol, host[host_ind].symbol
        donor_atom_number, host_atom_number = atomic_numbers[donor_atom_symbol], atomic_numbers[host_atom_symbol]
        donor_radius, host_radius = covalent_radii[donor_atom_number], covalent_radii[host_atom_number]

        if radius is None:
            radius = host_radius + donor_radius
        if cutoff is None:
            cutoff = host_radius

        pos = find_best_place(host, host_ind, radius, cutoff)
        return pos

    def get_donor_vec(self):
        """
        finds direction of lone pair electrons on a donor atom
        :return: vector direction of donor atoms
        """
        nl = NeighborList(natural_cutoffs(self), self_interaction=False, bothways=True)
        nl.update(self)
        # gets neighbors of donor atom and adds the vectors from neighbor to donor
        # for most donor atoms this is roughly in the proper binding direction
        donor_neighbors = nl.get_neighbors(donor_ind)[0]  # neighbor's index
        donor_vec = [0, 0, 0]
        for i in donor_neighbors:
            a = self.get_distance(i, donor_ind, vector=True)
            donor_vec = donor_vec + a
        donor_vec = donor_vec / np.linalg.norm(donor_vec)  # normalizes donor vec
        return donor_vec

    def place_ads(self, host, donor_ind=None, host_ind=None, pos=None):
        """
        Places adsorbate in host according to specified parameters
        :param pos: vector, the position to place adsorbate's donor atom
        :param host_ind: integer, index of site in host which adsorbate will be bound
        :param donor_ind: integer, index of donor atom on adsorbate
        :return: atoms object with adsorbate in host
        """

        if donor_ind is None:
            donor_ind = pick_donor(self)
        if host_ind is None:
            host_ind = pick_host_atom(host)
        if pos is None:
            pos = pick_pos(self, host, donor_ind, host_ind)

        dummy_host = host + Atom('H', position=pos)  # add dummy hydrogen atom to get distances to host atoms
        vec = dummy_host.get_distance(-1, host_ind, mic=True, vector=True)
        donor_vec = get_donor_vec(self)  # get the direction of lone pairs on donor atom

        # rotate the ads into binding direction, move the ads to proper pos and combine
        ads2 = self.copy()  # to avoid making changes to orig adsorbate
        ads2.rotate(donor_vec, vec)
        ads2.translate(pos - ads2.get_positions()[donor_ind])
        host_new = host + ads2
        return host_new

    def integrate_ads(self):
        """
        tricky one
        :return:
        """
        ...


# testing. NOT FUNCTIONAL!
if __name__ == '__main__':
    from ase.io import read
    from ase.visualize import view
    from ase.build import molecule

    host = read('BEA.cif')
    host[185].symbol = 'Sn'  # for visualization
    ads = molecule('CH3OH')
    new_host.place_ads(ads, host)
    view(new_host)
