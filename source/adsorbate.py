import numpy as np
from random import random
from scipy.optimize import minimize
from sklearn.cluster import MeanShift
from ase import Atom, Atoms
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.data import atomic_numbers, covalent_radii
import source.zeotype


class Adsorbate(Atoms):
    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None,  host_zeotype=None, ads_to_zeotype_index_map=None,
                 new_positions=None):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms, charges, scaled_positions, cell,
                         pbc, celldisp, constraint, calculator, info, velocities)

        self.host_zeotype = host_zeotype
        self.ads_to_zeotype_index_map = ads_to_zeotype_index_map
        self.new_positions = new_positions

    def min_distance(self, ads_position) -> float:
        """
        minimum distance from atom in a host to adsorbate
        :param ads_position: np.array x,y,z of adsorbate position
        :return: float, minimum distance between an atom in host zeotype and adsorbate position
        """
        assert self.host_zeotype is not None, "Cannot find min distance when host zeotype is None"

        dummy_atom = Atom('H', position=ads_position)
        dummy_host = self.host_zeotype + dummy_atom
        min_distance = min(dummy_host.get_distances(-1, [i for i in range(len(self.host_zeotype))], mic=True))
        return min_distance

    def avg_distance(self, ads_position):
        """
        average distance from position to all host atoms
        :param ads_position: np.array x,y,z of adsorbate position
        :return: float, average distance of host zeotype atoms to adsorbate position
        """
        assert self.host_zeotype is not None, "Cannot find average distance when host zeotype is None"

        dummy_atom = Atom('H', position=ads_position)
        dummy_host = self.host_zeotype + dummy_atom
        avg_distance = np.average(dummy_host.get_distances(-1, [i for i in range(len(self.host_zeotype))], mic=True))
        return avg_distance

    def find_void(self, void_position_guess):
        """
        finds nearest empty region in host zeotype
        :param void_position_guess: An initial guess for the center of the void as a np.array with x,y,z position
        :return: np.array position of center of empty void
        """

        # TODO: Find a way to speed this up

        assert self.host_zeotype is not None, "Cannot find void position when host zeotype is None"

        fun_to_min = lambda pos: -1 * self.min_distance(pos)  # 1 param function for scipy.minimize
        ans = minimize(fun_to_min, void_position_guess)
        return ans.x

    @staticmethod
    def sphere_sample(radius, num_pts=500):
        """
        generates random positions on the surface of a sphere of certain radius
        :param radius: radius of sphere surface to sample
        :param num_pts: number of points to try
        :return: list of x,y,z positions on surface of sphere
        """
        position_list = []
        for _ in range(int(num_pts)):
            # see https://stackoverflow.com/questions/33976911/
            # generate-a-random-sample-of-points-distributed-on-the-surface-of-a-unit-sphere/33977530#33977530
            # for discussion on this algorithm

            vec = np.random.normal(0, 1, 3)  # select three random points (if normal dist no skip needed)
            vec /= np.linalg.norm(vec)  # normalize vector
            vec *= radius  # lengthen vector to desired radius
            position_list.append(list(vec))

        return position_list

    def get_viable_positions(self, index, radius, cutoff, num_pts=None):
        """
        finds positions near host atom far enough from other framework atoms.
        :param index: index of host atom at center
        :param radius: radius around host atom to test points
        :param cutoff: minimum distance from other host atoms allowed for test points
        :param num_pts: number of points to try
        :return: list. positions of points which meet cutoff criteria
        """
        assert (radius > cutoff), "radius larger than cutoff distance"
        assert self.host_zeotype is not None, "host zeotype cannot be none"

        guess_positions = self.sphere_sample(radius, num_pts)
        host_pos = self.host_zeotype.get_positions()[index]
        viable_positions = []

        for pos in guess_positions:
            dist = self.min_distance(pos + host_pos)
            if dist > cutoff:
                viable_positions.append(pos + host_pos)

        return viable_positions

    def get_viable_pos_cluster_centers(self, index, radius, cutoff, num_pts=None):
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

        viable_pos = self.get_viable_positions(index, radius, cutoff, num_pts)
        ms = MeanShift(bin_seeding=True)
        ms.fit(np.array(viable_pos))
        cluster_centers = ms.cluster_centers_
        return cluster_centers

    def find_best_place(self, index, radius, cutoff, num_pts=500):
        """
        picks the best location to place an adsorbate around the host atom
        :param index: index of host atom at center
        :param radius: radius around host atom to test points
        :param cutoff: minimum distance from other host atoms allowed for test points
        :param num_pts: number of points to try
        :return: array. x,y,z position of best location
        """

        best_positions = self.get_viable_pos_cluster_centers(index, radius, cutoff, num_pts)
        best_avg_dist = 0
        for pos in best_positions:
            void = self.find_void(pos)
            void_avg_dist = self.avg_distance(void)  # average dist at void nearest to pos
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

        ads_symbols = {}
        for atom in self:
            if atom.symbol not in ads_symbols:
                ads_symbols[atom.symbol] = atom.index  # add first occurrence of atom symbol to dict

        donor_symbols = ['N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'C', 'H']
        for donor_symbol in donor_symbols:
            if donor_symbol in ads_symbols:
                ads_index = ads_symbols[donor_symbol]  # picks first index of atom with specified symbol
                return ads_index
        else:
            raise ValueError("Cannot find viable donor atom in adsorbate")

    def pick_host_atom(self) -> int:
        """
        picks element with largest atomic number higher than 14 (Si) if one exists. If none exist,
         pick Al as a donor atom.
        :return: index of host atom
        """
        sym_index_map, count = self.host_zeotype.count_elements()
        max_atomic_num_sym = max(sym_index_map.keys(), key=lambda x: atomic_numbers[x])
        if atomic_numbers[max_atomic_num_sym] > atomic_numbers['Si']:
            return sym_index_map[max_atomic_num_sym][0]
        else:
            try:
                return sym_index_map['Al'][0]
            except KeyError:
                print("No Al in host zeotype")
                raise

    def pick_ads_position(self, donor_ind, host_ind, radius=None, cutoff=None):
        """
        Finds a good position to add adsorbate to host
        :param donor_ind: index of donor atom on adsorbate
        :param host_ind: index of host binding site
        :param radius: distance between host atom and donor atom
        :param cutoff: minimum distance donor atom must be from host atoms
        :return: vector of best position to place adsorbate
        """
        donor_atom_symbol, host_atom_symbol = ads[donor_ind].symbol, self.host_zeotype[host_ind].symbol
        donor_atom_number, host_atom_number = atomic_numbers[donor_atom_symbol], atomic_numbers[host_atom_symbol]
        donor_radius, host_radius = covalent_radii[donor_atom_number], covalent_radii[host_atom_number]

        if radius is None:
            radius = host_radius + donor_radius
        if cutoff is None:
            cutoff = host_radius

        pos = self.find_best_place(host_ind, radius, cutoff)
        return pos

    def get_donor_vec(self, donor_index):
        """
        finds direction of lone pair electrons on an adsorbate donor atom
        :return: vector direction of donor atoms
        """
        nl = NeighborList(natural_cutoffs(self), self_interaction=False, bothways=True)
        nl.update(self)
        # gets neighbors of donor atom and adds the vectors from neighbor to donor
        # for most donor atoms this is roughly in the proper binding direction
        donor_neighbors = nl.get_neighbors(donor_index)[0]  # neighbor's index
        donor_vec = [0, 0, 0]
        for i in donor_neighbors:
            a = self.get_distance(i, donor_index, vector=True)
            donor_vec = donor_vec + a

        donor_vec = donor_vec / np.linalg.norm(donor_vec)  # normalizes donor vec
        return donor_vec

    def position_ads(self, donor_ind=None, host_ind=None, pos=None):
        """
        Rotates and positions adsorbate according to specified parameters
        if no parameters are provided, a reasonable position is found
        :param pos: vector, the position to place adsorbate's donor atom
        :param host_ind: integer, index of site in host which adsorbate will be bound
        :param donor_ind: integer, index of donor atom on adsorbate
        :return: atoms object with adsorbate in host
        """

        if donor_ind is None:
            donor_ind = self.pick_donor()
        if host_ind is None:
            host_ind = self.pick_host_atom()
        if pos is None:
            pos = self.pick_ads_position(donor_ind, host_ind)

        dummy_host = self.host_zeotype + Atom('H',
                                                             position=pos)  # add dummy hydrogen atom to get distances to host atoms
        vec = dummy_host.get_distance(-1, host_ind, mic=True, vector=True)
        donor_vec = self.get_donor_vec(donor_ind)  # get the direction of lone pairs on donor atom

        # rotate the ads into binding direction, move the ads to proper pos and combine
        self.rotate(donor_vec, vec)
        self.translate(pos - self.get_positions()[donor_ind])

    def integrate_ads(self, position_ads=False):
        """
        append adsorbate to main atom
        :return:
        """
        assert self not in self.host_zeotype.adsorbates, "cannot integrate single adsorbate object twice"
        if position_ads:
            self.position_ads()

        self.host_zeotype.extend(self)
        self.set_ads_zeotype_map()
        self.host_zeotype.adsorbates.append(self)


    def set_ads_zeotype_map(self):
        """
        Get the mapping between an adsorbate and the zeolite where it was added
        :return: a zeotype to adsorbate index dictionary
        """
        ads_position_index_map = {}
        for atom in self:
            ads_position_index_map[str(atom.position)] = atom.index

        zeotype_position_index_map = {}
        for atom in self.host_zeotype:
            zeotype_position_index_map[str(atom.position)] = atom.index

        ads_to_zeotype_index_map = {}

        for pos, index in ads_position_index_map.items():
            ads_to_zeotype_index_map[index] = zeotype_position_index_map[pos]

        self.ads_to_zeotype_index_map = ads_to_zeotype_index_map

    def remove_ads(self):
        """
        also tricky
        :return:
        """
        ...

if __name__ == '__main__':
    from ase.io import read
    from ase.visualize import view
    from ase.build import molecule

    zeotype = source.zeotype.Zeotype.build_from_cif_with_labels('BEA.cif')
    print(zeotype[30].symbol)
    zeotype[20].symbol = 'Sn'  # for visualization
    mol = molecule('CH3OH')
    ads = Adsorbate(mol)
    ads.host_zeotype = zeotype
    ads.position_ads()
    view(ads)
    view(zeotype)
    ads.integrate_ads()
    view(zeotype)
