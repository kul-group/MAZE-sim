import warnings
import numpy as np
from ase import Atom, Atoms
from ase.data import atomic_numbers, covalent_radii
from ase.neighborlist import NeighborList, natural_cutoffs
from scipy.optimize import minimize
from sklearn.cluster import MeanShift


class Adsorbate(Atoms):
    """
    This is an adsorbate class which represents an adsorbate
    """
    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None,  host_zeotype=None, name='', description=''):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms, charges, scaled_positions, cell,
                         pbc, celldisp, constraint, calculator, info, velocities)

        assert '_' not in description, 'cannot add _ to description'
        if isinstance(symbols, Adsorbate):
            if host_zeotype is None:
                host_zeotype = symbols.host_zeotype
            if description == '':
                description = symbols.description
            if name == '':
                name = symbols.name

        self.host_zeotype = host_zeotype
        self.description = description
        self.name = name

    def min_distance(self, ads_position) -> float:
        """minimum distance from atom in a host to adsorbate :param
        ads_position: np.array x,y,z of adsorbate position :return: float,
        minimum distance between an atom in host MAZE-sim and adsorbate position

        Args:
            ads_position:
        """
        assert self.host_zeotype is not None, "Cannot find min distance when host MAZE-sim is None"

        dummy_atom = Atom('H', position=ads_position)
        dummy_host = self.host_zeotype + dummy_atom
        min_distance = min(dummy_host.get_distances(-1, [i for i in range(len(self.host_zeotype))], mic=True))
        return min_distance

    def avg_distance(self, ads_position):
        """average distance from position to all host atoms :param ads_position:
        np.array x,y,z of adsorbate position :return: float, average distance of
        host MAZE-sim atoms to adsorbate position

        Args:
            ads_position:
        """
        assert self.host_zeotype is not None, "Cannot find average distance when host MAZE-sim is None"

        dummy_atom = Atom('H', position=ads_position)
        dummy_host = self.host_zeotype + dummy_atom
        avg_distance = np.average(dummy_host.get_distances(-1, [i for i in range(len(self.host_zeotype))], mic=True))
        return avg_distance

    def find_void(self, void_position_guess):
        """finds nearest empty region in host MAZE-sim :param
        void_position_guess: An initial guess for the center of the void as a
        np.array with x,y,z position :return: np.array position of center of
        empty void

        Args:
            void_position_guess:
        """

        # TODO: Find a way to speed this up

        assert self.host_zeotype is not None, "Cannot find void position when host MAZE-sim is None"

        fun_to_min = lambda pos: -1 * self.min_distance(pos)  # 1 param function for scipy.minimize
        ans = minimize(fun_to_min, void_position_guess)
        return ans.x

    @staticmethod
    def sphere_sample(radius, num_pts=500):
        """generates random positions on the surface of a sphere of certain
        radius :param radius: radius of sphere surface to sample :param num_pts:
        number of points to try :return: list of x,y,z positions on surface of
        sphere

        Args:
            radius:
            num_pts:
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
        """finds positions near host atom far enough from other framework atoms.
        :param index: index of host atom at center :param radius: radius around
        host atom to tests points :param cutoff: minimum distance from other host
        atoms allowed for tests points :param num_pts: number of points to try
        :return: list. positions of points which meet cutoff criteria

        Args:
            index:
            radius:
            cutoff:
            num_pts:
        """
        assert (radius > cutoff), "radius larger than cutoff distance"
        assert self.host_zeotype is not None, "host MAZE-sim cannot be none"

        guess_positions = self.sphere_sample(radius, num_pts)
        host_pos = self.host_zeotype.get_positions()[index]
        viable_positions = []

        for pos in guess_positions:
            dist = self.min_distance(pos + host_pos)
            if dist > cutoff:
                viable_positions.append(pos + host_pos)

        return viable_positions

    def get_viable_pos_cluster_centers(self, index, radius, cutoff, num_pts=None):
        """finds positions near host atom far enough from other framework atoms,
        clusters these viable positions and returns the centers of these
        clusters. If number of points is too small will return error :param
        index: index of host atom at center :param radius: radius around host
        atom to tests points :param cutoff: minimum distance from other host
        atoms allowed for tests points :param num_pts: number of points to try
        :return: list. center positions of clusters of points which meet
        criteria

        Args:
            index:
            radius:
            cutoff:
            num_pts:
        """

        viable_pos = self.get_viable_positions(index, radius, cutoff, num_pts)
        ms = MeanShift(bin_seeding=True)
        ms.fit(np.array(viable_pos))
        cluster_centers = ms.cluster_centers_
        return cluster_centers

    def find_best_place(self, index, radius, cutoff, num_pts=500):
        """picks the best location to place an adsorbate around the host atom
        :param index: index of host atom at center :param radius: radius around
        host atom to tests points :param cutoff: minimum distance from other host
        atoms allowed for tests points :param num_pts: number of points to try
        :return: array. x,y,z position of best location

        Args:
            index:
            radius:
            cutoff:
            num_pts:
        """

        best_positions = self.get_viable_pos_cluster_centers(index, radius, cutoff, num_pts)
        best_avg_dist = 0
        for pos in best_positions:
            void = self.find_void(pos)
            void_avg_dist = self.avg_distance(void)  # average dist at void nearest to pos
            best_pos = None
            if void_avg_dist > best_avg_dist:
                best_avg_dist = void_avg_dist  # selects pos with largest nearby void
                best_pos = pos

        if best_pos is None:
            assert False, 'No good positions at specified index'
        return best_pos

    def pick_donor(self):
        """finds atom in adsorbate most likely to bind to metal Heuristic: N > O
        > P > S > X > C > H :param ads: adsorbate atoms object :return: index of
        donor atom in adsorbate
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
        """picks element with largest atomic number higher than 14 (Si) if one exists. If none exist,
            pick Al as a donor atom.

        Returns:
            index of host atom
        """
        sym_index_map, count = self.host_zeotype.count_elements()
        max_atomic_num_sym = max(sym_index_map.keys(), key=lambda x: atomic_numbers[x])
        if atomic_numbers[max_atomic_num_sym] > atomic_numbers['Si']:
            return sym_index_map[max_atomic_num_sym][0]
        else:
            try:
                return sym_index_map['Al'][0]
            except KeyError:
                print("No Al in host MAZE-sim")
                raise

    def pick_ads_position(self, donor_ind, host_ind, radius=None, cutoff=None):
        """Finds a good position to add adsorbate to host :param donor_ind:
        index of donor atom on adsorbate :param host_ind: index of host binding
        site :param radius: distance between host atom and donor atom :param
        cutoff: minimum distance donor atom must be from host atoms :return:
        vector of best position to place adsorbate

        Args:
            donor_ind:
            host_ind:
            radius:
            cutoff:
        """
        donor_atom_symbol, host_atom_symbol = self[donor_ind].symbol, self.host_zeotype[host_ind].symbol
        donor_atom_number, host_atom_number = atomic_numbers[donor_atom_symbol], atomic_numbers[host_atom_symbol]
        donor_radius, host_radius = covalent_radii[donor_atom_number], covalent_radii[host_atom_number]

        if radius is None:
            radius = host_radius + donor_radius
        if cutoff is None:
            cutoff = host_radius

        pos = self.find_best_place(host_ind, radius, cutoff)
        return pos

    def get_donor_vec(self, donor_index):
        """finds direction of lone pair electrons on an adsorbate donor atom
        :return: vector direction of donor atoms

        Args:
            donor_index:
        """
        nl = NeighborList(natural_cutoffs(self), self_interaction=False, bothways=True)
        nl.update(self)
        # gets neighbors of donor atom and adds the vectors from neighbor to donor
        # for most donor atoms this is roughly in the proper binding direction
        donor_neighbors = nl.get_neighbors(donor_index)[0]  # neighbor's index
        donor_vec = np.array([0, 0, 0])
        for i in donor_neighbors:
            a = self.get_distance(i, donor_index, vector=True)
            donor_vec = donor_vec + a
        if np.linalg.norm(donor_vec) == 0:
            warnings.warn("donor vector with magnitude 0 found, providing default  vector")
            return np.array([1, 0, 0])

        donor_vec = donor_vec / np.linalg.norm(donor_vec)  # normalizes donor vec
        return donor_vec

    def position_ads(self, donor_ind=None, host_ind=None, pos=None) -> "Adsorbate":
        """Rotates and positions adsorbate according to specified parameters if
        no parameters are provided, a reasonable position is found :param pos:
        vector, the position to place adsorbate's donor atom :param host_ind:
        integer, index of site in host which adsorbate will be bound :param
        donor_ind: integer, index of donor atom on adsorbate :return: atoms
        object with adsorbate in host

        Args:
            donor_ind:
            host_ind:
            pos:
        """

        if donor_ind is None:
            donor_ind = self.pick_donor()
        if host_ind is None:
            host_ind = self.pick_host_atom()
        if pos is None:
            pos = self.pick_ads_position(donor_ind, host_ind)

        dummy_host = self.host_zeotype + Atom('H', position=pos)  # add dummy hydrogen atom to get distances to host atoms
        vec = dummy_host.get_distance(-1, host_ind, mic=True, vector=True)
        donor_vec = self.get_donor_vec(donor_ind)  # get the direction of lone pairs on donor atom

        new_self = self.__class__(self)
        # rotate the ads into binding direction, move the ads to proper pos and combine
        new_self.rotate(donor_vec, vec)
        new_self.translate(pos - self.get_positions()[donor_ind])
        return new_self


