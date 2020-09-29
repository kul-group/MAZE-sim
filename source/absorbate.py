import numpy as np
from random import random
from scipy.optimize import minimize
from sklearn.cluster import MeanShift
from ase import Atoms, Atom
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

        self.host_zeotype = host  # host atoms object. FIX. turn 'host' -> 'self'
        self.adsorbate = ads      # adsorbate
        self.indices_of_ads_in_host = None
        self.new_positions = None

    def position_atoms(self):
        """
        1. get the positions

        2. assign the positions
        Takes in the other optional arguments that you have

        rotate atoms in adsorbate so that they are correctly positioned to be added to the
        the host zeotype
        :return:
        """
        ...

    @staticmethod  # should this be here?
    def min_dist(pos):
        '''
        minimum distance from position to a host atom
        :param pos: vector x,y,z position
        :return: float, minimum distance
        '''
        dummy_atom = Atom('H', position=pos)
        dummy_host = host + dummy_atom
        return (min(dummy_host.get_distances(-1, [i for i in range(len(host))], mic=True)))

    @staticmethod  # should this be here?
    def avg_dist(pos):
        '''
        average distance from position to all host atoms
        :param pos: vector x,y,z position
        :return: float, average distance
        '''
        dummy_atom = Atom('H', position=pos)
        dummy_host = host + dummy_atom
        return (np.average(dummy_host.get_distances(-1, [i for i in range(len(host))], mic=True)))

    @staticmethod  # should this be here?
    def find_void(pos):
        '''
        finds nearest empty region in host
        :param pos: vector x,y,z position
        :return: vector position of center of empty void
        '''
        guess = pos
        func = lambda pos: -1 * min_dist(pos, host)  # 1 param function for scipy.minimize
        ans = minimize(func, guess)
        return (ans.x)

    @staticmethod  # should this be here?
    def sphere_sample(radius, num_pts=None):
        '''
        generates random positions on the surface of a sphere of certain radius
        :param radius: radius of sphere surface to sample
        :param num_pts: number of points to try
        :return: list of x,y,z positions on surface of sphere
        '''
        if num_pts == None:
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
        return (vect_list)

    def generate_poss_place(host, index, radius, cutoff, num_pts=None):
        '''
        finds positions near host atom far enough from other framework atoms.
        :param index: index of host atom at center
        :param radius: radius around host atom to test points
        :param cutoff: minimum distance from other host atoms allowed for test points
        :param num_pts: number of points to try
        :return: list. positions of points which meet cutoff criteria
        '''
        assert (radius > cutoff)

        guess_pos = sphere_sample(radius, num_pts)
        host_pos = host.get_positions()[index]
        viable_pos = []

        for pos in guess_pos:
            dist = min_dist(pos + host_pos, host)
            if dist > cutoff:
                viable_pos.append(pos + host_pos)
        return(viable_pos)

    def get_place_clusters(host, index, radius, cutoff, num_pts=None):
        '''
        finds positions near host atom far enough from other framework atoms,
        clusters these viable positions and returns the centers of these clusters.
        If number of points is too small will return error
        :param index: index of host atom at center
        :param radius: radius around host atom to test points
        :param cutoff: minimum distance from other host atoms allowed for test points
        :param num_pts: number of points to try
        :return: list. center positions of clusters of points which meet criteria
        '''

        viable_pos = generate_poss_place(host, index, radius, cutoff, num_pts)
        ms = MeanShift(bin_seeding=True)
        ms.fit(np.array(viable_pos))
        cluster_centers = ms.cluster_centers_
        return (cluster_centers)

    def integrate_ads(self):
        """
        tricky one
        :return:
        """
        ...
