from ase import Atoms
from .zeotype import Zeotype


class Adsorbate(Zeotype):
    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, silent: bool = False, zeolite_type: str = '',
                 t_site_to_atom_indices=None, atom_indices_to_t_site=None):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms, charges, scaled_positions, cell,
                         pbc, celldisp, constraint, calculator, info, velocities, silent, zeolite_type,
                         t_site_to_atom_indices, atom_indices_to_t_site)

        self.host_zeotype = host
        self.adsorbate = ads
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
    #@staticmethod  # should this be here?
    def min_dist(pos):
        '''
        minimum distance from position to a host atom
        :param pos: vector x,y,z position
        :return: float, minimum distance
        '''
        dummy_atom = Atom('H', position=pos)
        dummy_host = host + dummy_atom
        return (min(dummy_host.get_distances(-1, [i for i in range(len(host))], mic=True)))

    # @staticmethod  # should this be here?
    def avg_dist(pos):
        '''
        average distance from position to all host atoms
        :param pos: vector x,y,z position
        :return: float, average distance
        '''
        dummy_atom = Atom('H', position=pos)
        dummy_host = host + dummy_atom
        return (np.average(dummy_host.get_distances(-1, [i for i in range(len(host))], mic=True)))

    def integrate_ads(self):
        """
        tricky one
        :return:
        """
        ...
