from ase import Atoms
from .zeotype import Zeotype


class Absorbate(Zeotype):
    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, silent: bool = False, zeolite_type: str = '',
                 t_site_to_atom_indices=None, atom_indices_to_t_site=None):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms, charges, scaled_positions, cell,
                         pbc, celldisp, constraint, calculator, info, velocities, silent, zeolite_type,
                         t_site_to_atom_indices, atom_indices_to_t_site)

        self.host_zeotype = None
        self.indices_of_abs_in_host = None
        self.new_positions = None

    def position_atoms(self):
        """
        1. get the positions

        2. assign the positions
        Takes in the other optional arguments that you have

        rotate atoms in absorbate so that they are correctly positioned to be added to the
        the host zeotype
        :return:
        """
        ...

    def integrate_abs(self):
        """
        tricky one
        :return:
        """
        ...

    def test(self):
        """

        :return:
        """
        ...
