from maze.zeotypes import ImperfectZeotype
from ase import Atoms
from typing import Union, Tuple

#TODO: At some point make a common ABC for the ExtraFrameworkAtoms and Adsorbate
class ExtraFrameworkAtoms(Atoms):
    """
    This class represents the extraframework atoms that will be added to the
    """
    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None,  host_zeotype=None, name='', description=''):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms, charges, scaled_positions, cell,
                         pbc, celldisp, constraint, calculator, info, velocities)

        assert '_' not in description, 'cannot add _ to description'
        if isinstance(symbols, self.__class__):
            if host_zeotype is None:
                host_zeotype = symbols.host_zeotype
            if description == '':
                description = symbols.description
            if name == '':
                name = symbols.name

        self.host_zeotype = host_zeotype
        self.description = description
        self.name = name


class ExtraFramework(ImperfectZeotype):
    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, site_to_atom_indices=None, atom_indices_to_site=None,
                 additions=None):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms,
                         charges, scaled_positions, cell, pbc, celldisp, constraint,
                         calculator, info, velocities, site_to_atom_indices,
                         atom_indices_to_site, additions)

    def integrate_efa(self, efa: Atoms) -> Tuple['ExtraFramework', ExtraFrameworkAtoms]:
        """
        Integrate extraframework atoms into the extraframework and return new
        extraframework and the extraframework atoms
        :param efa:  Extraframework atoms needed for
        :type efa: Atoms, ExtraFrameworkAtoms
        :return:  Extraframework with efa integrated and ExtraFrameworkAtoms object
        :rtype: Tuple['ExtraFramework', ExtraFrameworkAtoms]
        """
        efa_name = 'extraframework_atoms'
        efa = ExtraFrameworkAtoms(efa)
        new_self = self.add_atoms(efa, efa_name, short_description=efa.description)
        efa.name = new_self.additions[efa_name][-1]
        efa.host_zeotype = new_self
        return new_self, efa
