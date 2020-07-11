from typing import List, Dict
from ase import Atoms
from ase.neighborlist import natural_cutoffs, NeighborList

class Zeotype(Atoms):
    """
    This is a Zeotype class that inherits from Atoms. It represents a Zeolite.
    """
    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, silent=False, zeolite_type: str = ''):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms, charges, scaled_positions,
                         cell, pbc, celldisp, constraint, calculator, info, velocities)

        self.zeolite_type = zeolite_type
        self.sites: List[str] = []

    @staticmethod
    def build_from_atoms(a: Atoms, silent=False) -> "Zeotype":
        """
        :param a: Atoms object that Zeotype object will be used to create the Zeotype object
        :param silent: currently does nothing, but will set if output is displayed
        :return: Zeotype object created from Atoms object
        """
        return Zeotype(a.symbols, a.positions, a.numbers, a.tags, a.momenta, a.masses, a.magmoms,
                       a.charges, a.scaled_positions, a.cell, a.pbc, a.celldisp, a.constraint,
                       a.calculator, a.info, a.velocities, silent)

        super() #currently this code does nothing
        self.atoms = atoms
        self.types = {'num': 0, 'unique': [], 'indices': {}, 'count': {}}
        indices, count, indices_atomtype, count_atomtype, atomtype_list = self.analyze_zeolite_atoms()
        self.types['indices'] = indices_atomtype
        self.types['unique'] = np.unique(atomtype_list)
        self.types['num'] = len(np.unique(atomtype_list))
        self.types['count'] = count_atomtype
        self.types['list'] = atomtype_list
        self.silent = silent

    def get_sites(self) -> List[str]:
        """
        :return: returns Zeotype sites
        """
        return self.sites

    def get_zeolite_type(self) -> List[str]:
        """
        :return: returns zeotype sites as a list
        """
        return self.zeolite_type

    def get_atom_types(self) -> Dict[str, List[int]]:
        """
        :return: Returns a dictionary of atom types where the key consists of the atom category
        (framework, absorbate, extraframework or other) followed by -atom chemical symbol. For
        example a Sn atom is a framework atom so the key is "framework-Sn". The value of the
        returned dictionary is a list of the indices of all of the atoms that belong to the
        category.
        """

        type_dict: Dict[str, List[int]] = {}

        # Getting neighbor list
        nl = NeighborList(natural_cutoffs(self), bothways=True, self_interaction=False)
        nl.update(self)

        for atom in self:  # iterate through atom objects in zeotype
            # labels framework atoms
            print(atom.symbol)
            if atom.symbol in ['Sn', 'Al', 'Si']:
                label = 'framework-%s' % atom.symbol
                if label in type_dict.keys():
                    type_dict[label].append(atom.index)
                else:
                    type_dict[label] = [atom.index]

            # label carbon and nitrogen as adsorbate
            if atom.symbol in ['C', 'N']:
                label = 'adsorbate-%s' %atom.symbol
                if label in type_dict.keys():
                    type_dict[label].append(atom.index)
                else:
                    type_dict[label] = [atom.index]

            # label Cu, Ni as extraframework
            if atom.symbol in ['Cu', 'Ni']:
                label = 'extraframework-%s' % atom.symbol
                if label in type_dict.keys():
                    type_dict[label].append(atom.index)
                else:
                    type_dict[label] = [atom.index]

            # all others get 'other' type
            else:
                type_dict['other'].append(atom.index)

        return type_dict

# testing
if __name__ == '__main__':
    from ase.io import read
    b = read('BEA.cif')
    z = Zeotype(b)
    z.get_atom_types()



