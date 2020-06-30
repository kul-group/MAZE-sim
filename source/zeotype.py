from typing import List
from ase import Atoms
from ase.neighborlist import natural_cutoffs, NeighborList

class Zeotype(Atoms):
    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, silent=False, zeolite_type: str = ''):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms, charges, scaled_positions,
                         cell, pbc, celldisp, constraint, calculator, info, velocities)

        self.zeolite_type = zeolite_type
        self.sites: List[str] = []

    @staticmethod
    def build_from_atoms(a: Atoms, silent=False) -> "Zeotype":
        return Zeotype(a.symbols, a.positions, a.numbers, a.tags, a.momenta, a.masses, a.magmoms,
                       a.charges, a.scaled_positions, a.cell, a.pbc, a.celldisp, a.constraint,
                       a.calculator, a.info, a.velocities, silent)

        super()
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
        return self.sites

    def get_zeolite_type(self) -> List[str]:
        return self.zeolite_type

    def type_atoms(self):
        type_dict = {}

        # Getting neighbor list
        nl = NeighborList(natural_cutoffs(self), bothways=True, self_interaction=False)
        nl.update(self)

        for i in self:
            # labels framework atoms
            if i.symbol in ['Sn', 'Al', 'Si']:
                label = 'framework-%s' %i.symbol
                if label in type_dict.keys():
                    type_dict[label].append(i.index)
                else:
                    type_dict[label] = [i.index]

        return(type_dict)


if __name__ == '__main__':
    from ase.io import read
    b = read('BEA.cif')
    z = Zeotype(b)
    print(z.type_atoms())


