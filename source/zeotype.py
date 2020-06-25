from ase import Atoms


class Zeotype(Atoms):
    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, silent=False):
        super().__init__(self, symbols, positions, numbers, tags, momenta, masses, magmoms, charges, scaled_positions,
                         cell, pbc, celldisp, constraint, calculator, info, velocities, silent)

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

    # sams branch test
    # test
