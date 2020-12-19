# types atoms in a zeolite
from typing import List
from ase import Atoms
from ase.neighborlist import natural_cutoffs, NeighborList
from collections import defaultdict
def type_atoms(atoms_obj):
    """
    Args:
        atoms_obj:
    """
    type_dict = defaultdict(list)

    # Getting neighbor list
    nl = NeighborList(natural_cutoffs(atoms_obj), bothways=True, self_interaction=False)
    nl.update(atoms_obj)

    for i in atoms_obj:
        # labels framework atoms
        if i.symbol in ['Sn', 'Al', 'Si']:
            label = 'framework-%s' % i.symbol
        # label carbon and nitrogen as adsorbate
        elif i.symbol in ['C', 'N']:
            label = 'adsorbate-%s' % i.symbol

        # label Cu, Ni as extraframework
        elif i.symbol in ['Cu', 'Ni']:
            label = 'extraframework-%s' % i.symbol

        # label oxygens as framework, adsorbate, or bound adsorbate
        elif i.symbol == 'O':
            neigh_ind = nl.get_neighbors(i.index)
            if len(neigh_ind) == 2:
                neigh_sym1, neigh_sym2 = atoms_obj[neigh_ind[0][0]].symbol, atoms_obj[neigh_ind[0][1]].symbol
                # assign framework oxygens
                if neigh_sym1 in ['Sn', 'Al', 'Si'] and neigh_sym2 in ['Sn', 'Al', 'Si']:
                    label = 'framework-%s' % i.symbol
                # assign bound adsorbate oxygen
                elif neigh_sym1 in ['H', 'C', 'N'] and neigh_sym2 in ['Sn', 'Al', 'Si']:
                    label = 'bound-adsorbate-%s' % i.symbol
                elif neigh_sym2 in ['H', 'C', 'N'] and neigh_sym1 in ['Sn', 'Al', 'Si']:
                    label = 'bound-adsorbate-%s' % i.symbol
                # assign rest as adsorbate oxygen
                else:
                    label = 'adsorbate-%s' % i.symbol
            # assign rest as adsorbate oxygen
            else:
                label = 'adsorbate-%s' % i.symbol
        # all others get 'none' type
        else:
            label = "None"
        type_dict[label].append(i.index)

    return dict(type_dict)

# testing
if __name__ == '__main__':
    from ase.io import read
    b = read('BEA.cif')
    print(type_atoms(b))