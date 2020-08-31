# types atoms in a zeolite
from typing import List
from ase import Atoms
from ase.neighborlist import natural_cutoffs, NeighborList

def type_atoms(atoms_obj):
    type_dict = {'none': []}

    # Getting neighbor list
    nl = NeighborList(natural_cutoffs(atoms_obj), bothways=True, self_interaction=False)
    nl.update(atoms_obj)

    for i in atoms_obj:
        # labels framework atoms
        if i.symbol in ['Sn', 'Al', 'Si']:
            label = 'framework-%s' % i.symbol
            if label in type_dict.keys():
                type_dict[label].append(i.index)
            else:
                type_dict[label] = [i.index]

        # label carbon and nitrogen as adsorbate
        if i.symbol in ['C', 'N']:
            label = 'adsorbate-%s' % i.symbol
            if label in type_dict.keys():
                type_dict[label].append(i.index)
            else:
                type_dict[label] = [i.index]

        # label Cu, Ni as extraframework
        if i.symbol in ['Cu', 'Ni']:
            label = 'extraframework-%s' % i.symbol
            if label in type_dict.keys():
                type_dict[label].append(i.index)
            else:
                type_dict[label] = [i.index]

        # label oxygens as framework, adsorbate, or bound adsorbate
        if i.symbol == 'O':
            neigh_ind = nl.get_neighbors(i.index)
            if len(neigh_ind) == 2:
                neigh_sym1, neigh_sym2 = atoms_obj[neigh_ind[0][0]].symbol, atoms_obj[neigh_ind[0][1]].symbol
                # assign framework oxygens
                if neigh_sym1 in ['Sn', 'Al', 'Si'] and neigh_sym2 in ['Sn', 'Al', 'Si']:
                    label = 'framework-%s' % i.symbol
                    if label in type_dict.keys():
                        type_dict[label].append(i.index)
                    else:
                        type_dict[label] = [i.index]
                # assign bound adsorbate oxygen
                if neigh_sym1 in ['H', 'C', 'N'] and neigh_sym2 in ['Sn', 'Al', 'Si']:
                    label = 'bound-adsorbate-%s' % i.symbol
                    if label in type_dict.keys():
                        type_dict[label].append(i.index)
                    else:
                        type_dict[label] = [i.index]

                if neigh_sym2 in ['H', 'C', 'N'] and neigh_sym1 in ['Sn', 'Al', 'Si']:
                    label = 'bound-adsorbate-%s' % i.symbol
                    if label in type_dict.keys():
                        type_dict[label].append(i.index)
                    else:
                        type_dict[label] = [i.index]
                # assign rest as adsorbate oxygen
                else:
                    label = 'adsorbate-%s' % i.symbol
                    if label in type_dict.keys():
                        type_dict[label].append(i.index)
                    else:
                        type_dict[label] = [i.index]
            # assign rest as adsorbate oxygen
            else:
                label = 'adsorbate-%s' % i.symbol
                if label in type_dict.keys():
                    type_dict[label].append(i.index)
                else:
                    type_dict[label] = [i.index]

        # all others get 'none' type
        else:
            type_dict['none'].append(i.index)
    return (type_dict)

# testing
if __name__ == '__main__':
    from ase.io import read
    b = read('BEA.cif')
    print(type_atoms(b))