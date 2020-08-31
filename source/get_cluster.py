from ase.neighborlist import NeighborList, natural_cutoffs
from ase import Atoms
import numpy as np
from ase.io import read, write

# Obtains the indicies of a cluster in a zeolite given center atom and cluster size

def get_cluster(atoms, index, size)
    '''
    :param atoms: atoms object 
    :param size: number of T-sites in cluster, integer
    :param index: index at center of cluster, integer
    :return: list of indicies included in the cluster
    '''

    if type(index) != int:
        raise ValueError('index must be an integer')
    if type(size) != int:
        raise ValueError('size must be an integer')
    try:
        atoms[index]
    except:
        raise ValueError('atoms must be atoms object containing atoms[index]')

    nl = NeighborList(natural_cutoffs(atoms), self_interaction=False, bothways=True)
    nl.update(atoms)
    atoms_cluster = Atoms()
    atoms_cluster.append(atoms[index])
    atoms.set_tags(0)
    atoms[index].tag = 1




# testing
if __name__ == '__main__':
    from ase.io import read
    b = read('BEA.cif')
    a = get_cluster(b, 100, 5)