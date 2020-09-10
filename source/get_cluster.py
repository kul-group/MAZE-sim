from ase.neighborlist import NeighborList, natural_cutoffs
from ase import Atoms
import numpy as np
from ase.io import read, write

# Obtains the indicies of a cluster in a zeolite given center atom and cluster size

def get_cluster(atoms, index, size):
    '''
    :param atoms: atoms object 
    :param size: number of atoms in cluster (not number of T-sites!), integer
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
    clust_ind = []

    # uses tagging system to find the cluster
    while len(clust_ind) < size:
        for i in atoms:
            if i.tag == 1:
                adj = nl.get_neighbors(i.index)[0]
                i.tag = 2  # change tag so atom not expanded in future
                for j in adj:
                    if atoms[j].tag == 0:
                        if len(clust_ind) < size:
                            clust_ind.append(j)
                            atoms[j].tag = 1 # neighbors of newly added atoms will be expanded next

    return(clust_ind)



# testing
if __name__ == '__main__':
    from ase.io import read
    b = read('BEA.cif')
    a = get_cluster(b, 101, 5)
    print(a)
