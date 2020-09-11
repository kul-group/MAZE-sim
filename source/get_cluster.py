from ase.neighborlist import NeighborList, natural_cutoffs

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
    atoms.set_tags(0)
    atoms[index].tag = 1
    clust_ind = [index]

    while len(clust_ind) < size:
        next_expand = [i.index for i in atoms if i.tag == 1]  # list of next atoms to expand
        for j in next_expand:
            atoms[j].tag = 2  # change tag so atom not expanded in future
            for k in nl.get_neighbors(j)[0]:
                if atoms[k].tag == 0:
                    atoms[k].tag = 1  # tag atoms which are adjacent to cluster atoms, but not part of cluster yet
                    clust_ind.append(k)

    return(clust_ind)

# testing
if __name__ == '__main__':
    from ase.io import read
    b = read('BEA.cif')
    a = get_cluster(b, 101, 5)
    print(a)
