from ase.neighborlist import NeighborList, natural_cutoffs

# Obtains the indicies of a cluster in a zeolite given center atom and cluster size

def get_cluster(atoms, index, size):
    """
    Args:
        atoms: atoms object
        index: index at center of cluster, integer
        size: min number of atoms in cluster (not number of T-sites!), integer

    Returns:
        list of indicies included in the cluster
    """

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

    # uses tagging to find indicies of cluster
    # intentionally overshoots the input size to make a symmetric cluster

    while len(clust_ind) < size:
        adj_lists = [nl.get_neighbors(j)[0] for j in clust_ind] # list of lists, lists are neighbors of ea. cluster atom
        for ind_list in adj_lists:
            for ind in ind_list:
                if atoms[ind].tag == 0:
                    clust_ind.append(ind)
                    atoms[ind].tag = 1  # prevents repeat addition of an index

    return(clust_ind)

# testing
if __name__ == '__main__':
    from ase.io import read
    b = read('BEA.cif')
    a = get_cluster(b, 101, 5)
    print(a)
