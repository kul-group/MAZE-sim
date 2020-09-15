# functions for capping a zeolite cluster with Si-OH
# assumes already has neighbor list
import numpy as np
from ase import Atom

def needs_neighb(clust, index):
    # checks if an atom in a cluster has missing bonds
    # only checks O, Si, and common metal atoms
    a = clust[index]  # atom object
    if a.symbol == 'O'and len(nl.get_neighbors(index)[0]) < 2:
        return(True)
    elif a.symbol in ['Si', 'Sn', 'Al', 'Ga', 'B'] and len(nl.get_neighbors(index)[0]) < 4:
        return(True)
    else:
        return(False)

def get_cap_ox(clust):
    ox_inds = [i.index for i in clust if i.symbol == 'O']
    ox_to_cap = [j for j in ox_inds if needs_neighb(clust, j)]
    return(ox_to_cap)

def get_cap_si(clust):
    si_inds = [i.index for i in clust if i.symbol == 'Si']
    si_to_cap = [j for j in si_inds if needs_neighb(clust, j)]
    return(si_to_cap)

def add_cap_ox(clust):
    new_clust = clust
    cap_inds = get_cap_si(clust)
    for ind in cap_inds:
        while len(nl.get_neighbors(ind)[0]) < 4:
            neighb = nl.get_neighbors(ind)[0][0]  # first index in the list of neighbor indicies
            direction = clust.get_positions()[ind] - clust.get_positions()[neighb]  # vector from neighbor to Si
            ox_pos = clust.get_positions()[ind] + 1.6*(clust.get_positions()[ind] + direction) / np.linalg.norm(direction)
            new_ox = Atom('O', position = ox_pos)
            new_clust.append(new_ox)
            nl.update(clust)
    return (new_clust)

def add_cap_h(clust):
    nl.update(clust)
    new_clust = clust
    cap_inds = get_cap_ox(clust)
    for ind in cap_inds:
        neighb = nl.get_neighbors(ind)[0][0]  # first index in the list of neighbor indicies
        direction = clust.get_positions()[ind] - clust.get_positions()[neighb]  # vector from neighbor to oxygen
        h_pos = clust.get_positions()[ind] + (clust.get_positions()[ind] + direction)/np.linalg.norm(direction)
        new_h = Atom('H', position=h_pos)
        new_clust.append(new_h)
    return(new_clust)

def cap_clust(clust):
    si_capped = add_cap_ox(clust)     # cluster with capped si
    ox_capped = add_cap_h(si_capped)  # cluster with capped o
    return(ox_capped)

# testing
if __name__ == '__main__':
    from ase.io import read
    from ase.visualize import view
    from zeotype import Cluster
    from ase.neighborlist import NeighborList, natural_cutoffs
    b = read('BEA.cif')
    c = Cluster(b, index=100, cluster_size=5)
    nl = NeighborList(natural_cutoffs(c), bothways=True, self_interaction=False)
    nl.update(c)
    capped_cluster = cap_clust(c)
    view(capped_cluster)


