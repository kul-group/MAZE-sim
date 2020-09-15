# functions for capping a zeolite cluster with Si-OH
# assumes already
def needs_neighb(clust, index):
    # checks if an atom in a cluster has missing bonds
    # only checks O, Si, and common metal atoms
    a = clust[index] # atom object
    if a.symbol == 'O'and len(nl.get_neighbors(index)[0]) < 2:
        return(True)
    elif a.symbol in ['Si', 'Sn', 'Al', 'Ga', 'B'] and len(nl.get_neighbors(index)[0]) < 4:
        return(True)
    else:
        return(False)

def get_cap_ox(clust, index):
    ox_inds = [i for i in clust if i.symbol == 'O']
    ox_to_cap = [j for j in ox_inds if needs_neighb(index)]
    return(ox_to_cap)

def get_cap_si(clust, index):
    si_inds = [i for i in clust if i.symbol == 'Si']
    si_to_cap = [j for j in si_inds if needs_neighb(index)]
    return(si_to_cap)

def add_cap_ox(clust, index):

def add_cap_h(clust, index):
    new_clust = clust
    cap_inds = get_cap_ox(clust, index)
    for i in cap_inds:
        neighb = nl.get_neighbors(index)[0]

        new_clust.append(new_h)

    return(new_clust)

def cap_clust(clust):
    si_capped = add_cap_ox(clust)
    ox_capped = add_cap_h(si_capped)
    return(ox_capped)

# testing
if __name__ == '__main__':
    from ase.io import read
    b = read('BEA.cif')

