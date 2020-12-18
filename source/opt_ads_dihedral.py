import zeotype  # import doesn't work
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.constraints import FixAtoms
from ase.constraints import FixBondLengths

def constrain_host(host, index):
    '''
    constrains a zeotype so only adsorbate, central atom, and adjacent oxygens can move.
    :param host: host zeotype object
    :param index: index of central atom in host zeotype
    :return: constrained host zeotype object
    '''

    host_ind_to_constrain = host.get_cluster_indices(host, index=index, max_neighbors=1)  # ind of framework atoms
    c = FixAtoms(host_ind_to_constrain)
    host.set_constraint(c)  # constrain everything except adsorbate, central atom, and adjacent oxygen
    return host  # should apply change to a copy??

def constrain_ads_bonds(host):
    '''
    Constrains the distance between atoms in adsorbate to adjacent adsorbate atoms
    :param host: host zeotype with adsorbate
    :return: host zeotype with constrained adsorbate
    '''
    atom_types = host.get_atom_types()
    ads_inds = [atom_types[i] for i in atom_types.keys() if 'adsorbate' in i]  # gets the indicies of adsorbate atoms

    nl = NeighborList(natural_cutoffs(host), self_interaction=False, bothways=True)
    nl.update(host)
    bonds_to_fix = []

    for ind in ads_inds:
        neigb_list = [j for j in nl.get_neighbors(ind)[0] if j in ads_inds]  # indicies of adsorbate atom neighbors
        for neighb_ind in neigb_list:
            # makes pairs of adjacent atom indicies, small index first
            if neighb_ind > ind:
                bond_pair = [ind, neighb_ind]
            else:
                bond_pair = [neighb_ind, ind]
            # adds bond pair to list iff not already there
            if bond_pair not in bonds_to_fix:
                bonds_to_fix.append(bond_pair)
            else:
                pass

    c = FixBondLengths(bonds_to_fix)
    host.set_constraint(c)  # should apply change to a copy??
    return host

# testing
if __name__ == '__main__':
    from ase.io import read
    from ase.visualize import view

    atoms = read('BEA.cif')
    host = Zeotype(atoms)
    a = constrain_host(host, 185)
    b = constrain_ads_bonds(a)
    view(b)











