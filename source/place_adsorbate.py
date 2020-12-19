from ase.neighborlist import NeighborList, natural_cutoffs
from placement import find_best_place
from ase.data import atomic_numbers, covalent_radii
from ase import Atom
import numpy as np

def pick_donor(ads):
    """finds atom in adsorbate most likely to bind to metal Heuristic: N > O > P
    > S > X > C > H :param ads: adsorbate atoms object :return: index of donor
    atom in adsorbate

    Args:
        ads:
    """

    donors = ['N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'C', 'H']
    ads_symbols = [k.symbol for k in ads]
    for symbol in donors:
        if symbol in ads_symbols:
            ads_ind = [s.index for s in ads if s.symbol == symbol][0] # picks first index of atom with specified symbol
            return(ads_ind)
    else:
        return ValueError("Cannot find donor atom in ads")

def pick_host_atom(host):
    """picks element with largest atomic number higher than 14 (Si), else, Al
    :param host: host atoms object :return: index of host atom

    Args:
        host:
    """

    atom_nums = host.get_atomic_numbers()
    max_num = max(atom_nums)
    symbol = Atom(max_num).symbol  # symbol of atom with highest atomic number
    host_ind = [i.index for i in host if i.symbol == symbol][0]  # picks first index of atom with specified symbol
    if max_num == 14:
        if 13 in atom_nums:
            al_ind = [j.index for j in host if j.symbol == 'Al'][0]
            return (al_ind)  # return index of an Al atom if no heavy atoms in host
        else:
            return ValueError("No heavy atoms other than Si in host")
    else:
        return(host_ind)

def pick_pos(ads, host, donor_ind, host_ind, radius=None, cutoff=None):
    """Finds a good position to add adsorbate to host :param ads: adsorbate
    atoms object :param host: host atoms object :param donor_ind: index of donor
    atom on adsorbate :param host_ind: index of host binding site :param radius:
    distance between host atom and donor atom :param cutoff: minimum distance
    donor atom must be from host atoms :return: vector of best position to place
    adsorbate

    Args:
        ads:
        host:
        donor_ind:
        host_ind:
        radius:
        cutoff:
    """
    donor_atom_symbol, host_atom_symbol = ads[donor_ind].symbol, host[host_ind].symbol
    donor_atom_number, host_atom_number = atomic_numbers[donor_atom_symbol], atomic_numbers[host_atom_symbol]
    donor_radius, host_radius = covalent_radii[donor_atom_number], covalent_radii[host_atom_number]

    if radius == None:
        radius = host_radius + donor_radius
    if cutoff == None:
        cutoff = host_radius

    pos = find_best_place(host, host_ind, radius, cutoff)
    return(pos)

def place_ads(ads, host, donor_ind=None, host_ind=None, pos=None):
    """Places adsorbate in host according to specified parameters :param pos:
    vector, the position to place adsorbate's donor atom :param host_ind:
    integer, index of site in host which adsorbate will be bound :param
    donor_ind: integer, index of donor atom on adsorbate :param ads: adsorbate
    atoms object :param host: host atoms object :param elements: list :return:
    atoms object with adsorbate in host

    Args:
        ads:
        host:
        donor_ind:
        host_ind:
        pos:
    """

    if donor_ind == None:
        donor_ind = pick_donor(ads)
    if host_ind == None:
        host_ind = pick_host_atom(host)
    if pos == None:
        pos = pick_pos(ads, host, donor_ind, host_ind)

    dummy_host = host + Atom('H', position=pos)  # add dummy hydrogen atom to get distances to host atoms
    vec = dummy_host.get_distance(-1, host_ind, mic=True, vector=True)
    nl = NeighborList(natural_cutoffs(ads), self_interaction=False, bothways=True)
    nl.update(ads)
    # gets neighbors of donor atom and adds the vectors from neighbor to donor
    # for most donor atoms this is roughly in the proper binding direction
    donor_neighbors = nl.get_neighbors(donor_ind)[0] # neighbor's index
    donor_vec = [0, 0, 0]
    for i in donor_neighbors:
        a = ads.get_distance(i, donor_ind, vector=True)
        donor_vec = donor_vec + a
    donor_vec = donor_vec/np.linalg.norm(donor_vec) # normalizes donor vec
    # rotate the ads into binding direction, move the ads to proper pos and combine
    ads2 = ads.copy() # to avoid making changes to orig adsorbate
    ads2.rotate(donor_vec, vec)
    ads2.translate(pos - ads2.get_positions()[donor_ind])
    host_new = host + ads2
    return (host_new)

# testing
if __name__ == '__main__':
    from ase.io import read
    from ase.visualize import view
    from ase.build import molecule
    host = read('BEA.cif')
    host[185].symbol = 'Sn'  # for visualization
    ads = molecule('CH3OH')
    #new_host = place_ads(ads, host, 1, 185, [5.6,8.3,15.2]) # example run with nice parameters
    new_host = place_ads(ads, host)
    view(new_host)