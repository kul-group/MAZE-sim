from ase.neighborlist import NeighborList, natural_cutoffs
from ase import Atom
import numpy as np

# needs to:
# 1. choose host/donor atom
# 2. point molecule along donor vec
# 3. move donor atom to correct position OR choose a good position
# 4. try out different orientations and test them by moving donor to pos
# 5. add donor and ads in orientation
# 6. provide option to automatically choose placement position
# 7. provide option to automatically shape the adsorbate: constrain bondlength and angles then minimize contact w/ pore

def pick_donor(ads):
    # finds atom in adsorbate most likely to bind to metal
    # Heristic: N > O > P > S > X > C > H
    donors = ['N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'C', 'H']
    ads_symbols = [k.symbol for k in ads]
    for symbol in donors:
        if symbol in ads_symbols:
            ads_ind = [s.index for s in ads if s.symbol == symbol][0]
            return(ads_ind)
    else:
        return ValueError("Cannot find donor atom in ads")

def pick_host_atom(host):
    # returns host atom index to add an adsorbate to
    # picks element with largest atomic number higher than 14 (Si) or Al
    atom_nums = host.get_atomic_numbers()
    max_num = max(atom_nums)
    symbol = Atom(max_num).symbol # symbol of atom with highest atomic number
    host_ind = [i.index for i in host if i.symbol == symbol][0]
    if max_num == 14:
        if 13 in atom_nums:
            al_ind = [j.index for j in host if j.symbol == 'Al'][0]
            return (al_ind)
        else:
            return ValueError("No heavy atoms other than Si in host")
    else:
        return(host_ind)



def pick_pos(host, host_ind):
    # chooses best position to add adsorbate


def place_ads(pos, ads, host, host_ind=None, donor_ind=None,):
    '''
    :param pos: vector, the position to place adsorbate's donor atom
    :param host_ind: integer, index of site in host which adsorbate will be bound
    :param donor_ind: integer, index of donor atom on adsorbate
    :param ads: atoms object
    :param host: atoms object
    :param elements: list
    :return:
    '''

    if host_ind == None:
        host_ind = pick_host_atom(host)
    if donor_ind == None:
        donor_ind = pick_donor(ads)
    dummy_atom = Atom('H', position=pos)
    dummy_host = host + dummy_atom
    vec = dummy_host.get_distance(-1, host_ind, mic=True, vector=True)
    host_indicies, ads_indicies = [k.index for k in host], [j.index for j in ads]
    nl = NeighborList(natural_cutoffs(ads), self_interaction=False, bothways=True)
    nl.update(ads)
# gets neighbors of donor atom and adds the vectors from neighbor to donor
# for these atoms this is roughly in the proper binding direction
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
    new_host = place_ads([5.6,8.3,15.2], ads, host) # example run with nice parameters
    view(new_host)