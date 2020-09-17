from ase.io import read
from ase.neighborlist import NeighborList, natural_cutoffs
from ase import Atom
import numpy as np
'''
Takes position the adsorbate needs to be placed as a vector,
the site at which the donor atom must
face (as the index of that atom in the host lattice),
the atoms object for the adsorbate and the atoms object for the host lattice
'''
# needs to choose donor atom and then 1. point molecule along donor vec 2. move donor atom to zero position
# 3. try out different orientations and test them by moving donor to pos 4. add donor and ads in orientation
def place_ads(pos, host_ind, donor_ind, ads, host):
    '''
    :param pos: vector, the position to place adsorbate's donor atom
    :param host_ind: integer, index of site in host which adsorbate will be bound
    :param donor_ind: integer, index of donor atom on adsorbate
    :param ads: atoms object
    :param host: atoms object
    :param elements: list
    :return:
    '''
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
    ads2.translate([0.0, 0.0, 0.0]) #- ads2.get_positions()[donor_no]) # moves donor to zero position to prepare for orientation
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
    new_host = place_ads([5.6,8.3,15.2], 185, 1, ads, host) # example run with nice parameters
    view(new_host)
