from source.zeotype import Zeotype
from ase.visualize import view
import os
from ase.neighborlist import natural_cutoffs, NeighborList
import copy
from ase.io import read

if __name__ == '__main__':
    cif_path = os.path.join(os.getcwd(), 'source', 'BEA.cif')
    b = read(cif_path)
    z = Zeotype(b)
    index = z.add_cluster(187, 7)
    view(z)
    view(z.clusters[index])
    w = copy.deepcopy(z)
    w.clusters[0].cap_atoms()
    view(w.clusters[0])
    # nl = NeighborList(natural_cutoffs(z), self_interaction=False, bothways=True)
    # nl.update(z)
    #print(z.clusters[0].count_elements())
    #z.clusters[0]._get_cap_O_pos(4)
    #print(type(z[0].position))

