from source.cif_download import build_zeolite_from_code
from source.zeotype import Zeotype
from ase.visualize import view
import os
from ase.neighborlist import natural_cutoffs, NeighborList
import copy
#from ase.io import read
from ase.io import cif
import ase


if __name__ == '__main__':

    cif_path = os.path.join(os.getcwd(), 'source', 'BEA.cif')
    zeotype = Zeotype.build_from_cif(cif_path)
    print(len(zeotype.atom_sites_label))
    print(len(zeotype.get_chemical_symbols()))
    # for a in zeotype:
    #    print(a)
    #print(c.info['_atom_site_label'])
    #print(c.get_tags())
    #d = ase.io.cif.
    #a_list = []

    #for a in c:
    #    print(type(a))
    #    print(a)
    #     a_list.append(a)
    # my_atoms = ase.Atoms(c)
    # print(my_atoms.get_tags())
    #print(c)
    #print(c.get_tags())
    # z = Zeotype(b)
    # index = z.add_cluster(187, 7)
    # view(z)
    #view(z.clusters[index])
    # w = copy.deepcopy(z)
    # w.clusters[0].cap_atoms(verbose=True)
    # w.clusters[0].set_pbc([True, True, True])
    # w.clusters[0].wrap()
    #view(w.clusters[0])
    # nl = NeighborList(natural_cutoffs(z), self_interaction=False, bothways=True)
    # nl.update(z)
    #print(z.clusters[0].count_elements())
    #z.clusters[0]._get_cap_O_pos(4)
    #print(type(z[0].position))

