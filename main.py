from source.cif_download import build_zeolite_from_code
from source.zeotype import Zeotype
from ase.visualize import view
import os
from ase.neighborlist import natural_cutoffs, NeighborList
import copy
#from ase.io import read
from ase.io import cif
import ase
import ase.data



if __name__ == '__main__':

    cif_path = os.path.join(os.getcwd(), 'source', 'BEA.cif')
    #zeotype = Zeotype.build_from_cif(cif_path)
    z = Zeotype.build_from_cif_with_labels(cif_path)
    z.add_cluster(3,10)
    z.add_cluster(30, 10)
    print(z.t_site_to_atom_indices)
    # make a list of avalible chem symbols
    # available_chem_symbols = copy.deepcopy(ase.data.chemical_symbols)
    # for symbol in set(zeotype.get_chemical_symbols()):
    #     pop_index = available_chem_symbols.index(symbol)
    #     item = available_chem_symbols.pop(pop_index)
    #



    #available_chem_symbols.pop()
    # print(len(zeotype.atom_sites_label))
    # print(len(zeotype.get_c
    # hemical_symbols()))
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

