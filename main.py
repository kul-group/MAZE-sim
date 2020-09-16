from source.zeotype import Zeotype
from ase.visualize import view
import os
from ase.neighborlist import natural_cutoffs, NeighborList
import copy
#from ase.io import read
from ase.io import cif
import ase

def my_read_cif(fileobj, store_tags=False, primitive_cell=False,
             subtrans_included=True, fractional_occupancies=True,
             reader='ase'):

    blocks = ase.io.cif.parse_cif(fileobj, reader)
    tags = blocks[0][1]
    # Find all CIF blocks with valid crystal data
    atoms = ase.io.cif.tags2atoms(tags, store_tags, primitive_cell,
                               subtrans_included,
                               fractional_occupancies=fractional_occupancies)

    return atoms
    #         if store_tags:
    #             try:
    #                 atoms.set_tags(atoms.info['_atom_Site_type_symbol'])
    #             except KeyError:
    #                 pass
    #         return atoms
    #         images.append(atoms)
    #     except KeyError:
    #         pass
    # for atoms in images[0]:
    #     yield atoms


def my_read_cif_2(fileobj, store_tags=False, primitive_cell=False,
             subtrans_included=True, fractional_occupancies=True,
             reader='ase'):
    blocks = ase.io.cif.parse_cif(fileobj, reader)
    # Find all CIF blocks with valid crystal data
    images = []
    for name, tags in blocks:
        try:
            atoms = ase.io.cif.tags2atoms(tags, store_tags, primitive_cell,
                               subtrans_included,
                               fractional_occupancies=fractional_occupancies)
            images.append(atoms)
        except KeyError:
            pass
    for atoms in images:
        yield atoms

def build_zeolite_from_cif(fileobj):
    atoms_gen = my_read_cif_2(fileobj, store_tags=True)
    #atoms = my_read_cif_2(fileobj, store_tags=True) # next(atoms_gen) # this might be good enough
    atoms = next(atoms_gen)
    my_zeotype = Zeotype(atoms)
    try:
        my_zeotype.atom_sites_label = atoms.info['_atom_site_label']
    except KeyError:
        print("No atom site labels in CIF")
        pass
    return my_zeotype

if __name__ == '__main__':
    cif_path = os.path.join(os.getcwd(), 'source', 'BEA.cif')
    zeotype = Zeotype.build_zeolite_from_cif(cif_path)
    #b = read(cif_path, store_tags=True)
    #c = my_read_cif(cif_path, store_tags=True)
    #zeotype = Zeotype(c)
    #zeotype.atom_sites_label = c.info['_atom_site_label']
    #zeotype = build_zeolite_from_cif(cif_path)
    print(zeotype.atom_sites_label)
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

