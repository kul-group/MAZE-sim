from typing import List, Dict, Tuple
from ase import Atoms, Atom
from ase.neighborlist import natural_cutoffs, NeighborList
from collections import defaultdict
import copy
import numpy as np
from ase.io import cif
import ase


class Zeotype(Atoms):
    """
    This is a Zeotype class that inherits from Atoms. It represents a Zeolite.
    """
    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, silent: bool = False, zeolite_type: str = '',
                 atom_site_labels=None):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms, charges, scaled_positions,
                         cell, pbc, celldisp, constraint, calculator, info, velocities)

        self.zeolite_type = zeolite_type
        self.sites: List[str] = []
        self.clusters = []
        self.silent = silent
        self.neighbor_list = NeighborList(natural_cutoffs(self), bothways=True, self_interaction=False)
        self.neighbor_list.update(self)
        self.atom_sites_label = atom_site_labels

    @classmethod
    def build_zeolite_from_cif(cls, fileobj):
        atoms_gen = cls.read_cif_with_info(fileobj, store_tags=True)
        atoms = next(atoms_gen)
        zeotype = Zeotype(atoms)
        try:
            zeotype.atom_sites_label = atoms.info['_atom_site_label']
        except KeyError:
            print("No atom site labels in CIF")
            pass
        return zeotype

    @staticmethod
    def read_cif_with_info(fileobj, store_tags=False, primitive_cell=False,
                      subtrans_included=True, fractional_occupancies=True,
                      reader='ase') -> Atoms:
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


    def get_sites(self) -> List[str]:
        """
        :return: returns Zeotype sites
        """
        return self.sites

    def get_zeolite_type(self) -> str:
        """
        :return: returns zeotype sites as a list
        """
        return self.zeolite_type

    def get_atom_types(self) -> Dict[str, List[int]]:
        """
        :return: Returns a dictionary of atom types where the key consists of the atom category
        (framework, absorbate, extraframework or other) followed by -atom chemical symbol. For
        example a Sn atom is a framework atom so the key is "framework-Sn". The value of the
        returned dictionary is a list of the indices of all of the atoms that belong to the
        category.
        """

        type_dict: Dict[str, List[int]] = defaultdict(list)  # default dict sets the default value of dict to []

        nl = self.neighbor_list
        for atom in self:  # iterate through atom objects in zeotype
            # labels framework atoms
            if atom.symbol in ['Sn', 'Al', 'Si']:
                label = 'framework-%s' % atom.symbol
            elif atom.symbol in ['C', 'N']:
                label = 'adsorbate-%s' % atom.symbol
            elif atom.symbol in ['Cu', 'Ni', 'Fe', 'Cr']:
                label = 'extraframework-%s' % atom.symbol
            elif atom.symbol == 'O':
                neigh_ind = nl.get_neighbors(atom.index)
                if len(neigh_ind) == 2:
                    neigh_sym1, neigh_sym2 = self[neigh_ind[0][0]].symbol, self[neigh_ind[0][1]].symbol
                    # assign framework oxygens
                    if neigh_sym1 in ['Sn', 'Al', 'Si'] and neigh_sym2 in ['Sn', 'Al', 'Si']:
                        label = 'framework-%s' % atom.symbol
                    # assign bound adsorbate oxygen
                    elif neigh_sym1 in ['H', 'C', 'N'] and neigh_sym2 in ['Sn', 'Al', 'Si']:
                        label = 'bound-adsorbate-%s' % atom.symbol
                    elif neigh_sym2 in ['H', 'C', 'N'] and neigh_sym1 in ['Sn', 'Al', 'Si']:
                        label = 'bound-adsorbate-%s' % atom.symbol
                    # assign rest as adsorbate oxygen
                    else:
                        label = 'adsorbate-%s' % atom.symbol
                else:  # is this correct
                    label = 'other'
            else:
                label = 'other'

            type_dict[label].append(atom.index)

        return dict(type_dict)

    def count_elements(self) -> Tuple[Dict['str', List[int]], Dict['str', int]]:
        """
        :return:
        """
        indices: Dict['str', List['int']] = defaultdict(list)  # indices of the elements grouped by type
        count: Dict['str', int] = defaultdict(lambda: 0)  # number of elements of each type
        for atom in self:
            element = atom.symbol
            indices[element].append(atom.index)
            count[element] += 1
        return indices, count

    @staticmethod
    def count_atomtypes(atomtype_list) -> Tuple[Dict['str', List[int]], Dict['str', int]]:
        """
        Counts the number of different atoms of each type in a list of atom symbols
        :param atomtype_list: A list of atom chemical symbols
        :return: A tuple of two dictionaries the first containing a mapping between the chemical
        symbol and the indices in the symbol and the second containing a mapping between the
        chemical symbol and the number of symbols of that type in the input list
        """
        indices: Dict['str', List['int']] = defaultdict(list)  # indices of the elements grouped by type
        count: Dict['str', int] = defaultdict(lambda: 0)  # number of elements of each type
        for i, element in enumerate(atomtype_list):
            indices[element].append(i)
            count[element] += 1
        return indices, count  # TODO: Combine with count_elements method

    def add_cluster(self, index: int, size: int) -> int:
        """
        Generates a Cluster of atoms around the specified index. The number of atoms in the cluster
        is given by the size parameter.
        :param index: index of the central atom in the cluster
        :param size: number of atoms in the final cluster
        :return: index of the cluster in the zeotype cluster array
        """
        new_cluster = Cluster.build_from_zeolite(self, index, size)
        self.clusters.append(new_cluster)
        return len(self.clusters) - 1  # returns final index in the clusters list

    def remove_cluster(self, index: int):
        """
        :param index: Index of cluster to remove from zeotype list
        :return: None
        """
        self.clusters.pop(index)

    def integrate_cluster(self, cluster_index: int):
        cluster = self.clusters[cluster_index]
        index_map = cluster.zeotype_to_cluster_index_map

        for key, value in index_map.items():
            # same unit cell same oxygen atom positions
            self[key].symbol = cluster[value].symbol
            self[key].position = cluster[value].position
            self[key].tag = cluster[value].tag
            self[key].momentum = cluster[value].momentum
            self[key].mass = cluster[value].mass
            self[key].magmom = cluster[value].magmom
            self[key].charge = cluster[value].charge


class Cluster(Zeotype):  # TODO include dynamic inheritance

    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, silent: bool = False, zeolite_type: str = '',
                 parent_zeotype=None, zeotype_to_cluster_index_map=None, neighbor_list=None):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms,
                         charges, scaled_positions, cell, pbc, celldisp, constraint,
                         calculator, info, velocities, silent, zeolite_type)

        self.parent_zeotype = parent_zeotype
        self.zeotype_to_cluster_index_map = zeotype_to_cluster_index_map
        self.neighbor_list = neighbor_list

    @staticmethod
    def build_from_zeolite(parent_zeotype: Zeotype, index: int, cluster_size: int) -> "Cluster":
        cluster_indices = Cluster._get_cluster_indices(parent_zeotype, index, cluster_size)
        cluster_atoms = parent_zeotype[cluster_indices]
        new_cluster = Cluster(cluster_atoms)
        new_cluster.parent_zeotype = parent_zeotype
        new_cluster.zeotype_to_cluster_index_map = \
            new_cluster._get_new_cluster_mapping(new_cluster.parent_zeotype, cluster_atoms, cluster_indices)

        new_cluster.neighbor_list = NeighborList(natural_cutoffs(new_cluster), bothways=True, self_interaction=False)
        new_cluster.neighbor_list.update(new_cluster)

        return new_cluster

    def cap_atoms(self, cap_atoms_dict=None, bonds_needed=None, verbose=False):
        """each bare Si atom needs 4 Si atoms in oxygen and each of those oxygen needs two neighbors
        A lot of these clusters we add a bunch of H to the ends of the O because this a chemically plausable
        way to cap the oxygens. For example, when a zeolite is growing in solution it starts off as SiOH4
        1. Go through all Si in cluster and check to see if they have 4 bonds
        2.  If they don't have four bonds add an oxygen in a plausable direction
        3. go through all oxygens and find the ones that don't have two bonds
        4. If oxygens don't have two bonds add a hydrogen in a reasonable location (1 A distance)
        """
        if cap_atoms_dict is None:
            cap_atoms_dict = self.build_cap_atoms_dict(bonds_needed=bonds_needed)
        if verbose:
            print('atom caps: {symbol, [position arrays]}', cap_atoms_dict)

        # append cap atoms self
        for symbol, pos_list in cap_atoms_dict.items():
            for pos in pos_list:
                self.append(Atom(symbol, position=pos))

    def build_cap_atoms_dict(self, bonds_needed=None):
        """
        Builds a dictionary of the cap atom positions
        :param bonds_needed: a dict mapping atom symbol to number of bonds needed
        :return: a dictionary of atom symbols and positions of cap atoms
        """
        cap_atoms_dict: Dict[str, List[int]] = defaultdict(list)
        if bonds_needed is None:
            bonds_needed = {'O': 2, 'Si': 4, 'Sn': 4, 'Al': 4, 'Ga': 4, 'B': 4}
        indices, count = self.count_elements()
        for si_index in indices['Si']:
            if self.needs_cap(si_index, bonds_needed['Si']):
                pos = self.get_oxygen_cap_pos(si_index, bonds_needed['Si'])
                cap_atoms_dict['O'].append(pos)
        for o_index in indices['O']:
            if self.needs_cap(o_index, bonds_needed['O']):
                pos = self.get_hydrogen_cap_pos(si_index)
                cap_atoms_dict['H'].append(pos)

        return dict(cap_atoms_dict)

    def needs_cap(self, atom_index: int, bonds_needed: int) -> bool:
        return len(self.neighbor_list.get_neighbors(atom_index)[0]) < bonds_needed

    def get_oxygen_cap_pos(self, index, bonds_needed):
        # while len(self.neighbor_list.get_neighbors(index)[0]) < bonds_needed:
        neighbor = self.neighbor_list.get_neighbors(index)[0][0]  # first index in the list of neighbor indicies
        direction = self.get_positions()[index] - self.get_positions()[neighbor]  # vector from neighbor to Si
        oxygen_pos = self.get_positions()[index] + (self.get_positions()[index] + direction) / np.linalg.norm(direction)
        return oxygen_pos
        # self.neighbor_list = NeighborList(natural_cutoffs(self), bothways=True, self_interaction=False)
        # self.neighbor_list.update(self)

    def get_hydrogen_cap_pos(self, index):
        neighbor = self.neighbor_list.get_neighbors(index)[0][0]  # first index in the list of neighbor indicies
        direction = self.get_positions()[index] - self.get_positions()[neighbor]  # vector from neighbor to oxygen
        hydrogen_pos = self.get_positions()[index] + (self.get_positions()[index] + direction) / np.linalg.norm(
            direction)
        return hydrogen_pos

    @staticmethod
    def _get_cluster_indices(zeolite, index: int, size: int):
        nl = NeighborList(natural_cutoffs(zeolite), self_interaction=False, bothways=True)
        nl.update(zeolite)

        cluster_indices = set()
        new_cluster_indices = set([index])

        while True:
            current_cluster_indices = set()
            for cluster_index in new_cluster_indices:
                cluster_indices.add(cluster_index)
                if len(cluster_indices) >= size:
                    return list(cluster_indices)
                for new_index in nl.get_neighbors(cluster_index)[0]:
                    current_cluster_indices.add(new_index)
            new_cluster_indices = current_cluster_indices

    @staticmethod
    def _get_new_cluster_mapping(zeolite, cluster, indices: List[int]):
        cluster_position_index_map = {}
        for atom in cluster:
            cluster_position_index_map[str(atom.position)] = atom.index

        zeotype_index_position_map = {}
        for i in indices:
            zeotype_index_position_map[i] = str(zeolite[i].position)
        zeotype_to_cluster_index_map = {}

        for key in zeotype_index_position_map.keys():
            zeotype_to_cluster_index_map[key] = \
                cluster_position_index_map[zeotype_index_position_map[key]]

        return zeotype_to_cluster_index_map


# testing
if __name__ == '__main__':
    from ase.io import read

    b = read('BEA.cif')
    z = Zeotype(b)
    z.add_cluster(35, 3)
    print([a for a in z.clusters[0]])
    z.integrate_cluster(0)
