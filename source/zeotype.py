from typing import List, Dict, Tuple
from ase import Atoms, Atom
from ase.neighborlist import natural_cutoffs, NeighborList
from collections import defaultdict
import copy
import numpy as np

class Zeotype(Atoms):
    """
    This is a Zeotype class that inherits from Atoms. It represents a Zeolite.
    """

    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, silent: bool = False, zeolite_type: str = ''):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms, charges, scaled_positions,
                         cell, pbc, celldisp, constraint, calculator, info, velocities)

        self.zeolite_type = zeolite_type
        self.sites: List[str] = []
        self.clusters = []
        self.silent = silent
        self.neighbor_list = NeighborList(natural_cutoffs(self), bothways=True, self_interaction=False)
        self.neighbor_list.update(self)

    @staticmethod
    def build_from_atoms(a: Atoms, silent=False) -> "Zeotype":
        """
        :param a: Atoms object that Zeotype object will be used to create the Zeotype object
        :param silent: currently does nothing, but will set if output is displayed
        :return: Zeotype object created from Atoms object
        """
        new_zeolite = Zeotype(silent=silent)
        atom_zeolite = copy.deepcopy(a)
        atom_zeolite.__class__ = Zeotype
        new_zeolite.__dict__.update(atom_zeolite.__dict__)
        return new_zeolite

    def get_sites(self) -> List[str]:
        """
        :return: returns Zeotype sites
        """
        return self.sites

    def get_zeolite_type(self) -> List[str]:
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

        for atom in self:  # iterate through atom objects in zeotype
            # labels framework atoms
            if atom.symbol in ['Sn', 'Al', 'Si']:
                label = 'framework-%s' % atom.symbol
            elif atom.symbol in ['C', 'N']:
                label = 'adsorbate-%s' % atom.symbol
            elif atom.symbol in ['Cu', 'Ni', 'Fe', 'Cr']:
                label = 'extraframework-%s' % atom.symbol
            else:
                label = 'other'

            type_dict[label].append(atom.index)

        return type_dict

    def count_elements(self) -> Tuple[Dict['str', List[int]], Dict['str', int]]:
        """
        Stores the indices and counts the number of each element in the
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
        new_cluster = Cluster(self, index, size)
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
    def __init__(self, parent_zeotype: Zeotype, index: int, cluster_size: int):
        super().__init__()
        cluster_indices = self._get_cluster_indices(parent_zeotype, index, cluster_size)
        cluster_atoms = parent_zeotype[cluster_indices]

        new_self = copy.deepcopy(cluster_atoms)
        new_self.__class__ = Cluster
        self.__dict__.update(new_self.__dict__)

        self.parent_zeotype = parent_zeotype
        self.zeotype_to_cluster_index_map = \
            self._get_new_cluster_mapping(self.parent_zeotype, cluster_atoms, cluster_indices)

        self.neighbor_list = NeighborList(natural_cutoffs(self), bothways=True, self_interaction=False)
        self.neighbor_list.update(self)

    def append(self, atom):
        """Append atom to end."""
        self.extend(atom)

    def extend(self, other):
        """Extend atoms object by appending atoms from *other*."""
        if isinstance(other, Atom):
            other = Atoms([other])

        n1 = len(self)
        n2 = len(other)

        for name, a1 in self.arrays.items():
            a = np.zeros((n1 + n2,) + a1.shape[1:], a1.dtype)
            a[:n1] = a1
            if name == 'masses':
                a2 = other.get_masses()
            else:
                a2 = other.arrays.get(name)
            if a2 is not None:
                a[n1:] = a2
            self.arrays[name] = a

        for name, a2 in other.arrays.items():
            if name in self.arrays:
                continue
            a = np.empty((n1 + n2,) + a2.shape[1:], a2.dtype)
            a[n1:] = a2
            if name == 'masses':
                a[:n1] = self.get_masses()[:n1]
            else:
                a[:n1] = 0

            self.set_array(name, a)

    def cap_atoms(self):
        """each bare Si atom needs 4 Si atoms in oxygen and each of those oxygen needs two neighbors
        A lot of these clusters we add a bunch of H to the ends of the O because this a chemically plausable
        way to cap the oxygens. For example, when a zeolite is growing in solution it starts off as SiOH4
        1. Go through all Si in cluster and check to see if they have 4 bonds
        2.  If they don't have four bonds add an oxygen in a plausable direction
        3. go through all oxygens and find the ones that don't have two bonds
        4. If oxygens don't have two bonds add a hydrogen in a reasonable location (1 A distance)
        """
        bonds_needed = {'O': 2, 'Si': 4, 'Sn': 4, 'Al': 4, 'Ga': 4, 'B': 4}
        self.neighbor_list.update(self)  # check to see if it runs without this
        indices, count = self.count_elements()

        for si_index in indices['Si']:
            if self.needs_cap(si_index, bonds_needed['Si']):
                self.add_oxygen_cap(si_index, bonds_needed['Si'])
        for o_index in indices['O']:
            if self.needs_cap(o_index, bonds_needed['O']):
                self.add_hydrogen_cap(o_index)

    def needs_cap(self, atom_index: int, bonds_needed: int) -> bool:
        return len(self.neighbor_list.get_neighbors(atom_index)[0]) < bonds_needed

    def add_oxygen_cap(self, index, bonds_needed):
        while len(self.neighbor_list.get_neighbors(index)[0]) < bonds_needed:
            neighbor = self.neighbor_list.get_neighbors(index)[0][0]  # first index in the list of neighbor indicies
            direction = self.get_positions()[index] - self.get_positions()[neighbor]  # vector from neighbor to Si
            oxygen_pos = self.get_positions()[index] + (self.get_positions()[index] + direction) / np.linalg.norm(direction)
            new_oxygen = Atom('O', position=oxygen_pos)
            self.append(new_oxygen)
            self.neighbor_list = NeighborList(natural_cutoffs(self), bothways=True, self_interaction=False)
            self.neighbor_list.update(self)

    def add_hydrogen_cap(self, index):
        neighbor = self.neighbor_list.get_neighbors(index)[0][0]  # first index in the list of neighbor indicies
        direction = self.get_positions()[index] - self.get_positions()[neighbor]  # vector from neighbor to oxygen
        h_pos = self.get_positions()[index] + (self.get_positions()[index] + direction)/np.linalg.norm(direction)
        new_h = Atom('H', position=h_pos)
        self.append(new_h)


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
