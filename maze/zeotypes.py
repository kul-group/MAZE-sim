import copy
from collections import defaultdict
from typing import List, Dict, Tuple, Iterable, Optional, Union

import ase
import ase.data
import numpy as np
from ase import Atoms
from ase.io import cif
from ase.neighborlist import natural_cutoffs, NeighborList
from ase.visualize import view

from maze.adsorbate import Adsorbate
from maze.index_mapper import IndexMapper


class Silanol():
    """
    This Silanol class represents a Silanol Group (Si-O-H) It can be used to find Silanol nests in the Zeolite
    """

    def __init__(self, parent_zeotype: 'Zeotype', Si_index: int, O_index: int, H_index: int,
                 Si_neighbor_list: List[int]):
        """
       Create a Silanol object
        :param parent_zeotype: The parent zeolite, where the silanol group is located
        :param Si_index: the index of the Si atom
        :param O_index: the index of the O atom
        :param H_index: the index of the H atom
        :param Si_neighbor_list: The neighbor list of the Si atom, useful for other algos
        """
        self.parent_zeotype = parent_zeotype
        self.Si_index = Si_index
        self.O_index = O_index
        self.H_index = H_index
        self.Si_neighbor_list = Si_neighbor_list

    def __str__(self) -> str:
        """
        This ensures the str rep of ta Silanol object is nice
        :return: str rep of object
        """
        return f'Si: {self.Si_index} O: {self.O_index} H: {self.H_index}'

    def __repr__(self) -> str:
        """
        Makes repr same as str
        :return: str rep of object
        """
        return self.__str__()


class Zeotype(Atoms):
    """
    A class that inherits from ase.Atoms, which represents an unmodified maze. If a maze is built from a cif
    file from the zeolite structure database, then the atoms are tagged according to their unique site and the
    dictionaries site_to_atom_indices and atom_indices_to_site are filled. If not these dictionaries are set to none.
    This class contains a bunch of static methods for identifying different types of atoms in the zeolite as well as
    as methods to build an imperfect maze from a parent Zeotype class. Imperfect zeotypes have additional
    functionality and are dependent on a parent maze class.
    """

    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, site_to_atom_indices=None, atom_indices_to_site=None,
                 additions=None, _is_zeotype=True):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms, charges, scaled_positions,
                         cell, pbc, celldisp, constraint, calculator, info, velocities)

        # To Be compatible with ASE's atoms object, this code has to have the functionality to build from
        # an atom's object, a Zeotype object, or a sub class of a Zeotype object
        # The following if statements take care of this functionality depending on
        # if the object being built is a maze or if symbols is a maze (there are four unique paths)

        if isinstance(symbols, Zeotype):  # if symbols is a maze or maze subclass
            self.additions = copy.deepcopy(symbols.additions)
            if _is_zeotype:  # if the object being built is a maze
                self.site_to_atom_indices = symbols.site_to_atom_indices
                self.atom_indices_to_site = symbols.atom_indices_to_site
                self.name = 'parent'  # must be parent to agree with index mapper
                self.index_mapper = IndexMapper(self.get_indices(self))
                self.parent_zeotype = self

            else:  # if the object being built is a subtype of maze
                self.parent_zeotype = symbols.parent_zeotype
                self.site_to_atom_indices = None
                self.atom_indices_to_site = None
                self.index_mapper = symbols.index_mapper
                self.name = self.index_mapper.get_unique_name(type(self).__name__)  # use name
                self.index_mapper.add_name(self.name, symbols.name, self._get_old_to_new_map(symbols, self))

        else:  # if symbols is not a maze or maze child class
            if _is_zeotype:
                self.name = 'parent'  # must be parent for code to work properly
                self.index_mapper = IndexMapper(self.get_indices(self))
                self.parent_zeotype = self
            else:
                self.name = None
                self.index_mapper = None
                self.parent_zeotype = None

            self.site_to_atom_indices = site_to_atom_indices
            self.atom_indices_to_site = atom_indices_to_site
            self.additions = copy.deepcopy(additions) if additions else defaultdict(list)

        self.update_nl()

    @staticmethod
    def get_available_symbols(atom_list: List[str]) -> List[str]:
        """
        :param atom_list: A list of atom symbols to be excluded from returned list
        :return: a list of all possible atom symbols that do not include items from
        the original list
        """
        available_chem_symbols = copy.deepcopy(ase.data.chemical_symbols)
        for symbol in set(atom_list):
            pop_index = available_chem_symbols.index(symbol)
            available_chem_symbols.pop(pop_index)

        return available_chem_symbols

    @classmethod
    def build_from_cif_with_labels(cls, filepath: str) -> "Zeotype":
        """
        Takes a filepath/fileobject of a cif file and returns a Zeotype or child class with T-sites
        labeled as specified in the cif file. The dictionaries for T site labels are also filled in
        these map from the T-site to a list of atom indices. The atom indices to t sites maps atom
        indices to T-sites.

        :param filepath: Filepath of the cif file
        :return: Zeotype with labeled sites
        """
        atoms, site_to_atom_indices, atom_indices_to_site = cls._read_cif_note_sites(filepath)
        zeotype = cls(atoms)
        zeotype.site_to_atom_indices = site_to_atom_indices
        zeotype.atom_indices_to_site = atom_indices_to_site
        return zeotype

    @staticmethod
    def _read_cif_note_sites(fileobj, store_tags=False, primitive_cell=False,
                             subtrans_included=True, fractional_occupancies=True,
                             reader='ase') -> Tuple[Atoms, Dict[str, int], Dict[int, str]]:
        """
        The helper function used by build_from_cif_with_labels. This loads a CIF file
        using ase.io.cif.parse_cif and then finds the T-site information. After finding
        the T-site information it then replaces each atom in the T-site with another
        atom type not found in the cif file. The crystal is then generated and then
        the resulting atoms object atoms are replaced by the original atoms. The mapping
        between the T-site names and the atom objects in the final atoms object are
        recored in two dictionaries.

        :param fileobj: CIF file location
        :param store_tags: store the tags in resulting atoms object
        :param primitive_cell: An option for the reader
        :param subtrans_included: An option for the reader
        :param fractional_occupancies: an option for the final crystal
        :param reader: the reader used (ase works others have not been tested)
        :return: atoms, site_to_atom_indices, atom_indices_to_site
        """
        blocks = ase.io.cif.parse_cif(fileobj, reader)  # read CIF file into dictionary
        b_dict = blocks[0][1]  # get the dictionary from CIF file
        # find the atoms that are T sites
        # t_indices = [i for i in range(0, len(b_dict["_atom_site_label"])) if 'T' in b_dict["_atom_site_label"][i]]

        # replace elements with replacement symbol
        element_to_T_site = {}
        sym_to_original_element = {}
        possible_syms = Zeotype.get_available_symbols(
            b_dict["_atom_site_type_symbol"])  # only get symbols not in CIF file
        for i in range(len(b_dict["_atom_site_label"])):  # label all atoms
            sym = possible_syms.pop()
            sym_to_original_element[sym] = b_dict["_atom_site_type_symbol"][i]
            b_dict["_atom_site_type_symbol"][i] = sym
            element_to_T_site[sym] = b_dict["_atom_site_label"][i]

        site_to_atom_indices = defaultdict(list)
        atom_indices_to_site = {}

        images = []
        atoms = None
        for name, tags in blocks:
            atoms = ase.io.cif.tags2atoms(tags, store_tags, primitive_cell,
                                          subtrans_included,
                                          fractional_occupancies=fractional_occupancies)
            images.append(atoms)

            for i in range(len(atoms)):  # replace substitute atoms with original Si and get indices of replaced atoms
                if atoms[i].symbol in element_to_T_site.keys():
                    sym = atoms[i].symbol
                    key = element_to_T_site[sym]
                    site_to_atom_indices[key].append(i)
                    atom_indices_to_site[i] = key
                    atoms[i].tag = ase.data.atomic_numbers[sym]
                    atoms[i].symbol = sym_to_original_element[sym]

        return atoms, dict(site_to_atom_indices), atom_indices_to_site

    def get_imperfect_zeolite(self) -> 'ImperfectZeotype':
        """
        :return: An imperfect Zeotype constructed from the current Zeotype
        """
        return ImperfectZeotype(self)

    def update_nl(self, mult: int = 1) -> None:
        """
        Builds and updates neighborlist
        :param mult: The mult (multiply) parameter for natural cutoffs (Default 1)
        :return: None
        """
        self.neighbor_list = NeighborList(natural_cutoffs(self, mult=mult), bothways=True, self_interaction=False)
        self.neighbor_list.update(self)

    def get_hetero_atoms(self, hetero_atoms_list: Optional[List[str]] = None) -> List[int]:
        """
        :return: Returns a list of all of the hetero-atoms in the maze
        """
        if not hetero_atoms_list:
            hetero_atoms_list = ['Sn', 'Hf', 'Zr', 'Ge', 'Ti']

        indices_list = []
        for atom in self:
            if atom.symbol in hetero_atoms_list:
                indices_list.append(atom.index)

        return indices_list

    def get_atom_types(self) -> Dict[str, List[int]]:
        """
        :return: Returns a dictionary of atom types where the key consists of the atom category
        (framework, adsorbate, extraframework or other) followed by -atom chemical symbol. For
        example a Sn atom is a framework atom so the key is "framework-Sn". The value of the
        returned dictionary is a list of the indices of all of the atoms that belong to the
        category.
        """

        type_dict: Dict[str, List[int]] = defaultdict(list)  # default dict sets the default value of dict to []

        nl = self.neighbor_list
        for atom in self:  # iterate through atom objects in maze
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
        :return: a dictionary where the key is the element symbol and the value is the number in the maze
        """
        indices: Dict['str', List['int']] = defaultdict(list)  # indices of the elements grouped by type
        count: Dict['str', int] = defaultdict(lambda: 0)  # number of elements of each type
        for atom in self:
            element = atom.symbol
            indices[element].append(atom.index)
            count[element] += 1
        return indices, count

    @staticmethod
    def count_atomtypes(atomtype_list: List[str]) -> Tuple[Dict['str', List[int]], Dict['str', int]]:
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

    def get_cluster(self, index: int = 0, max_size: int = 0, max_neighbors: int = 0, cluster_indices=None) -> int:
        """
        Generates a Cluster of atoms around the specified index. The number of atoms in the cluster
        is given by the size parameter.

        :param cluster_indices: Indices of cluster to make
        :param max_size: max size of cluster
        :param index: index of the central atom in the cluster
        :param max_neighbors: number of neighbors from the host atom in the final cluster
        :return: index of the cluster in the maze cluster array
        """
        if cluster_indices is None:
            cluster_indices = Cluster.get_cluster_indices(self, index, max_size, max_neighbors)

        new_cluster = Cluster(self)
        indices_to_delete = self.get_indices_compliment(self, cluster_indices)
        new_cluster = new_cluster.delete_atoms(indices_to_delete)
        od = OpenDefect.build_from_indices(self, cluster_indices)
        # iz = self.get_imperfect_zeolite()
        # iz.delete_atoms(cluster_indices)
        # self.clusters.append(new_cluster)
        return new_cluster, od

    def find_silanol_groups(self) -> List[Silanol]:
        """
        Finds all of the silanol groups in the Zeotype

        :return: A list of Silanol groups
        """
        silanol_list = []
        for atom in self:
            if atom.symbol == 'Si':
                for neighbor_index in self.neighbor_list.get_neighbors(atom.index)[0]:
                    if self[neighbor_index].symbol == 'O':
                        for next_neighbor_index in self.neighbor_list.get_neighbors(self[neighbor_index].index)[0]:
                            if self[next_neighbor_index].symbol == 'H':
                                # found one!
                                silanol = Silanol(self, atom.index, neighbor_index, next_neighbor_index,
                                                  self.neighbor_list.get_neighbors(atom.index)[0])
                                silanol_list.append(silanol)
        return silanol_list

    def find_silanol_nest_T_sites(self) -> List[int]:
        """
        Finds all of the T sites that are in silanol nests

        :return: A list of T sites in silanol nests
        """
        sites_list = []
        sil_list = self.find_silanol_groups()
        if self.atom_indices_to_site is None:

            for sil in sil_list:
                sites_list.append(sil.Si_index)
                for i in sil.Si_neighbor_list:
                    if 'Sn' == self[i].symbol:
                        sites_list.append(i)
        else:
            for sil in sil_list:

                if 'T' in self.atom_indices_to_site(sil.Si_index):
                    sites_list.append(sil.Si_index)
                else:
                    for index in sil.Si_neighbor_list:
                        if 'Sn' == self[index].symbol:
                            sites_list.append(index)

        return sites_list

    @staticmethod
    def get_indices_compliment(zeotype: 'Zeotype', indices: Iterable[int]) -> List[int]:
        """
        Gets the compliment of indices in a Zeotype

        :param zeotype: Zeotype containing all indices
        :param indices: Indices to get the compliment of
        :return: Compliment of indices
        """
        return list(set([a.index for a in zeotype]) - set(indices))

    @staticmethod
    def get_indices(atoms_object: Atoms) -> List[int]:
        """
        Get the indices in an atoms object

        :param atoms_object: Atoms object to get Indices of
        :return: List of indices in Atoms object
        """
        return [a.index for a in atoms_object]

    def get_site_type(self, index: int) -> str:
        """
        Get the idenity of a site

        :param index: Index of site in self
        :return: Label for the site (comes from CIF) file
        """
        assert self.parent_zeotype.atom_indices_to_site is not None, 'atom_indices_to_site is None, cannot get label'
        pz_index = self.index_mapper.get_index(self, self.parent_zeotype, index)
        return self.parent_zeotype.atom_indices_to_site[pz_index]

    @staticmethod
    def _get_old_to_new_map(old: Atoms, new: Atoms) -> Dict[int, int]:
        """
        Get the index mapping between old and new self.
        his matching is done by position, so it is essential that the atom
        positions have not changed. It is also essential that new <= old

        :return: A mapping between old and new
        """
        new_index_to_position_map = {}
        for atom in new:
            new_index_to_position_map[atom.index] = str(atom.position)

        old_position_to_index_map = {}
        for atom in old:
            old_position_to_index_map[str(atom.position)] = atom.index

        old_to_new_map = {}

        for new_index, pos in new_index_to_position_map.items():
            try:
                old_to_new_map[old_position_to_index_map[pos]] = new_index
            except KeyError:
                continue

        return old_to_new_map

    def __del__(self) -> None:
        if self.index_mapper is not None:
            self.index_mapper.delete_name(self.name)


class ImperfectZeotype(Zeotype):
    """
    This imperfect maze inherits from Zeotype and requires a parent maze to function properly.
    This represents an imperfect Zeotype such as maze with some modification like an adsorbate addition
    or the removal a cluster.
    """

    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, site_to_atom_indices=None, atom_indices_to_site=None,
                 additions=None):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms,
                         charges, scaled_positions, cell, pbc, celldisp, constraint,
                         calculator, info, velocities, site_to_atom_indices, atom_indices_to_site,
                         additions, _is_zeotype=False)

    def register_self(self, source: Zeotype) -> None:
        """
        This method registers the current ImperfectZeotype with the index_mapper by using the common
        positions between the atoms in self and source.

        :param source: Source Zeotype or subclass that was used to build current ImperfectZeotype
            for the mapping to work correctly, the atom positions must be identical. This
        :return: None
        """
        name = self.index_mapper.get_unique_name(self.__name__)
        self.index_mapper.add_name(name, source.name, self._get_old_to_new_map(source, self))

    def build_cap_atoms(self, cap_atoms_dict: Dict[str, np.array]) -> Atoms:
        symbol_list = []
        position_list = []
        for symbol, pos_list in cap_atoms_dict.items():
            for pos in pos_list:
                symbol_list.append(symbol)
                position_list.append(pos)
        return ase.Atoms(symbol_list, positions=position_list)

    def cap_atoms(self, cap_description: str = '') -> 'ImperfectZeotype':
        """
        Cap all of the atoms in the ImpefectZeotype

        :param cap_description: A short description for the caps that will be added
        :return: A copy of self with the correct parameters added
        """
        # TODO: Find a way to map caps back to original
        self.update_nl()  # might not be needed
        self.parent_zeotype.update_nl()
        new_self = self
        o_cap_pos = self.build_O_atoms_cap_dict()
        if o_cap_pos:
            o_caps = self.build_cap_atoms(o_cap_pos)
            new_self = self.add_atoms(o_caps, 'o_caps', cap_description)

        new_self.update_nl()  # might not be needed
        new_self.parent_zeotype.update_nl()
        h_cap_pos = new_self.build_H_atoms_cap_dict()
        if h_cap_pos:
            h_caps = new_self.build_cap_atoms(h_cap_pos)
            new_self = new_self.add_atoms(h_caps, 'h_caps', cap_description)

        return new_self

    def remove_caps(self, cap_type: str = 'h_cap', cap_name: str = "cap") -> 'ImperfectZeotype':
        """
        Remove caps from an imperfect maze

        :param cap_type: The type of cap (h_cap, o_cap)
        :param cap_name: The name of the cap
        :return: A copy of self with the caps removed
        """
        assert cap_name in self.additions[cap_type], 'cap not in additions'
        indices_to_delete = self.index_mapper.get_overlap(self.name, cap_name)
        new_self = self.delete_atoms(indices_to_delete)
        new_self.additions[cap_type].remove(cap_name)
        return new_self

    def integrate_adsorbate(self, adsorbate: Atoms) -> Tuple['ImperfectZeotype', Adsorbate]:
        """
        Add an adsorbate into the imperfect maze

        :param adsorbate: Adsorbate object
        :param ads_name: name of the adsorbate
        :return: a copy of imperfect maze with the adsorabte added
        """
        ads_name = 'adsorbate'
        ads = Adsorbate(adsorbate)
        new_self = self.add_atoms(adsorbate, ads_name, short_description=ads.description)
        ads.name = new_self.additions[ads_name][-1]
        return new_self, ads

    def remove_adsorbate(self, adsorbate: Union[Adsorbate, str]) -> "ImperfectZeotype":
        """
        Removes an existing adsorbate from the ImperfectZeotype
        :param adsorbate: adsorbate object or str
        :return: new imperfect maze with item removed
        """
        ads_cat = 'adsorbate'
        if hasattr(adsorbate, "name"):
            ads_name = adsorbate.name
        else:
            assert isinstance(adsorbate, str), 'adsorbate must have .name attr or be str'
            ads_name = adsorbate
        assert adsorbate.name in self.additions[ads_cat], 'ads not in additions'
        indices_to_delete = self.index_mapper.get_overlap(self.name, ads_name)
        new_self = self.delete_atoms(indices_to_delete)
        new_self.additions[ads_cat].remove(ads_name)
        return new_self

    def integrate_other_zeotype(self, other: Zeotype):
        """
        Integrate another maze into the current maze

        :param other: the other maze to integrate
        :return: a new imperfect maze with the other maze integrated
        """
        new_self = self.__class__(self)
        atoms_to_add = ase.Atoms()
        for atom in other:
            self_index = new_self.index_mapper.get_index(other.name, new_self.name, atom.index)
            if self_index is not None:
                new_self.change_atom_properties(self_index, atom.index, other)
            else:
                atoms_to_add.extend(atom)

        new_self = new_self.add_atoms(atoms_to_add, atom_type='other_zeotype')

        return new_self

    def change_atom_properties(self, self_index: int, other_index: int, other: Atoms) -> None:
        """
        Change the atom properties

        :param self_index: index of self that needs to be changed
        :param other_index: index of other that holds the properties to change self
        :param other: other Atoms object that the other_indices corresponds to
        :return: None
        """
        self[self_index].symbol = other[other_index].symbol
        self[self_index].position = other[other_index].position
        self[self_index].tag = other[other_index].tag
        self[self_index].momentum = other[other_index].momentum
        self[self_index].mass = other[other_index].mass
        self[self_index].magmom = other[other_index].magmom
        self[self_index].charge = other[other_index].charge

    def build_H_atoms_cap_dict(self, bonds_needed: Optional[Dict[str, int]] = None):
        """
        Build a dictionary of hydrogen caps

        :param bonds_needed: Number of bonds needed for each atom
        :return: A list of positions of H caps
        """
        self.update_nl()
        if bonds_needed is None:
            bonds_needed = {'O': 2, 'Si': 4, 'Sn': 4, 'Al': 4, 'Ga': 4, 'B': 4}
        cap_atoms_dict: Dict[str, List[int]] = defaultdict(list)
        indices, count = self.count_elements()
        for o_index in indices['O']:
            if self.needs_cap(o_index, bonds_needed['O']):
                pos = self.get_H_pos(o_index)
                cap_atoms_dict['H'].append(pos)

        return dict(cap_atoms_dict)

    def build_O_atoms_cap_dict(self, bonds_needed: Optional[Dict[str, int]] = None):
        """
        Builds a dictionary of oxygen caps

        :param bonds_needed: Number of bonds needed for each atom
        :return: A list of oxygen cap positions
        """
        self.update_nl()
        if bonds_needed is None:
            bonds_needed = {'O': 2, 'Si': 4, 'Sn': 4, 'Al': 4, 'Ga': 4, 'B': 4}

        cap_atoms_dict: Dict[str, List[int]] = defaultdict(list)
        indices, count = self.count_elements()
        for si_index in indices['Si']:
            if self.needs_cap(si_index, bonds_needed['Si']):
                for i in range(bonds_needed['Si'] - len(self.neighbor_list.get_neighbors(si_index)[0])):
                    pos = self.get_oxygen_cap_pos(si_index)
                    cap_atoms_dict['O'].append(pos)

        return cap_atoms_dict

    def build_all_atoms_cap_dict(self, bonds_needed: Optional[Dict[str, int]] = None) -> Dict[str, np.array]:
        """
        Builds a dictionary for all atoms (not currently being used, but could be resurrected).

        :param bonds_needed: number of bonds needed for each atom type
        :return: A list of capping atoms to add
        """
        self.update_nl()
        if bonds_needed is None:
            bonds_needed = {'O': 2, 'Si': 4, 'Sn': 4, 'Al': 4, 'Ga': 4, 'B': 4}

        cap_atoms_dict: Dict[str, List[int]] = defaultdict(list)
        indices, count = self.count_elements()
        for si_index in indices['Si']:
            if self.needs_cap(si_index, bonds_needed['Si']):
                for i in range(bonds_needed['Si'] - len(self.neighbor_list.get_neighbors(si_index)[0])):
                    pos = self.get_oxygen_cap_pos(si_index)
                    cap_atoms_dict['O'].append(pos)

        for o_index in indices['O']:
            if self.needs_cap(o_index, bonds_needed['O']):
                pos = self.get_H_pos(o_index)
                cap_atoms_dict['H'].append(pos)

        return dict(cap_atoms_dict)

    def get_H_pos(self, atom_to_cap_self_i: int) -> np.array:
        """
        :param atom_to_cap_self_i: atom to cap index (self index)
        :return: hydrogen position
        """
        atom_to_cap_pi = self.index_mapper.get_index(self.name, self.parent_zeotype.name, atom_to_cap_self_i)
        if atom_to_cap_pi is None:
            site_pi = None
            print(f'atom_to_cap_self_i {atom_to_cap_self_i} does not map to parent maze')
        else:
            site_pi = self.find_missing_atom(atom_to_cap_pi, ['H', 'Si'])

        if site_pi is None:
            print(f'atom_to_cap_self_i {atom_to_cap_self_i} No matching atom in parent found using h_finder')
            return self.get_hydrogen_cap_pos_simple(atom_to_cap_self_i)

        direction = self.parent_zeotype.get_distance(atom_to_cap_pi, site_pi, mic=True, vector=True)
        hydrogen_pos = self.get_positions()[atom_to_cap_self_i] + direction / np.linalg.norm(direction)
        return hydrogen_pos

    def find_missing_atom(self, oxygen_atom_to_cap_pi, atom_symbol_list) -> int:
        """
        :param atom_symbol_list: The symbols to look for in the parent maze
        :param oxygen_atom_to_cap_pi:
        :return: parent atom Si index or H index
        """
        nl = self.parent_zeotype.neighbor_list.get_neighbors(oxygen_atom_to_cap_pi)[0]
        for atom_index in nl:
            if self.index_mapper.get_index(self.parent_zeotype.name, self.name, atom_index) is None:
                if self.parent_zeotype[atom_index].symbol in atom_symbol_list:
                    return atom_index

    def delete_atoms(self, indices_to_delete) -> 'ImperfectZeotype':
        """Delete atoms from imperfect maze by returning a copy with atoms deleted

        :param indices_to_delete: Indices of atoms in current maze to delete
        :return: a copy of self with atoms deleted
        """
        new_self_a = ase.Atoms(self)
        del new_self_a[indices_to_delete]
        new_self = self.__class__(new_self_a)
        self.set_attrs_source(new_self, self)
        old_to_new_map = self._get_old_to_new_map(self, new_self)
        self.index_mapper.register(self.name, new_self.name, old_to_new_map)
        return new_self

    @staticmethod
    def set_attrs_source(new_z: 'ImperfectZeotype', source: Zeotype) -> None:
        """
        Set the attributes of a new imperfect maze to that of its source

        :param new_z: Newly created zeolite without attributes set
        :param source: the source from which new_z was created
        :return: None
        """
        new_z.parent_zeotype = source.parent_zeotype
        new_z.index_mapper = source.index_mapper
        new_z.additions = copy.deepcopy(source.additions)
        new_z.name = new_z.index_mapper.get_unique_name(type(new_z).__name__)

    def _change_atoms(self, operation, *args, **kwargs) -> 'ImperfectZeotype':
        """
        Applies a custom function to the atoms and returns a new self

        :param operation:
        :param args:
        :param kwargs:
        :return:
        """
        new_self = self.__class__(self)
        operation(new_self, *args, **kwargs)
        old_to_new_map = self._get_old_to_new_map(self, new_self)
        self.index_mapper.add_name(new_self.name, self.name, old_to_new_map)
        return new_self

    def add_atoms(self, atoms_to_add: Atoms, atom_type: str, short_description: str = '') -> 'ImperfectZeotype':
        """
        Adds additional atoms to current imperfect maze

        :param atoms_to_add: The new atoms to add
        :param atom_type: A str describing the type of atoms that are being added
        :param short_description: A short description of the atoms that are being added (optional)
        :return: A new imperfect maze with the atoms added
        """
        # add atoms to indexer
        # register atoms_to_add with index mapper
        atom_indices = [atom.index for atom in atoms_to_add]
        if short_description:
            new_atom_name = self.index_mapper.get_unique_name(atom_type) + '_' + short_description
        else:
            new_atom_name = self.index_mapper.get_unique_name(atom_type)
        self.index_mapper.add_atoms(new_atom_name, atom_indices)

        # create a new self
        new_self_a = ase.Atoms(self)
        new_self_a.extend(atoms_to_add)
        new_self = self.__class__(new_self_a)
        self.set_attrs_source(new_self, self)

        # map new self to new atoms object and self
        self_to_new_self_map = self._get_old_to_new_map(self, new_self)
        atoms_to_new_self_map = self._get_old_to_new_map(atoms_to_add, new_self)

        # combine maps into a single map to main indexer
        main_to_new_self_map = {}
        # first map to self
        self_to_main_map = self.index_mapper.get_reverse_main_index(self.name)
        for self_i, new_self_i in self_to_new_self_map.items():
            main_index = self_to_main_map[self_i]
            main_to_new_self_map[main_index] = new_self_i
        # then map to atoms
        atoms_to_main_map = self.index_mapper.get_reverse_main_index(new_atom_name)
        for atoms_i, new_self_i in atoms_to_new_self_map.items():
            main_index = atoms_to_main_map[atoms_i]
            main_to_new_self_map[main_index] = new_self_i
        # finally register object
        self.index_mapper.register_with_main(new_self.name, main_to_new_self_map)
        # add new_atom_name to additions list
        new_self.additions[atom_type].append(new_atom_name)

        return new_self

    def remove_addition(self, addition_name, addition_type) -> 'ImperfectZeotype':
        """
        Removes an addition to the maze

        :param addition_name: name of the addition
        :param addition_type: the type of additon (h_cap, o_cap, ect.)
        :return: A new maze with the additional atoms remvoed
        """
        addition_to_self_map = self.index_mapper.get_name1_to_name2_map(addition_name, self.name)
        to_delete = list(addition_to_self_map.values())
        new_self = self.delete_atoms(to_delete)
        new_self.additions[addition_type].remove(addition_name)
        return new_self

    def create_silanol_defect(self, site_index) -> 'ImperfectZeotype':
        """
        Creates a silanol defect by deleting an atom and capping the resulting imperfect maze

        :param site_index:
        :return: An imperfect maze with a silanol defect
        """
        return self.delete_atoms([site_index]).cap_atoms()

    def needs_cap(self, atom_index: int, bonds_needed: int) -> bool:
        """
        Finds if an atom needs a cap

        :param atom_index: index of atom to check for cap
        :param bonds_needed: number of bonds an atom needs
        :return: boolean stating if an atom needs a cap
        """
        return len(self.neighbor_list.get_neighbors(atom_index)[0]) < bonds_needed

    def get_oxygen_cap_pos(self, atom_to_cap_self_i) -> np.array:

        """
        Find a position of an oxygen cap

        :param self_index: index of atom needing cap
        :return: A position array of the oxygen cap position
        """
        atom_to_cap_pi = self.index_mapper.get_index(self.name, self.parent_zeotype.name, atom_to_cap_self_i)
        site_pi = self.find_missing_atom(atom_to_cap_pi, ['O'])

        # parend_index = self.index_mapper.get_index(self.name, self.parent_zeotype.name, self_index)
        if site_pi is None:
            print("no matching parent oxygen found using old method")
            self_index = atom_to_cap_self_i
            self.update_nl()
            neighbor = self.neighbor_list.get_neighbors(self_index)[0][
                -1]  # first index in the list of neighbor indicies
            direction = self.get_positions()[self_index] - self.get_positions()[neighbor]  # vector from neighbor to Si
            oxygen_pos = self.get_positions()[self_index] + 1.6 * direction / np.linalg.norm(direction)
            return oxygen_pos
        else:
            return self.parent_zeotype.get_positions()[site_pi]

    def get_hydrogen_cap_pos_simple(self, index) -> np.array:
        """
        Finds the position of a hydrogen cap position

        :param index: index of hydrogen cap
        :return: the hydrogen position to add the cap too
        """
        neighbor = self.neighbor_list.get_neighbors(index)[0][0]  # first index in the list of neighbor indicies
        direction = self.get_positions()[index] - self.get_positions()[neighbor]  # vector from neighbor to oxygen
        hydrogen_pos = self.get_positions()[index] + direction / np.linalg.norm(direction)
        return hydrogen_pos


class OpenDefect(ImperfectZeotype):
    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, site_to_atom_indices=None, atom_indices_to_site=None,
                 additions=None):
        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms,
                         charges, scaled_positions, cell, pbc, celldisp, constraint,
                         calculator, info, velocities, site_to_atom_indices,
                         atom_indices_to_site, additions)

    @classmethod
    def build_from_indices(cls, parent_zeotype: Zeotype, indices_to_delete: Iterable[int]) -> 'OpenDefect':
        new_od = cls(parent_zeotype)
        new_od.set_attrs_source(new_od, parent_zeotype)

        old_to_new_map = Zeotype._get_old_to_new_map(parent_zeotype, new_od)
        new_od.index_mapper.add_name(new_od.name, parent_zeotype.name, old_to_new_map)
        new_od = new_od.delete_atoms(indices_to_delete)
        return new_od


class Cluster(ImperfectZeotype):  # TODO include dynamic inheritance and

    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, site_to_atom_indices=None, atom_indices_to_site=None,
                 additions=None):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms,
                         charges, scaled_positions, cell, pbc, celldisp, constraint,
                         calculator, info, velocities, site_to_atom_indices,
                         atom_indices_to_site, additions)

    @staticmethod
    def build_from_zeolite(parent_zeotype: Zeotype, index: int, max_cluster_size: int, max_neighbors: int,
                           cluster_indices: Iterable[int] = None) -> "Cluster":
        if cluster_indices is None:
            cluster_indices = Cluster.get_cluster_indices(parent_zeotype, index, max_cluster_size, max_neighbors)
        cluster = Cluster(parent_zeotype)
        to_delete = cluster.get_indices_compliment(cluster, cluster_indices)
        cluster = cluster.delete_atoms(to_delete)
        return cluster

    @staticmethod
    def get_oh_cluster_multi_t_sites(zeolite: Zeotype, t_sites: int) -> List[int]:
        """
        get an OH cluster with multiple T sites
        :param zeolite: The maze from which to extract the cluster
        :param t_sites: the central t site
        :return: A list of indices of the cluster
        """
        all_indces = set()
        for t_site in t_sites:
            all_indces.update(Cluster.get_oh_cluster_indices(zeolite, t_site))

        return list(all_indces)

    @staticmethod
    def get_oh_cluster_indices(zeolite: Zeotype, t_site: int) -> List[int]:
        """
        Create a cluster that only includes one central T site and then Oxygen
        and Hydrogen atoms. This is different than the other cluster selection
        methods that take in other
        :param zeolite: The zeolite from which the cluster indices will be drawn
        :param t_site: The index of the T site around which the cluster will be built
        :return: The indices of the new cluster
        """

        nl = NeighborList(natural_cutoffs(zeolite), self_interaction=False, bothways=True)
        nl.update(zeolite)

        all_indices = set([t_site])
        oxygen_indices = set()
        for index in nl.get_neighbors(t_site)[0]:
            if zeolite[index].symbol == "O":
                oxygen_indices.add(index)

        si_indices = set()
        for oxygen_index in oxygen_indices:
            for index in nl.get_neighbors(oxygen_index)[0]:
                if zeolite[index].symbol == "Si":
                    si_indices.add(index)

        for si_index in si_indices:
            for index in nl.get_neighbors(si_index)[0]:
                if zeolite[index].symbol == "O":
                    oxygen_indices.add(index)

        all_indices = all_indices.union(oxygen_indices).union(si_indices)

        return list(all_indices)

    @staticmethod
    def get_cluster_indices(zeolite, index: int, max_size: int, max_neighbors: int) -> List[int]:
        """
        get the indices of a cluster from a zeolite when specifying the
        center atom index and size of the cluster

        :param zeolite: the zeolite from which to build the cluster
        :param index: the centeral atom index
        :param max_size: the max number of atoms in the final cluster
        :param max_neighbors: the max number of neighbors from the starting cluster
        :return: a list of indices
        """
        nl = NeighborList(natural_cutoffs(zeolite), self_interaction=False, bothways=True)
        # instead of using zeolite use only the T sites
        # Sn Ti Hf, Si , Al, Zn T sites
        # Look at the
        # The remove functions should take all the T site elements or what T sites
        # what to remove should be different from the function that actually removes the T sites
        # User

        # if I give it the indices of 5 T sites. I remove 5 Si atoms I should create 5 * 4 O-H bonds
        #
        nl.update(zeolite)
        cluster_indices = set()
        new_cluster_indices = set([index])

        for _ in range(max_neighbors):
            current_cluster_indices = set()
            for cluster_index in new_cluster_indices:
                cluster_indices.add(cluster_index)
                if len(cluster_indices) >= max_size:
                    return list(cluster_indices)
                for new_index in nl.get_neighbors(cluster_index)[0]:
                    current_cluster_indices.add(new_index)
            new_cluster_indices = current_cluster_indices

        return list(cluster_indices)

    @staticmethod
    def get_cluster_indices_multi_T_site(zeolite, T_indices: Iterable[int], max_size: int, max_neighbors: int) -> List[
        int]:
        """
        get the indices of a cluster from a zeolite when specifying the
        center atom index and size of the cluster

        :param zeolite: the zeolite from which to build the cluster
        :param index: the centeral atom index
        :param max_size: the max number of atoms in the final cluster
        :param max_neighbors: the max number of neighbors from the starting cluster
        :return: a list of indices
        """
        nl = NeighborList(natural_cutoffs(zeolite), self_interaction=False, bothways=True)
        nl.update(zeolite)
        cluster_indices = set()
        new_cluster_indices = set(T_indices)

        for _ in range(max_neighbors):
            current_cluster_indices = set()
            for cluster_index in new_cluster_indices:
                cluster_indices.add(cluster_index)
                if len(cluster_indices) >= max_size:
                    return list(cluster_indices)
                for new_index in nl.get_neighbors(cluster_index)[0]:
                    if new_index not in T_indices:  # don't add T sites to current cluster indices
                        current_cluster_indices.add(new_index)
            new_cluster_indices = current_cluster_indices

        return list(cluster_indices)


# testing
if __name__ == '__main__':
    from ase.io import read

    b = read('/data/sam/sn_ti-periodic.traj')
    z = Zeotype(b)
    view(b)
    si_list = z.find_silanol_groups()
    print(si_list)
    sites_list = []
    for sil in si_list:
        # if 'T' in z.atom_indices_to_site(sil.Si_index):
        # sites_list.append(sil.Si_index)
        # else:
        for index in sil.Si_neighbor_list:
            # if 'T' in z.atom_indices_to_site(index):
            if 'Sn' == z[index].symbol:
                sites_list.append(index)
    print(list(set(sites_list)))
