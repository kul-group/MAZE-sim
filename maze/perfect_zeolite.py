import copy
import os
from collections import defaultdict
from typing import List, Dict, Tuple, Iterable, Optional
import uuid

import ase
import ase.data
from ase import Atoms
from ase.io import cif
from ase.neighborlist import natural_cutoffs, NeighborList
from packaging.version import parse as parse_version
import pkg_resources
from maze.index_mapper import IndexMapper
from maze.cif_download import download_cif

class PerfectZeolite(Atoms):
    """
    A class that inherits from ase.Atoms, which represents an 'perfect' zeolite. If a zeolite is built from a cif
    file from the zeolite structure database, then the atoms are tagged according to their unique site and the
    dictionaries site_to_atom_indices and atom_indices_to_site are filled. If not these dictionaries are set to none.
    This class contains a bunch of static methods for identifying different types of atoms in the zeolite as well as
    as methods to build an imperfect MAZE-sim from a parent Zeotype class. Imperfect zeotypes have additional
    functionality and are dependent on a parent MAZE-sim class.
    """
    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, site_to_atom_indices=None, atom_indices_to_site=None,
                 additions=None, _is_zeotype=True, ztype=None):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms, charges, scaled_positions,
                         cell, pbc, celldisp, constraint, calculator, info, velocities)

        # To Be compatible with ASE's atoms object, this code has to have the functionality to build from
        # an atom's object, a Zeotype object, or a sub class of a Zeotype object
        # The following if statements take care of this functionality depending on
        # if the object being built is a Zeotype or if symbols is a Zeotype (there are four unique paths)
        if isinstance(symbols, PerfectZeolite):  # if symbols is Zeotype or Zeotype subclass
            self.additions = copy.deepcopy(symbols.additions)
            if _is_zeotype:  # if the object being built is a Zeotype
                self._site_to_atom_indices = symbols._site_to_atom_indices
                self._atom_indices_to_site = symbols._atom_indices_to_site
                self.ztype = 'parent'
                self.name = 'parent'  # must be parent to agree with index mapper
                self.index_mapper = IndexMapper(self.get_indices(self))
                self.parent_zeotype = self

            else:  # if the object being built is a subtype of Zeotype
                self.parent_zeotype = symbols.parent_zeotype
                self._site_to_atom_indices = None
                self._atom_indices_to_site = None
                self.index_mapper = symbols.index_mapper
                if ztype is None:
                    self.ztype = symbols.ztype if symbols.ztype != 'parent' and symbols.ztype else type(self).__name__
                else:
                    self.ztype = ztype
                self.name = self.index_mapper.get_unique_name(self.ztype)
                # add name and register with index mapper
                self.index_mapper.register(symbols.name, self.name, self._get_old_to_new_map(symbols, self))

        else:  # if symbols is not a Zeotype or Zeotype child class
            if _is_zeotype:  # if Zeotype object is being built  # TODO: get rid of _is_zeotype attribute
                self.ztype = 'parent'
                self.name = 'parent'  # must be parent for code to work properly
                self.index_mapper = IndexMapper(self.get_indices(self))
                self.parent_zeotype = self
            else:   # building a non-parent zeolite
                # make a parent zeolite of self
                parent = PerfectZeolite(symbols)
                self.parent_zeotype = parent
                self.index_mapper = parent.index_mapper

                if ztype is not None:
                    self.ztype = ztype
                else:
                    self.ztype = type(self).__name__  # use type for ztype by default

                self.name = self.index_mapper.get_unique_name(self.ztype)  # use name
                self.index_mapper.register(parent.name, self.name, self._get_old_to_new_map(parent, self))

            self._site_to_atom_indices = site_to_atom_indices
            self._atom_indices_to_site = atom_indices_to_site
            self.additions = copy.deepcopy(additions) if additions else defaultdict(list)

        self.unique_id = str(uuid.uuid4())
        self.update_nl()

    @property
    def site_to_atom_indices(self) -> Dict:
        my_site_to_atom_indices = {}
        for site, indices in self.parent_zeotype._site_to_atom_indices.items():
            my_indices = []
            for index in indices:
                index = self.index_mapper.main_index[index][self.name]
                if index is not None:
                    my_indices.append(index)
            my_site_to_atom_indices[site] = my_indices

        return my_site_to_atom_indices

    @property
    def atom_indices_to_sites(self) -> Dict:
        my_atom_indices_to_sites = {}
        for index, site in self.parent_zeotype._atom_indices_to_site.items():
            my_index = self.index_mapper.main_index[index][self.name]
            if my_index is not None:
                my_atom_indices_to_sites[my_index] = site

        return my_atom_indices_to_sites


    @staticmethod
    def _get_available_symbols(atom_list: List[str]) -> List[str]:
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
    def build_from_cif_with_labels(cls, filepath: str, **kwargs) -> "PerfectZeolite":
        """
        Takes a filepath/fileobject of a cif file and returns a Zeotype or child class with T-sites
        labeled as specified in the cif file. The dictionaries for T site labels are also filled in
        these map from the T-site to a list of atom indices. The atom indices to t sites maps atom
        indices to T-sites.

        :param filepath: Filepath of the cif file
        :return: Zeotype with labeled sites
        """
        if parse_version(pkg_resources.get_distribution('ase').version) < parse_version('3.21.0'):
            # if using an older versino of ase
            cif_reader = cls._read_cif_note_sites
        else:
            # if using a newer version of ase
            cif_reader = cls._read_cif_note_siteJan2021Update

        atoms, site_to_atom_indices, atom_indices_to_site = cif_reader(filepath, **kwargs)
        zeotype = cls(atoms)
        zeotype._site_to_atom_indices = site_to_atom_indices
        zeotype._atom_indices_to_site = atom_indices_to_site
        return cls(zeotype)


    @staticmethod
    def _read_cif_note_siteJan2021Update(fileobj: str, store_tags=False, primitive_cell=False,
                                         subtrans_included=True, fractional_occupancies=True,
                                         reader='ase') -> Tuple[Atoms, Dict[str, int], Dict[int, str]]:

        """
        The helper function used by build_from_cif_with_labels when using an ASE version
        3.21.0 or higher. This loads a CIF file using ase.io.cif.parse_cif and then
        finds the T-site information. After finding the T-site information it then replaces
        each atom in the T-site with another  atom type not found in the cif file. The
        crystal is then generated and then the resulting atoms object atoms are replaced
        by the original atoms. The mapping between the T-site names and the atom objects
        in the final atoms object are recored in two dictionaries.

        A note about capability:
        ASE Version 3.21.0 released Jan 18, 2021 refactored the cif reading functions on
        on which the older function relies. The major change (from this function's perspective)
        is that ase.io.cif.parse_cif now returns a generator rather than a list, which contains
        a CIFBlock object.


        :param fileobj: CIF file location
        :param store_tags: store the tags in resulting atoms object
        :param primitive_cell: An option for the reader
        :param subtrans_included: An option for the reader
        :param fractional_occupancies: an option for the final crystal
        :param reader: the reader used (ase works others have not been tested)
        :return: atoms, site_to_atom_indices, atom_indices_to_site

        """

        blocks = ase.io.cif.parse_cif(fileobj, reader)  # get blocks generator
        cif = list(blocks)[0]  # get cif object

        # replace elements with replacement symbol
        element_to_T_site = {}
        sym_to_original_element = {}
        possible_syms = PerfectZeolite._get_available_symbols(
            cif._tags["_atom_site_type_symbol"])  # only get symbols not in CIF file
        for i in range(len(cif._tags["_atom_site_label"])):  # label all atoms
            sym = possible_syms.pop()
            sym_to_original_element[sym] = cif._tags["_atom_site_type_symbol"][i]
            cif._tags["_atom_site_type_symbol"][i] = sym
            element_to_T_site[sym] = cif._tags["_atom_site_label"][i]

        site_to_atom_indices = defaultdict(list)
        atom_indices_to_site = {}
        atoms = cif.get_atoms(store_tags, primitive_cell, subtrans_included, fractional_occupancies)

        for i in range(len(atoms)):  # replace substitute atoms with original Si and get indices of replaced atoms
            if atoms[i].symbol in element_to_T_site.keys():
                sym = atoms[i].symbol
                key = element_to_T_site[sym]
                site_to_atom_indices[key].append(i)
                atom_indices_to_site[i] = key
                atoms[i].tag = ase.data.atomic_numbers[sym]
                atoms[i].symbol = sym_to_original_element[sym]

        return atoms, dict(site_to_atom_indices), atom_indices_to_site

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
        possible_syms = PerfectZeolite._get_available_symbols(
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

    def update_nl(self, mult: float = 1) -> None:
        """
        Builds and updates neighborlist
        :param mult: The mult (multiply) parameter for natural cutoffs (Default 1)
        :return: None
        """
        self.neighbor_list = NeighborList(natural_cutoffs(self, mult=mult), bothways=True, self_interaction=False)
        self.neighbor_list.update(self)

    def get_hetero_atoms(self, hetero_atoms_list: Optional[List[str]] = None) -> List[int]:
        """
        :return: Returns a list of all of the hetero-atoms in the MAZE-sim
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
        for atom in self:  # iterate through atom objects in MAZE-sim
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
        :return: a dictionary where the key is the element symbol and the value is the number in the MAZE-sim
        """
        indices: Dict['str', List['int']] = defaultdict(list)  # indices of the elements grouped by type
        count: Dict['str', int] = defaultdict(lambda: 0)  # number of elements of each type
        for atom in self:
            element = atom.symbol
            indices[element].append(atom.index)
            count[element] += 1
        return dict(indices), dict(count)

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

    @staticmethod
    def get_indices_compliment(zeotype: 'PerfectZeolite', indices: Iterable[int]) -> List[int]:
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

    def extend(self, other):
        """
        This extends the current Zeotype with additional atoms
        :param other: atoms like object to extend with
        :type other: Atoms
        :return: None
        :rtype: None
        """
        self_length = len(self)
        try:
            other_indices = [self_length + i for i in range(0, len(other))]
        except TypeError:
            other_indices = [self_length + 1]


        self.index_mapper.extend(self.name, other_indices)
        super().extend(other)

    def pop(self, index: int = -1):
        """
        This removes
        :param index: index to pop
        :type index: int
        :return: Atom
        :rtype: Atom
        """

        self.index_mapper.pop(index)
        return super().pop(index)

    def get_site_type(self, index: int) -> str:
        """
        Get the identity of a site

        :param index: Index of site in self
        :return: Label for the site (comes from CIF) file
        """
        assert self.parent_zeotype._atom_indices_to_site is not None, 'atom_indices_to_site is None, cannot get label'
        pz_index = self.index_mapper.get_index(self.name, self.parent_zeotype.name, index)
        return self.parent_zeotype._atom_indices_to_site[pz_index]

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

    def __copy__(self):
        return self.__class__(self)

    def __del__(self) -> None:
        if hasattr(self, 'index_mapper') and self.index_mapper is not None:
            self.index_mapper.delete_name(self.name)

    @classmethod
    def make(cls, iza_code: str, data_dir='data'):
        """
        Builds an Zeolite from iza code
        :param iza_zeolite_code: zeolite iza code
        :type iza_zeolite_code: str
        :return: An imperfect zeolite class or subclass
        :rtype: cls
        """
        iza_code.capitalize()
        cif_path = os.path.join('data', iza_code + '.cif')
        if not os.path.exists(cif_path):
            download_cif(iza_code, data_dir)
        parent = PerfectZeolite.build_from_cif_with_labels(cif_path)
        return cls(parent)

    def build_additions_map(self):
        """
        Build a serializable additions map that can be used to rebuild the Zeolite from file
        :return: additions map
        :rtype: Dict
        """
        additions_map = defaultdict(dict)
        for category, names in self.additions.items():
            for name in names:
                additions_map[category][name] = list(self.index_mapper.get_reverse_main_index(name).values())
        return dict(additions_map)
