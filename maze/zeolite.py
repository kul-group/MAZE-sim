import copy
from collections import defaultdict
from typing import List, Dict, Tuple, Iterable, Optional, Union, Callable
import warnings

import ase
import ase.data
import numpy as np
from ase import Atoms
from ase.neighborlist import natural_cutoffs, NeighborList
from maze.adsorbate import Adsorbate
from maze.perfect_zeolite import PerfectZeolite
from abc import ABC


class Zeolite(PerfectZeolite):
    """
    This Zeolite object inherits from PerfectZeolite. It represents a zeolite and consists of additional functionality
    for adding, removing and integrating other atoms.
    """

    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, site_to_atom_indices=None, atom_indices_to_site=None,
                 additions=None, ztype=None, cluster_maker=None):
        """
        Constructor for Zeolite object. The majority of the functionality comes from ase.Atoms object
        # TODO: write out docstring
        """

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms,
                         charges, scaled_positions, cell, pbc, celldisp, constraint,
                         calculator, info, velocities, site_to_atom_indices, atom_indices_to_site,
                         additions, _is_zeotype=False, ztype=ztype)

        # handle cluster maker
        if isinstance(symbols, self.__class__) and cluster_maker is None:  # build from IZ
            self.cluster_maker = copy.deepcopy(symbols.cluster_maker)
        elif cluster_maker is None:  # cluster maker argument passed and not built from IZ
            self.cluster_maker = DefaultClusterMaker()
        else:  # cluster maker argument passed
            self.cluster_maker = cluster_maker

    @staticmethod
    def build_cap_atoms(cap_atoms_dict: Dict[str, np.array]) -> Atoms:
        symbol_list = []
        position_list = []
        for symbol, pos_list in cap_atoms_dict.items():
            for pos in pos_list:
                symbol_list.append(symbol)
                position_list.append(pos)
        return ase.Atoms(symbol_list, positions=position_list)

    def cap_atoms(self, cap_description: str = '') -> 'Zeolite':
        """
        Cap all of the atoms in the zeolite

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

    def remove_caps(self, cap_type: str = 'cap', cap_name: Optional[str] = None) -> 'Zeolite':
        """
        Remove caps from an imperfect MAZE-sim, if no arguments are provided all
        caps are removed from the Zeolite

        :param cap_type: The type of cap (h_cap, o_cap or cap for both) partial matching will work
        :param cap_name: The name of the cap
        :return: A copy of self with the caps removed
        """
        cap_names_to_delete = []
        cap_names_to_delete_dict = defaultdict(list)
        for addition_type in self.additions.keys():
            if cap_type not in addition_type:
                continue
            if cap_name is None:  # select all names
                cap_names_to_delete.extend(self.additions[addition_type])
                cap_names_to_delete_dict[addition_type].extend(self.additions[addition_type])
            else:
                matching_names = [name for name in self.additions[addition_type] if cap_name in name]
                cap_names_to_delete.extend(matching_names)
                cap_names_to_delete_dict[addition_type].extend(matching_names)

        assert cap_names_to_delete, 'no matching caps found'

        indices_to_delete = []
        for cap_name in cap_names_to_delete:
            indices_to_delete.extend(self.index_mapper.get_overlap(self.name, cap_name))

        new_self = self.delete_atoms(indices_to_delete)
        for cap_type in cap_names_to_delete_dict.keys():
            for cap_name in cap_names_to_delete_dict[cap_type]:
                new_self.additions[cap_type].remove(cap_name)

        return new_self

    def integrate_adsorbate(self, adsorbate: Atoms) -> Tuple['Zeolite', Adsorbate]:
        """
        Add an adsorbate into the imperfect MAZE-sim

        :param adsorbate: Adsorbate object
        :param ads_name: name of the adsorbate
        :return: a copy of imperfect MAZE-sim with the adsorabte added
        """
        ads_name = 'adsorbate'
        ads = Adsorbate(adsorbate)
        new_self = self.add_atoms(adsorbate, ads_name, short_description=ads.description)
        ads.name = new_self.additions[ads_name][-1]
        return new_self, ads

    def remove_adsorbate(self, adsorbate: Union[Adsorbate, str]) -> "Zeolite":
        """
        Removes an existing adsorbate from the Zeolite
        :param adsorbate: adsorbate object or str
        :return: new imperfect MAZE-sim with item removed
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

    def integrate(self, other: PerfectZeolite) -> "Zeolite":
        """
        Integrate another MAZE-sim into the current MAZE-sim

        :param other: the other MAZE-sim to integrate
        :return: a new imperfect MAZE-sim with the other MAZE-sim integrated
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
        cap_atoms_dict: Dict[str, List[np.array]] = defaultdict(list)
        indices, count = self.count_elements()
        for o_index in indices['O']:
            if self.needs_cap(o_index, bonds_needed['O']):
                pos = self.get_H_pos(o_index)
                cap_atoms_dict['H'].append(np.array(pos))

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

        cap_atoms_dict: Dict[str, List[np.array]] = defaultdict(list)
        indices, count = self.count_elements()
        distict_positions = set()
        for si_index in indices['Si']:
            if self.needs_cap(si_index, bonds_needed['Si']):
                for i in range(bonds_needed['Si'] - len(self.neighbor_list.get_neighbors(si_index)[0])):
                    pos = self.get_oxygen_cap_pos(si_index)
                    distict_positions.add(tuple(pos))
        for pos in distict_positions:
            cap_atoms_dict['O'].append(np.array(pos))

        return dict(cap_atoms_dict)

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
            print(f'atom_to_cap_self_i {atom_to_cap_self_i} does not map to parent MAZE-sim')
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
        :param atom_symbol_list: The symbols to look for in the parent MAZE-sim
        :param oxygen_atom_to_cap_pi:
        :return: parent atom Si index or H index
        """
        nl = self.parent_zeotype.neighbor_list.get_neighbors(oxygen_atom_to_cap_pi)[0]
        for atom_index in nl:
            if self.index_mapper.get_index(self.parent_zeotype.name, self.name, atom_index) is None:
                if self.parent_zeotype[atom_index].symbol in atom_symbol_list:
                    return atom_index

    def delete_atoms(self, indices_to_delete) -> 'Zeolite':
        """Delete atoms from imperfect MAZE-sim by returning a copy with atoms deleted

        :param indices_to_delete: Indices of atoms in current MAZE-sim to delete
        :return: a copy of self with atoms deleted
        """
        # lambda a: del a[indices_to_delete]

        new_self_a = copy.deepcopy(ase.Atoms(self))
        del new_self_a[indices_to_delete]
        new_self = self.__class__(new_self_a, ztype=self.ztype)
        self.set_attrs_source(new_self, self)
        old_to_new_map = self._get_old_to_new_map(self, new_self)
        self.index_mapper.register(self.name, new_self.name, old_to_new_map)
        return new_self

    @staticmethod
    def set_attrs_source(new_z: 'Zeolite', source: "Zeolite") -> None:
        """
        Set the attributes of a new imperfect MAZE-sim to that of its source

        :param new_z: Newly created zeolite without attributes set
        :param source: the source from which new_z was created
        :return: None
        """
        new_z.parent_zeotype = source.parent_zeotype
        new_z.index_mapper = source.index_mapper
        new_z.additions = copy.deepcopy(source.additions)
        new_z.cluster_maker = copy.deepcopy(source.cluster_maker)
        # copy over potential energies
        # copy over forces
        # copy over other key atoms attributes

    def translate_self(self, displacement) -> "Zeolite":
        """
        Returns a copy of self with applied translation
        :param displacement: atomic positions
        :type displacement: Iterable
        :return: new self translatedß
        :rtype: Zeolite
        """
        new_self = self.apply(lambda a: a.translate(displacement))
        return new_self

    def wrap_self(self, **wrap_kw) -> "Zeolite":
        new_self = self.apply(lambda a: a.wrap(**wrap_kw))
        return new_self

    @staticmethod
    def _count_equal(list1: Iterable, list2: Iterable) -> bool:
        count_dict_1 = defaultdict(lambda: 0)
        count_dict_2 = defaultdict(lambda: 0)
        for el1, el2 in zip(list1, list2):
            count_dict_1[el1] += 1
            count_dict_2[el2] += 1

        if len(count_dict_1) != len(count_dict_2):
            return False

        for key in count_dict_1:
            if count_dict_1[key] != count_dict_2[key]:
                return False

        return True

    def apply(self, operation: Callable, *args, **kwargs) -> 'Zeolite':
        """
        Applies a custom function to the atoms and returns a new self.
        Custom operation cannot change tags or else errors will occur!

        :param operation: applies a custom function to the atoms
        :return: zeolite with applied modifications
        """
        new_self_a = copy.deepcopy(ase.Atoms(self.retag_self()))
        pre_operation_tags = new_self_a.get_tags()
        operation(new_self_a, *args, **kwargs)
        post_operation_tags = new_self_a.get_tags()
        assert self._count_equal(pre_operation_tags, post_operation_tags), 'tags cannot be changed by operation'
        new_self = self.__class__(new_self_a)
        new_self.register_with_parent(self.parent_zeotype, self.build_additions_map())
        return new_self

    def get_type(self, index: int) -> 'str':
        """
        Get the type of atom at a certain index
        :param index: the atom index to get the type of
        :type index: int
        :return: the name of the type of index
        :rtype: str
        """
        assert self.parent_zeotype._atom_indices_to_site is not None, 'Parent Zeotype site mapping missing'
        return self.parent_zeotype._atom_indices_to_site[
            self.index_mapper.get_index(self.name, self.parent_zeotype.name, index)]

    def get_site_types(self) -> List[str]:
        """
        Get a list of all of the types in the parent zeolite
        :return: List of the names of all types
        :rtype: List[str]
        """
        assert self.parent_zeotype._site_to_atom_indices is not None, 'Parent Zeotype site mapping missing'
        assert self.parent_zeotype._atom_indices_to_site is not None, 'Parent Zeotype site mapping missing'
        return list(self.parent_zeotype._site_to_atom_indices.keys())

    def find_type(self, atom_type_name: str) -> List[int]:
        """
        Find the type of the atom in the parent zeolite list

        :param atom_type_name: the name of a certain site
        :type atom_type_name: str
        :return: List of indices that match site description
        :rtype: List[int]
        """
        assert self.parent_zeotype._site_to_atom_indices is not None, 'Parent Zeotype site mapping missing'
        assert atom_type_name in self.parent_zeotype._site_to_atom_indices, 'atom type name not found'
        parent_list = self.parent_zeotype._site_to_atom_indices[atom_type_name]
        self_list = []
        for p_index in parent_list:
            self_index = self.index_mapper.get_index(self.parent_zeotype.name, self.name, p_index)
            if self_index is not None:
                self_list.append(self_index)

        return self_list

    @staticmethod
    def _check_unique_positions(position_list):
        """
        A helper function to check that the positions of the atoms are unique. An error is thrown
        if there are overlapping atoms. This slows down the code, because it is O(N^2).

        :param position_list (np.array[np.array]):  an array of positions
        :return: None

        """
        for i in range(0, len(position_list)):
            for j in range(i, len(position_list)):
                if i == j:
                    continue
                for x1, y1 in zip(position_list[i], position_list[j]):
                    if x1 != y1:
                        break
                else:  # if break statement not called two vectors are equal (which means overlapping atoms)
                    raise ValueError('An Atoms object with Overlapping atoms was passed to the add_atoms method. '
                                     'Overlapping atoms cause trouble and thus are not allowed. ')

    def add_atoms(self, atoms_to_add: Atoms, atom_type: str, short_description: str = '') -> 'Zeolite':
        """
        Adds additional atoms to current imperfect MAZE-sim

        :param atoms_to_add: The new atoms to add
        :param atom_type: A str describing the type of atoms that are being added
        :param short_description: A short description of the atoms that are being added (optional)
        :return: A new imperfect MAZE-sim with the atoms added
        """
        # add atoms to indexer
        # register atoms_to_add with index mapper
        self._check_unique_positions(atoms_to_add.get_positions())

        atom_indices = [atom.index for atom in atoms_to_add]
        # TODO: Find a way to map new atoms back to parent zeolite if that mapping exists
        if short_description:
            new_atom_name = self.index_mapper.get_unique_name(atom_type) + '_' + short_description
        else:
            new_atom_name = self.index_mapper.get_unique_name(atom_type)
        self.index_mapper.add_atoms(new_atom_name, atom_indices)

        # create a new self
        new_self_a = copy.deepcopy(ase.Atoms(self))
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

    def remove_addition(self, addition_name, addition_type) -> 'Zeolite':
        """
        Removes an addition to the MAZE-sim

        :param addition_name: name of the addition
        :param addition_type: the type of additon (h_cap, o_cap, ect.)
        :return: A new MAZE-sim with the additional atoms remvoed
        """
        addition_to_self_map = self.index_mapper.get_name1_to_name2_map(addition_name, self.name)
        to_delete = list(addition_to_self_map.values())
        new_self = self.delete_atoms(to_delete)
        new_self.additions[addition_type].remove(addition_name)
        return new_self

    def create_silanol_defect(self, site_index) -> 'Zeolite':
        """
        Creates a silanol defect by deleting an atom and capping the resulting imperfect MAZE-sim

        :param site_index:
        :return: An imperfect MAZE-sim with a silanol defect
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
            nl = self.neighbor_list.get_neighbors(self_index)[0]
            if len(nl) != 0:
                neighbor = nl[-1]  # first index in the list of neighbor indicies
                direction = self.get_positions()[self_index] - self.get_positions()[
                    neighbor]  # vector from neighbor to Si
                oxygen_pos = self.get_positions()[self_index] + 1.6 * direction / np.linalg.norm(direction)
            else:
                direction = np.array([0, 0, 1])
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

    def get_cluster(self, start_index: int, cluster_indices=None, **kwargs) -> Tuple["Zeolite", "Zeolite"]:
        """
        Generates a Cluster of atoms around the specified index. The number of atoms in the cluster
        is given by the size parameter.

        :param start_index: center index of cluster
        :type start_index: int
        :param cluster_indices: Indices of cluster to make
        :param max_size: max size of cluster
        :param index: index of the central atom in the cluster
        :param max_neighbors: number of neighbors from the host atom in the final cluster
        :return: index of the cluster in the MAZE-sim cluster array
        """
        if cluster_indices is None:
            self_2 = self
            cluster_indices = self.cluster_maker.get_cluster_indices(self_2, start_index, **kwargs)

        cluster = self.cluster_maker.get_cluster(self, cluster_indices)
        open_defect = self.cluster_maker.get_open_defect(self, cluster_indices)
        return cluster, open_defect

    def add_additions_map_to_self(self, additions_map: Dict) -> None:
        """
        Add an additions map to the current Zeolite. This function has side effects and modifies the current Zeolite.
        This is useful in building a Zeolite from a saved file and reloading the additions map into memory.
        :param additions_map: an additions map (format is specified in the additions map function)
        :type additions_map: Dict[Dict]
        :return: None
        :rtype: None
        """
        for key in additions_map:
            self.additions[key].extend(list(additions_map[key].keys()))

    def register_with_parent(self, parent_zeolite: PerfectZeolite, additions_map: Optional[Dict] = None) -> None:
        """
        Register a tagged zeolite with the parent zeolite/zeolite to register the zeolite with. For this to work
        the tags of the current zeolite must be equal to the parent indices.
        :param parent_zeolite: the parent zeolite to register the zeolite with.
        :type parent_zeolite: PerfectZeolite
        :param additions_map: Additions map containing addition names and their parent indices
        :type additions_map: Optional[Dict]
        :return: None
        :rtype: None
        """
        # test that all tags are unique
        if len(self.get_tags()) != len(set(self.get_tags())):
            raise RuntimeError('Tags must be unique to register self with parent_zeolite')

        main_to_new_map = {}
        for atom in self:
            main_to_new_map[atom.tag] = atom.index

        parent_zeolite.index_mapper.register_with_main(self.name, main_to_new_map)
        self.index_mapper = parent_zeolite.index_mapper
        self.parent_zeotype = parent_zeolite

        if additions_map is not None:  # register additions
            self.add_additions_map_to_self(additions_map)
            for key in additions_map:
                addition_main_to_new_map = {}
                for name, main_indices in additions_map[key].items():
                    addition_index = 0
                    for atom in self:
                        if atom.tag in main_indices:
                            addition_main_to_new_map[atom.tag] = addition_index
                            addition_index += 1
                    self.index_mapper.register_with_main(name, addition_main_to_new_map)

    def retag_self(self) -> None:
        """
        Add tags to the Zeolite that correspond to their main index in the index mapper
        :return: None
        :rtype: None
        """
        assert self.index_mapper is not None, 'cannot retag when index mapper is None'
        self_to_main_index_map = self.index_mapper.get_reverse_main_index(self.name)
        new_self = self.__class__(self)
        for atom in new_self:
            atom.tag = self_to_main_index_map[atom.index]

        return new_self


class ClusterMaker(ABC):
    """
    This is an abstract class that selects indices from a Zeolite to return a cluster.
    Multiple implementations of a cluster maker are possible.
    """

    def get_cluster_indices(self, zeolite: PerfectZeolite, start_site: int, **kwargs) -> List[int]:
        """
        Get the cluster indices of a zeo
        :param zeolite: zeolite to select the cluster from
        :type zeolite: Zeolite
        :param start_site: the start site (usually a T-site) from which to select the indices
        :type start_site: int
        :param kwargs: additional keyword arguments needed for other functions
        :type kwargs: Dict
        :return: list of indices
        :rtype: List[int]
        """
        raise NotImplementedError

    @staticmethod
    def get_cluster(zeolite: Zeolite, indices: Iterable[int], name="Cluster") -> Zeolite:
        """
        Get a cluster from a zeolite
        :param zeolite: Zeolite from which to extract the cluster
        :type zeolite: Zeolite
        :param indices: the indices of the zeolite to construct the cluster
        :type indices: Iterable[int]
        :param name: the ztype of the created cluster
        :type name: string
        :return: The created cluster of type Zeolite
        :rtype: Zeolite
        """
        cluster = Zeolite(zeolite, ztype=name)
        to_delete = cluster.get_indices_compliment(cluster, indices)
        cluster = cluster.delete_atoms(to_delete)
        return cluster

    @staticmethod
    def get_open_defect(zeolite: PerfectZeolite, indices: Iterable[int], name="Open Defect") -> Zeolite:
        """
        Get an opendefect object from a zeolite
        :param zeolite: zeolite from which to get the open defect
        :type zeolite: Zeolite
        :param indices: the indices of the cluster (these indices are going to be deleted)
        :type indices: Iterable[int]
        :param name: the ztype of the open defect that will be created
        :type name: str
        :return: The created opendefect Zeolite
        :rtype: Zeolite
        """
        new_od = Zeolite(zeolite, ztype=name)
        new_od = new_od.delete_atoms(indices)
        return new_od


class DefaultClusterMaker(ClusterMaker):
    """
    This is an implementation of a Cluster Maker that is used if no cluster maker object is specified
    """

    def __init__(self):
        pass

    def get_cluster_indices(self, zeolite: PerfectZeolite, start_site: int, **kwargs):
        return self.get_oh_cluster_indices(zeolite, start_site)

    @classmethod
    def get_oh_cluster_multi_t_sites(cls, zeolite: PerfectZeolite, t_sites: Iterable[int]) -> List[int]:
        """
        get an OH cluster with multiple T sites
        :param zeolite: The MAZE-sim from which to extract the cluster
        :param t_sites: the central t site
        :return: A list of indices of the cluster
        """
        all_indices = set()
        for t_site in t_sites:
            all_indices.update(cls.get_oh_cluster_indices(zeolite, t_site))
        return list(all_indices)

    @staticmethod
    def get_oh_cluster_indices(zeolite: PerfectZeolite, t_site: int) -> List[int]:
        """
        Create a cluster that only includes one central T site and then Oxygen
        and Hydrogen atoms. This is different than the other cluster selection
        methods that take in other
        :param zeolite: The zeolite from which the cluster indices will be drawn
        :param t_site: The index of the T site around which the cluster will be built
        :return: The indices of the new cluster
        """
        if zeolite.parent_zeotype._site_to_atom_indices is not None and 'T' not in zeolite.get_site_type(t_site):
            # short circuit to prevent get_site_type from throwing error
            warnings.warn(f'Cluster indices generator requires a T site for the initial '
                          f'start site. Index {t_site} might not be a T site')
        nl = NeighborList(natural_cutoffs(zeolite), self_interaction=False, bothways=True)
        nl.update(zeolite)

        all_indices = {t_site}
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
    def get_cluster_indices_all_atoms(zeolite, index: int, max_size: int, max_neighbors: int) -> List[int]:
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
        new_cluster_indices = {index}

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

        :param T_indices: Indices of the T site
        :type T_indices: Iterable[int]
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
