from typing import List, Dict, Tuple, Iterable
from ase import Atoms, Atom
from ase.neighborlist import natural_cutoffs, NeighborList
from collections import defaultdict
import copy
import numpy as np
from ase.io import cif
import ase
import ase.data
from ase.visualize import view
from source.index_mapper import IndexMapper


def get_available_symbols(atom_list):
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


class Silanol():
    def __init__(self, parent_zeotype, Si_index, O_index, H_index, Si_neighbor_list):
        self.parent_zeotype = parent_zeotype
        self.Si_index = Si_index
        self.O_index = O_index
        self.H_index = H_index
        self.Si_neighbor_list = Si_neighbor_list

    def __str__(self):
        return f'Si: {self.Si_index} O: {self.O_index} H: {self.H_index}'

    def __repr__(self):
        return self.__str__()


class Zeotype(Atoms):
    """
    This is a Zeotype class that inherits from Atoms. It represents a Zeolite.
    """

    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, silent: bool = False, zeolite_type: str = '',
                 t_site_to_atom_indices=None, atom_indices_to_t_site=None, build_index_mapper=True):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms, charges, scaled_positions,
                         cell, pbc, celldisp, constraint, calculator, info, velocities)

        if type(symbols) == Zeotype:  # copy over some, but not all, attributes of copying zeotype
            self.zeolite_type = symbols.zeolite_type
            self.silent = silent
            self.site_to_atom_indices = t_site_to_atom_indices
            self.atom_indices_to_site = atom_indices_to_t_site

        self.zeolite_type = zeolite_type
        self.sites: List[str] = []
        self.clusters: List[Cluster] = []
        self.adsorbates = []
        self.silent = silent
        self.neighbor_list = NeighborList(natural_cutoffs(self), bothways=True, self_interaction=False)
        self.neighbor_list.update(self)
        self.site_to_atom_indices = t_site_to_atom_indices
        self.atom_indices_to_site = atom_indices_to_t_site

        # things for index mapping
        self.index_mapper = None
        self.name = 'pristine'
        self.parent_zeotype = self
        if build_index_mapper:
            self_indices = [a.index for a in self]
            self.index_mapper = IndexMapper(self.name, self_indices)

    @classmethod
    def build_from_cif_with_labels(cls, filepath) -> "Zeotype":
        """
        Takes a filepath/fileobject of a cif file and returns a Zeotype or child class with T-sites
        labeled as specified in the cif file. The dictionaries for T site labels are also filled in
        these map from the T-site to a list of atom indices. The atom indices to t sites maps atom
        indices to T-sites.
        :param filepath:
        :return:
        """
        atoms, site_to_atom_indices, atom_indices_to_site = cls.read_cif_note_sites(filepath)
        zeotype = cls(atoms)
        zeotype.site_to_atom_indices = site_to_atom_indices
        zeotype.atom_indices_to_site = atom_indices_to_site
        return zeotype

    @staticmethod
    def read_cif_note_sites(fileobj, store_tags=False, primitive_cell=False,
                            subtrans_included=True, fractional_occupancies=True,
                            reader='ase'):
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
        possible_syms = get_available_symbols(b_dict["_atom_site_type_symbol"])  # only get symbols not in CIF file
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

    def update_nl(self):
        self.neighbor_list = NeighborList(natural_cutoffs(self), bothways=True, self_interaction=False)
        self.neighbor_list.update(self)

    def get_hetero_atoms(self, hetero_atoms_list=None) -> List[int]:
        """
        :return: Returns a list of all of the hetero-atoms in the zeotype
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
        :return: a dictionary where the key is the element symbol and the value is the number in the zeotype
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

    def add_cluster(self, index: int, max_size: int, max_neighbors: int) -> int:
        """
        Generates a Cluster of atoms around the specified index. The number of atoms in the cluster
        is given by the size parameter.
        :param max_size: max size of cluster
        :param index: index of the central atom in the cluster
        :param max_neighbors: number of neighbors from the host atom in the final cluster
        :return: index of the cluster in the zeotype cluster array
        """
        new_cluster = Cluster.build_from_zeolite(self, index, max_size, max_neighbors)
        self.clusters.append(new_cluster)
        return len(self.clusters) - 1  # returns final index in the clusters list

    def add_custom_cluster(self, cluster_indices: Iterable[int]):
        new_cluster = Cluster.build_from_zeolite(self, 0, 0, 0, cluster_indices=cluster_indices)
        self.clusters.append(new_cluster)
        return len(self.clusters) - 1  # returns final index in the clusters list

    def remove_cluster(self, index: int):
        """
        :param index: Index of cluster to remove from zeotype list
        :return: None
        """
        self.clusters[index] = None

    def find_silanol_groups(self):
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

    def find_silanol_nest_T_sites(self):
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

    def get_indices(self):
        return [a.index for a in self]

    def delete_atoms(self, indices_to_delete, iz_name='iz_1'):
        new_indices = [a.index for a in self if a.index not in indices_to_delete]
        iz = ImperfectZeotype(self[new_indices])
        iz.name = iz_name
        iz.parent_zeotype = self.parent_zeotype
        iz.add_iz_to_index_mapper()
        return iz

    def get_site_type(self, index):
        # pz_index = self[index_mapper.get_index(self.name, self.parent_zeotype.name, )
        pz_index = self.index_mapper.get_index(self, self.parent_zeotype, index)
        return self.parent_zeotype.atom_indices_to_site[pz_index]


class ImperfectZeotype(Zeotype):
    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, silent: bool = False, zeolite_type: str = '',
                 parent_zeotype=None, zeotype_to_cluster_index_map=None, neighbor_list=None, name='iz_1'):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms,
                         charges, scaled_positions, cell, pbc, celldisp, constraint,
                         calculator, info, velocities, silent, zeolite_type, zeotype_to_cluster_index_map,
                         neighbor_list, build_index_mapper=False)

        self.parent_zeotype = parent_zeotype
        self.index_mapper = self.parent_zeotype.index_mapper
        self.name = name

    def _get_pz_to_iz_map_by_pos(self):
        """
        Get the index mapping between a parent zeotype and imperfect zeotype.
        This matching is done by position, so it is essential that the atom
        positions have not changed after the iz was created.
        :return: a parent zeotype to imperfect zeotype index map
        """

        return self._get_old_to_new_map(self.parent_zeotype, self)

    @staticmethod
    def _get_old_to_new_map(old, new):
        """
        Get the index mapping between old and new self.
        his matching is done by position, so it is essential that the atom
        positions have not changed.
        """
        new_position_index_map = {}
        for atom in new:
            new_position_index_map[str(atom.position)] = atom.index

        old_index_position_map = {}
        for atom in old:
            old_index_position_map[atom.index] = str(atom.position)

        old_to_new_map = {}

        for key in old_index_position_map.keys():
            old_to_new_map[key] = new_position_index_map[old_index_position_map[key]]

        return old_to_new_map

    def add_iz_to_index_mapper(self):
        pz_to_iz_map = self._get_pz_to_iz_map_by_pos()
        self.index_mapper.add_name(self, self.name, self.parent_zeotype.name, pz_to_iz_map)

    def _delete_atoms(self, indices_to_delete):
        old_self = self.copy()
        for i in indices_to_delete:
            self.pop(i)
        self.index_mapper.delete_atoms(self.name, indices_to_delete)
        old_to_new_map = self._get_old_to_new_map(old_self, self)
        self.index_mapper.update_indices(self.name, old_to_new_map)

    def _add_atoms(self, atoms_to_add, atoms_name):
        old_self = self.copy()
        self.extend(atoms_to_add)
        new_atom_indices = list(set([a.index for a in self]) - set([a.index for a in old_self]))
        self.index_mapper.add_atoms(self.name, new_atom_indices)
        old_to_new_map = self._get_old_to_new_map(self, atoms_to_add)
        self.index_mapper.add_name(atoms_name, self.name, old_to_new_map)

    def

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

        for ads in cluster.adsorbates:
            # integrates all adsorbates, even non integrated ones in cluster
            ads_copy = ads.copy()
            ads_copy.host_zeotype = self
            ads_copy.integrate_ads()

    def integrate_adsorbate(self, adsorbate_index: int) -> None:
        self.adsorbates[adsorbate_index].integrate_ads()

    def cap_specific_atoms(self, indices_atom_to_cap, site_index):
        cap_atoms_dict = self.build_specific_cap_atoms_dict(indices_atom_to_cap, site_index)
        # append cap atoms self
        site_not_deleted = True
        for symbol, pos_list in cap_atoms_dict.items():
            for pos in pos_list:
                if site_not_deleted:  # replace first site with capping atom
                    self[site_index].symbol = symbol
                    self[site_index].position = pos
                    site_not_deleted = False
                else:
                    self.append(Atom(symbol, position=pos))

    def create_silanol_defect(self, site_index):
        atoms_to_cap = self.neighbor_list.get_neighbors(site_index)[0]
        self.cap_specific_atoms(atoms_to_cap, site_index)

    def get_hydrogen_cap_pos_site_dir(self, parent_zeolite, site_index, atom_to_be_capped_index):
        """

        :param parent_zeolite: The zeolite without the T site removed
        :param site_index: the index of the T site that has been removed
        :param atom_to_be_capped_index: the index of the atom where the cap is being added
        :return: The position of the new hydrogen atom
        """
        site_position = parent_zeolite[site_index].position
        direction = site_position - self.get_positions()[atom_to_be_capped_index]  # vector from neighbor to oxygen
        hydrogen_pos = self.get_positions()[atom_to_be_capped_index] + direction / np.abs(np.linalg.norm(direction))
        return hydrogen_pos

    def needs_cap(self, atom_index: int, bonds_needed: int) -> bool:
        """
        Finds if an atom needs a cap
        :param atom_index: index of atom to check for cap
        :param bonds_needed: number of bonds an atom needs
        :return: boolean stating if an atom needs a cap
        """
        return len(self.neighbor_list.get_neighbors(atom_index)[0]) < bonds_needed

    def build_specific_cap_atoms_dict(self, indices_atom_to_cap, site_index):
        """
        Builds a dictionary of the cap atom positions
        :param bonds_needed: a dict mapping atom symbol to number of bonds needed
        :return: a dictionary of atom symbols and positions of cap atoms
        """
        bonds_needed = {'O': 2, 'Si': 4, 'Sn': 4, 'Al': 4, 'Ga': 4, 'B': 4}
        cap_atoms_dict: Dict[str, List[int]] = defaultdict(list)
        self_copy = copy.deepcopy(self)
        self_copy.pop(site_index)
        self_copy.update_nl()
        indices, count = self_copy.count_elements()
        for si_index in indices['Si']:
            # tmp_si_index = si_index if si_index < site_index else site_index + 1  # because indices is messed up by del
            # if tmp_si_index not in indices_atom_to_cap:
            #     continue
            if self_copy.needs_cap(si_index, bonds_needed['Si']):
                for i in range(bonds_needed['Si'] - len(self_copy.neighbor_list.get_neighbors(si_index)[0])):
                    pos = self_copy.get_oxygen_cap_pos(si_index)
                    cap_atoms_dict['O'].append(pos)
                    self_copy.update_nl()

        for o_index in indices['O']:
            # tmp_o_index = o_index if o_index < site_index else o_index + 1  # because indices is messed up by del
            # if tmp_o_index not in indices_atom_to_cap:
            #     continue
            if self_copy.needs_cap(o_index, bonds_needed['O']):
                pos = self_copy.get_hydrogen_cap_pos_site_dir(self, site_index, o_index)
                cap_atoms_dict['H'].append(pos)
                self_copy.update_nl()

        return dict(cap_atoms_dict)

    def build_cap_atoms_dict(self, bonds_needed=None, hcap_fun=None):
        """
        Builds a dictionary of the cap atom positions
        :param bonds_needed: a dict mapping atom symbol to number of bonds needed
        :return: a dictionary of atom symbols and positions of cap atoms
        """
        if hcap_fun is None:
            hcap_fun = self.get_hydrogen_cap_pos

        cap_atoms_dict: Dict[str, List[int]] = defaultdict(list)
        if bonds_needed is None:
            bonds_needed = {'O': 2, 'Si': 4, 'Sn': 4, 'Al': 4, 'Ga': 4, 'B': 4}
        indices, count = self.count_elements()
        for si_index in indices['Si']:
            if self.needs_cap(si_index, bonds_needed['Si']):
                for i in range(bonds_needed['Si'] - len(self.neighbor_list.get_neighbors(si_index)[0])):
                    pos = self.get_oxygen_cap_pos(si_index)
                    cap_atoms_dict['O'].append(pos)
                    self.update_nl()

        for o_index in indices['O']:
            if self.needs_cap(o_index, bonds_needed['O']):
                pos = hcap_fun(o_index)
                cap_atoms_dict['H'].append(pos)
                self.update_nl()

        return dict(cap_atoms_dict)

    def get_oxygen_cap_pos(self, index):

        """
        Find a position of an oxygen cap
        :param index: index of atom needing cap
        :return:
        """
        # TODO: Remove bonds_needed from arg list

        neighbor = self.neighbor_list.get_neighbors(index)[0][-1]  # first index in the list of neighbor indicies
        direction = self.get_positions()[index] - self.get_positions()[neighbor]  # vector from neighbor to Si
        oxygen_pos = self.get_positions()[index] + 1.6 * direction / np.linalg.norm(direction)
        return oxygen_pos

    def get_hydrogen_cap_pos(self, index):
        """
        Finds the position of a hydrogen cap position
        :param index: index of hydrogen cap
        :return: the hydrogen position to add the cap too
        """
        neighbor = self.neighbor_list.get_neighbors(index)[0][0]  # first index in the list of neighbor indicies
        direction = self.get_positions()[index] - self.get_positions()[neighbor]  # vector from neighbor to oxygen
        hydrogen_pos = self.get_positions()[index] + direction / np.linalg.norm(direction)
        return hydrogen_pos


class Cluster(Zeotype):  # TODO include dynamic inheritance and

    def __init__(self, symbols=None, positions=None, numbers=None, tags=None, momenta=None, masses=None, magmoms=None,
                 charges=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None,
                 calculator=None, info=None, velocities=None, silent: bool = False, zeolite_type: str = '',
                 parent_zeotype=None, zeotype_to_cluster_index_map=None, neighbor_list=None):

        super().__init__(symbols, positions, numbers, tags, momenta, masses, magmoms,
                         charges, scaled_positions, cell, pbc, celldisp, constraint,
                         calculator, info, velocities, silent, zeolite_type, build_index_mapper=False)

        self.parent_zeotype = parent_zeotype
        self.zeotype_to_cluster_index_map = zeotype_to_cluster_index_map
        self.neighbor_list = neighbor_list

    @staticmethod
    def build_from_zeolite(parent_zeotype: Zeotype, index: int, max_cluster_size: int, max_neighbors: int,
                           cluster_indices: Iterable[int] = None) -> "Cluster":
        if cluster_indices is None:
            cluster_indices = Cluster._get_cluster_indices(parent_zeotype, index, max_cluster_size, max_neighbors)
        cluster_atoms = parent_zeotype[cluster_indices]
        new_cluster = Cluster(cluster_atoms)
        new_cluster.parent_zeotype = parent_zeotype
        new_cluster.zeotype_to_cluster_index_map = \
            new_cluster._get_new_cluster_mapping(new_cluster.parent_zeotype, cluster_atoms, cluster_indices)

        new_cluster.neighbor_list = NeighborList(natural_cutoffs(new_cluster), bothways=True, self_interaction=False)
        new_cluster.neighbor_list.update(new_cluster)

        return new_cluster

    def cap_atoms(self, cap_atoms_dict=None, bonds_needed=None, verbose=False, hcap_type='bond'):
        """each bare Si atom needs 4 Si atoms in oxygen and each of those oxygen needs two neighbors
        A lot of these clusters we add a bunch of H to the ends of the O because this a chemically plausable
        way to cap the oxygens. For example, when a zeolite is growing in solution it starts off as SiOH4
        1. Go through all Si in cluster and check to see if they have 4 bonds
        2.  If they don't have four bonds add an oxygen in a plausable direction
        3. go through all oxygens and find the ones that don't have two bonds
        4. If oxygens don't have two bonds add a hydrogen in a reasonable location (1 A distance)
        """

        if hcap_type == 'bond':
            hcap_fun = self.get_hydrogen_cap_pos
        else:
            hcap_fun = self.get_hydrogen_cap_pos_si_dir

        if cap_atoms_dict is None:
            cap_atoms_dict = self.build_cap_atoms_dict(bonds_needed=bonds_needed, hcap_fun=hcap_fun)
        if verbose:
            print('atom caps: {symbol, [position arrays]}', cap_atoms_dict)

        # append cap atoms self
        for symbol, pos_list in cap_atoms_dict.items():
            for pos in pos_list:
                self.append(Atom(symbol, position=pos))

    def get_hydrogen_cap_pos_si_dir(self, index):
        """
        Finds the position of a hydrogen cap position
        :param index: index of hydrogen cap
        :return: the hydrogen position to add the cap too
        """
        si_neighbor_position = self.find_si_neighbor(index)
        direction = si_neighbor_position - self.get_positions()[index]  # vector from neighbor to oxygen
        hydrogen_pos = self.get_positions()[index] + direction / np.abs(np.linalg.norm(direction))
        return hydrogen_pos

    @staticmethod
    def _get_cluster_indices(zeolite, index: int, max_size: int, max_neighbors: int) -> List[int]:
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
    def _get_new_cluster_mapping(zeolite: Zeotype, cluster: "Cluster", indices: List[int]):
        """
        Get the mapping between an unoptimized cluster and the zeolite from
        which the cluster was created
        :param zeolite: parent zeotype from which the cluster was created
        :param cluster: the child cluster that is being matched to the parrent zeotype
        :param indices: the indices in the parrent zeotype corresponding to the cluster
        :return: a zeotype to cluster index dictionary
        """
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

    def find_si_neighbor(self, cluster_index):
        inv_map = {v: k for k, v in self.zeotype_to_cluster_index_map.items()}
        zeotype_index = inv_map[cluster_index]
        for possible_index in self.parent_zeotype.neighbor_list.get_neighbors(zeotype_index)[0]:
            print(possible_index)
            if self.parent_zeotype[
                possible_index].symbol == 'Si' and possible_index not in self.zeotype_to_cluster_index_map.keys():
                return self.parent_zeotype.get_positions()[possible_index]

        # self.parent_zeotype.get_neighbors(index)

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

        for ads in cluster.adsorbates:
            # integrates all adsorbates, even non integrated ones in cluster
            ads_copy = ads.copy()
            ads_copy.host_zeotype = self
            ads_copy.integrate_ads()

    def integrate_adsorbate(self, adsorbate_index: int) -> None:
        self.adsorbates[adsorbate_index].integrate_ads()

    def cap_specific_atoms(self, indices_atom_to_cap, site_index):
        cap_atoms_dict = self.build_specific_cap_atoms_dict(indices_atom_to_cap, site_index)
        # append cap atoms self
        site_not_deleted = True
        for symbol, pos_list in cap_atoms_dict.items():
            for pos in pos_list:
                if site_not_deleted:  # replace first site with capping atom
                    self[site_index].symbol = symbol
                    self[site_index].position = pos
                    site_not_deleted = False
                else:
                    self.append(Atom(symbol, position=pos))

    def create_silanol_defect(self, site_index):
        atoms_to_cap = self.neighbor_list.get_neighbors(site_index)[0]
        self.cap_specific_atoms(atoms_to_cap, site_index)

    def get_hydrogen_cap_pos_site_dir(self, parent_zeolite, site_index, atom_to_be_capped_index):
        """

        :param parent_zeolite: The zeolite without the T site removed
        :param site_index: the index of the T site that has been removed
        :param atom_to_be_capped_index: the index of the atom where the cap is being added
        :return: The position of the new hydrogen atom
        """
        site_position = parent_zeolite[site_index].position
        direction = site_position - self.get_positions()[atom_to_be_capped_index]  # vector from neighbor to oxygen
        hydrogen_pos = self.get_positions()[atom_to_be_capped_index] + direction / np.abs(np.linalg.norm(direction))
        return hydrogen_pos

    def needs_cap(self, atom_index: int, bonds_needed: int) -> bool:
        """
        Finds if an atom needs a cap
        :param atom_index: index of atom to check for cap
        :param bonds_needed: number of bonds an atom needs
        :return: boolean stating if an atom needs a cap
        """
        return len(self.neighbor_list.get_neighbors(atom_index)[0]) < bonds_needed

    def build_specific_cap_atoms_dict(self, indices_atom_to_cap, site_index):
        """
        Builds a dictionary of the cap atom positions
        :param bonds_needed: a dict mapping atom symbol to number of bonds needed
        :return: a dictionary of atom symbols and positions of cap atoms
        """
        bonds_needed = {'O': 2, 'Si': 4, 'Sn': 4, 'Al': 4, 'Ga': 4, 'B': 4}
        cap_atoms_dict: Dict[str, List[int]] = defaultdict(list)
        self_copy = copy.deepcopy(self)
        self_copy.pop(site_index)
        self_copy.update_nl()
        indices, count = self_copy.count_elements()
        for si_index in indices['Si']:
            # tmp_si_index = si_index if si_index < site_index else site_index + 1  # because indices is messed up by del
            # if tmp_si_index not in indices_atom_to_cap:
            #     continue
            if self_copy.needs_cap(si_index, bonds_needed['Si']):
                for i in range(bonds_needed['Si'] - len(self_copy.neighbor_list.get_neighbors(si_index)[0])):
                    pos = self_copy.get_oxygen_cap_pos(si_index)
                    cap_atoms_dict['O'].append(pos)
                    self_copy.update_nl()

        for o_index in indices['O']:
            # tmp_o_index = o_index if o_index < site_index else o_index + 1  # because indices is messed up by del
            # if tmp_o_index not in indices_atom_to_cap:
            #     continue
            if self_copy.needs_cap(o_index, bonds_needed['O']):
                pos = self_copy.get_hydrogen_cap_pos_site_dir(self, site_index, o_index)
                cap_atoms_dict['H'].append(pos)
                self_copy.update_nl()

        return dict(cap_atoms_dict)

    def build_cap_atoms_dict(self, bonds_needed=None, hcap_fun=None):
        """
        Builds a dictionary of the cap atom positions
        :param bonds_needed: a dict mapping atom symbol to number of bonds needed
        :return: a dictionary of atom symbols and positions of cap atoms
        """
        if hcap_fun is None:
            hcap_fun = self.get_hydrogen_cap_pos

        cap_atoms_dict: Dict[str, List[int]] = defaultdict(list)
        if bonds_needed is None:
            bonds_needed = {'O': 2, 'Si': 4, 'Sn': 4, 'Al': 4, 'Ga': 4, 'B': 4}
        indices, count = self.count_elements()
        for si_index in indices['Si']:
            if self.needs_cap(si_index, bonds_needed['Si']):
                for i in range(bonds_needed['Si'] - len(self.neighbor_list.get_neighbors(si_index)[0])):
                    pos = self.get_oxygen_cap_pos(si_index)
                    cap_atoms_dict['O'].append(pos)
                    self.update_nl()

        for o_index in indices['O']:
            if self.needs_cap(o_index, bonds_needed['O']):
                pos = hcap_fun(o_index)
                cap_atoms_dict['H'].append(pos)
                self.update_nl()

        return dict(cap_atoms_dict)

    def get_oxygen_cap_pos(self, index):

        """
        Find a position of an oxygen cap
        :param index: index of atom needing cap
        :return:
        """
        # TODO: Remove bonds_needed from arg list

        neighbor = self.neighbor_list.get_neighbors(index)[0][-1]  # first index in the list of neighbor indicies
        direction = self.get_positions()[index] - self.get_positions()[neighbor]  # vector from neighbor to Si
        oxygen_pos = self.get_positions()[index] + 1.6 * direction / np.linalg.norm(direction)
        return oxygen_pos

    def get_hydrogen_cap_pos(self, index):
        """
        Finds the position of a hydrogen cap position
        :param index: index of hydrogen cap
        :return: the hydrogen position to add the cap too
        """
        neighbor = self.neighbor_list.get_neighbors(index)[0][0]  # first index in the list of neighbor indicies
        direction = self.get_positions()[index] - self.get_positions()[neighbor]  # vector from neighbor to oxygen
        hydrogen_pos = self.get_positions()[index] + direction / np.linalg.norm(direction)
        return hydrogen_pos


# testing
if __name__ == '__main__':
    from ase.io import read

    b = read('/Users/dda/Code/zeotype/data/sam/sn_ti-periodic.traj')
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
