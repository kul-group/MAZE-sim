import copy
from collections import defaultdict
from typing import Dict, List

import numpy as np
from ase import Atom

from source.zeotypes import Zeotype, ImperfectZeotype


class ExtraFramework(ImperfectZeotype):
    def find_si_neighbor(self, cluster_index):
        inv_map = {v: k for k, v in self.zeotype_to_cluster_index_map.items()}
        zeotype_index = inv_map[cluster_index]
        for possible_index in self.parent_zeotype.neighbor_list.get_neighbors(zeotype_index)[0]:
            if self.parent_zeotype[possible_index].symbol == 'Si' and\
                    possible_index not in self.zeotype_to_cluster_index_map.keys():
                return self.parent_zeotype.get_positions()[possible_index]




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
