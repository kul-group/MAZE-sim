from itertools import count
from typing import Iterable


class IndexMapper:
    """
    This class maps between the different atoms objects in a zeotype project.
    It is essential for keeping track of the various
    """

    def __init__(self, name: str, atom_indices: Iterable):
        """

        :param name: the str name of the atom-like object being used to initialize the Index Mapper.
        This should be a zeotype object
        :param atom_indices: A list of atom indices
        """
        self.main_index = {}
        self.i_max = 0
        self.names = [name]
        for main_index, atom_index in zip(count(), atom_indices):
            self.main_index[main_index] = {name: atom_index}
            self.i_max = main_index

    def _reverse_main_index(self, name):
        """

        :param name:
        :return:
        """
        name_to_main_dict = {}
        for main_index, value in self.main_index.items():
            name_index = value[name]
            if name_index is None:  # ignore none indices
                continue
            name_to_main_dict[name_index] = main_index

        return name_to_main_dict

    def add_name(self, new_name, old_name, old_name_to_new_name, new_atom_indices=None):
        """
        If new_name already exist, this will still work well
        :param new_name:
        :param old_name:
        :param old_name_to_new_name:
        :param new_atom_indices:
        :return:
        """
        if new_name not in self.names:
            self.names.append(new_name)

        old_name_to_main_dict = self._reverse_main_index(old_name)
        main_to_new_name_dict = {}
        for old_ind, main_ind in old_name_to_main_dict.items():
            main_to_new_name_dict[main_ind] = old_name_to_new_name.get(old_ind, None)

        for i in self.main_index.keys():
            self.main_index[i][new_name] = main_to_new_name_dict.get(i, None)

        if new_atom_indices is not None:
            self.add_atoms(new_name, new_atom_indices)

    def make_none_dict(self):
        none_dict = {}
        for name in self.names:
            none_dict[name] = None
        return none_dict

    def add_atoms(self, name, new_atom_indices):
        for index in new_atom_indices:
            none_dict = self.make_none_dict()  # could be slow
            none_dict[name] = index
            self.i_max += 1
            self.main_index[self.i_max] = none_dict

    def delete_atoms(self, name, atom_indices_to_delete):
        name_to_main_dict = self._reverse_main_index(name)
        for i in atom_indices_to_delete:
            self.main_index[name_to_main_dict[i]][name] = None

    def update_indices(self, name, old_to_new_map):
        for index, value in self.main_index.items():
            old_value = value[name]
            if old_value in old_to_new_map:
                value[name] = old_to_new_map[old_value]
            else:
                value[name] = None

    def get_index(self, sender_name, receiver_name, sender_index):
        for name_dict in self.main_index.values():
            my_name_dict = name_dict[sender_name]
            my_sender_index = sender_index
            if name_dict[sender_name] == sender_index:
                return name_dict[receiver_name]

    def get_overlap(self, name1, name2):
        overlap_indices_name1 = []
        for name_dict in self.main_index.values():
            if name_dict[name1] is not None and name_dict[name2] is not None:
                overlap_indices_name1.append(name_dict[name1])

        return overlap_indices_name1
