from itertools import count
from typing import Iterable, Dict, List


class IndexMapper:
    """
    This class maps between the different atoms objects in a maze project.
    It is essential for keeping track of the various
    """

    id = 0

    @staticmethod
    def get_id():
        """
        Get a unique id
        :return:
        """
        IndexMapper.id += 1
        return str(IndexMapper.id)

    @staticmethod
    def get_unique_name(name: str):
        """
        Get a unique name
        :param name: name
        :return: name + _ + unique id number
        """

        return name + '_' + IndexMapper.get_id()

    def __init__(self, atom_indices: Iterable) -> None:
        """
        :param atom_indices: A list of atom indices from a Zeotype
        """

        self.main_index = {}
        self.i_max = 0
        self.names = ['parent']
        for main_index, atom_index in zip(count(), atom_indices):
            self.main_index[main_index] = {'parent': atom_index}
            self.i_max = main_index

    def get_reverse_main_index(self, name: str) -> Dict[int, int]:
        """
        Reverses an index so that the name indices are used as a key for the main index
        :param name: name of the item to make the key
        :return: a reverse index map where the indices of name are the key
        """
        return self._reverse_main_index(name)

    def _reverse_main_index(self, name: str) -> Dict[int, int]:
        """
        Internal function for making a reverse index map
        :param name: name of the item to make the key
        :return:a reverse index map where the indices of name are the key
        """

        name_to_main_dict = {}
        for main_index, value in self.main_index.items():
            try:
                name_index = value[name]
            except KeyError:
                print(f'name not found at index {main_index}')
                continue

            if name_index is None:  # ignore none indices
                continue
            name_to_main_dict[name_index] = main_index

        return name_to_main_dict

    def get_name1_to_name2_map(self, name1: str, name2: str) -> Dict[int, int]:
        """
        Gets a map between the indices of name1 and the indices of name2
        :param name1: name whose indices are the key in the map
        :param name2: name whose indices are the value in the map
        :return: name1.index -> name2.index map
        """

        name1_to_name2_map = {}
        for row in self.main_index.values():
            if row[name1] is not None:
                name1_to_name2_map[name1] = row[name2]

        return name1_to_name2_map

    def register_with_main(self, new_name: str, main_to_new_map: Dict[int, int]) -> None:
        """
        Register a new object by using the main index in mapping.
        :param new_name: name of the new object being registered
        :param main_to_new_map:
        :return:
        """

        self.names.append(new_name)
        for main_i in self.main_index.keys():
            self.main_index[main_i][new_name] = main_to_new_map.get(main_i, None)

    def register(self, old_name: str, new_name: str, old_to_new_map: Dict[int, int]) -> None:
        """
        Register a new object with the indexer
        :param old_name: name of object known by the indexer
        :param new_name: name of new object being registered
        :param old_to_new_map: a index mapping between the old map and the new map
        note that additional atoms in new will be added to the index
        :return:
        """
        assert new_name not in self.names, f'Error: {new_name} has already been registered'
        assert old_name in self.names, f'Error: {old_name} has not been registered'
        self.names.append(new_name)

        old_name_to_main_dict = self._reverse_main_index(old_name)
        main_to_new_name_dict = {}
        for old_ind, main_ind in old_name_to_main_dict.items():
            main_to_new_name_dict[main_ind] = old_to_new_map.get(old_ind, None)

        for i in self.main_index.keys():
            self.main_index[i][new_name] = main_to_new_name_dict.get(i, None)

    def add_name(self, new_name: str, old_name: str, old_name_to_new_name: Dict[int, int],
                 new_atom_indices: Iterable[int] = None) -> None:
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

    def _make_none_dict(self) -> Dict[str, None]:
        """
        Get a dictionary full of None for each name
        :return: A dictionary full of None for each known name
        """
        none_dict = {}
        for name in self.names:
            none_dict[name] = None
        return none_dict

    def add_atoms(self, name: str, new_atom_indices: Iterable[int]) -> None:
        """
        Add new atoms to the dictionary
        :param name: name of object with new atoms
        :param new_atom_indices: indices of the new atoms
        :return:
        """
        assert name not in self.names, 'name already exists'
        for index, value in self.main_index.items():
            value[name] = None
        self.names.append(name)

        for index in new_atom_indices:
            none_dict = self._make_none_dict()  # could be slow
            none_dict[name] = index
            self.i_max += 1
            self.main_index[self.i_max] = none_dict

        x = self.main_index

    def delete_atoms(self, name: str, atom_indices_to_delete: Iterable[int]) -> None:
        """

        :param name: name of with indices to delete
        :param atom_indices_to_delete: indices to delete
        :return: None
        """

        name_to_main_dict = self._reverse_main_index(name)
        for i in atom_indices_to_delete:
            self.main_index[name_to_main_dict[i]][name] = None

    def update_indices(self, name: str, old_to_new_map: Dict[int, int]) -> None:
        """
        Update the indices of an item
        :param name: name of indices to update
        :param old_to_new_map: old map to new map
        :return: None
        """

        for index, value in self.main_index.items():
            old_value = value[name]
            if old_value in old_to_new_map:
                value[name] = old_to_new_map[old_value]
            else:
                value[name] = None

    def get_index(self, sender_name: str, receiver_name: str, sender_index: int) -> int:
        """
        get the index of another object
        :param sender_name: name of the sender zeolite
        :param receiver_name: the name of the receiving zeolite
        :param sender_index: the index of the sender
        :return: the receiving index
        """

        for name_dict in self.main_index.values():
            my_name_dict = name_dict[sender_name]
            my_sender_index = sender_index
            if name_dict[sender_name] == sender_index:
                return name_dict[receiver_name]

    def get_overlap(self, name1: str, name2: str) -> List[int]:
        """
        Get the list of names that overlap
        :param name1: name of object 1
        :param name2: name of object 2
        :return: overlapping indices
        """
        overlap_indices_name1 = []
        for name_dict in self.main_index.values():
            if name_dict[name1] is not None and name_dict[name2] is not None:
                overlap_indices_name1.append(name_dict[name1])

        return overlap_indices_name1

    def delete_name(self, name: str) -> None:
        """
        Delete zeolite from the index
        :param name: The name of the zeolite to delete from the index
        :return: None
        """
        try:
            self.names.remove(name)
        except:
            print(name)  # TODO: Write custom error message
            return None
        for index, old_row in self.main_index.items():
            new_row = self._make_none_dict()
            for key in new_row.keys():
                new_row[key] = old_row[key]

            self.main_index[index] = new_row