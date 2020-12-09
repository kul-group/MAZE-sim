from itertools import count
class IndexMapper:
    """
    This class maps between the different atoms objects in a zeotype project.
    It is essential for
    """
    def __init__(self, name: str, atom_indices):
        self.main_index = {}
        self.i_max = 0
        self.names = [name]
        for main_index, atom_index in zip(count(), atom_indices):
            self.main_index[main_index] = {name: atom_index}
            self.i_max = main_index

    def reverse_main_index(self, name):
        name_to_main_dict = {}
        for main_index, value in self.main_index.items():
            name_index = value[name]
            if name_index is None:  # ignore none indices
                continue
            name_to_main_dict[name_index] = main_index

        return name_to_main_dict

    def add_name(self, new_name, old_name, old_name_to_new_name, new_atom_indices=None):
        self.names.append(new_name)
        old_name_to_main_dict = self.reverse_main_index(old_name)
        main_to_new_name_dict = {}
        for old_ind, main_ind in old_name_to_main_dict.items():
            main_to_new_name_dict[main_ind] = old_name_to_new_name.get(old_ind, default=None)

        for i in self.main_index.keys():
            self.main_index[i][new_name] = main_to_new_name_dict.get(i, default=None)

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

    def get_index(self, sender_name, receiver_name, sender_index):
        for name_dict in self.main_index.values():
            if name_dict[sender_name] is sender_index:
                return name_dict[receiver_name]

