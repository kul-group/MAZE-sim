from typing import List


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

    @staticmethod
    def find_silanol_groups(zeotype: "Zeotype") -> List["Silanol"]:
        """
        Finds all of the silanol groups in the Zeotype

        :return: A list of Silanol groups
        """
        silanol_list = []
        for atom in zeotype:
            if atom.symbol == 'Si':
                for neighbor_index in zeotype.neighbor_list.get_neighbors(atom.index)[-1]:
                    if zeotype[neighbor_index].symbol == 'O':
                        for next_neighbor_index in zeotype.neighbor_list.get_neighbors(zeotype[neighbor_index].index)[
                            -1]:
                            if zeotype[next_neighbor_index].symbol == 'H':
                                # found one!
                                silanol = Silanol(zeotype, atom.index, neighbor_index, next_neighbor_index,
                                                  zeotype.neighbor_list.get_neighbors(atom.index)[-1])
                                silanol_list.append(silanol)
        return silanol_list

    @staticmethod
    def find_silanol_nest_T_sites(zeotype: "Zeotype") -> List[int]:
        """
        Finds all of the T sites that are in silanol nests

        :return: A list of T sites in silanol nests
        """
        sites_list = []
        sil_list = Silanol.find_silanol_groups(zeotype)
        if zeotype._atom_indices_to_site is None:
            for sil in sil_list:
                sites_list.append(sil.Si_index)
                for i in sil.Si_neighbor_list:
                    if 'Sn' == zeotype[i].symbol:
                        sites_list.append(i)
        else:
            for sil in sil_list:

                if 'T' in zeotype._atom_indices_to_site(sil.Si_index):
                    sites_list.append(sil.Si_index)
                else:
                    for index in sil.Si_neighbor_list:
                        if 'Sn' == zeotype[index].symbol:
                            sites_list.append(index)

        return sites_list
