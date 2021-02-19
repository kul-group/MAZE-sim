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