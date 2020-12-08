from .zeotype import Zeotype
class ExtraFramework(Zeotype)
    ...

    def find_si_neighbor(self, cluster_index):
        inv_map = {v: k for k, v in self.zeotype_to_cluster_index_map.items()}
        zeotype_index = inv_map[cluster_index]
        for possible_index in self.parent_zeotype.neighbor_list.get_neighbors(zeotype_index)[0]:
            print(possible_index)
            if self.parent_zeotype[possible_index].symbol == 'Si' and possible_index not in self.zeotype_to_cluster_index_map.keys():
                return self.parent_zeotype.get_positions()[possible_index]

        #self.parent_zeotype.get_neighbors(index)