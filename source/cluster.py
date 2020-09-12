from .zeotype import Zeotype


class Cluster(Zeotype):
    def __init__(self, parent_zeotype: Zeotype, index: int, cluster_size: int):
        self.parent_zeotype = parent_zeotype

        cluster_indices = self._get_cluster_indices(self.parent_zeotype, index, cluster_size)
        cluster_atoms = self.parent_zeotype[cluster_indices]
        self.zeotype_to_cluster_index_map = \
            self._get_new_cluster_mapping(self.parent_zeotype, cluster_atoms, cluster_indices)

        pz = self.parent_zeotype
        super().__init__(self, pz.symbols, pz.positions, pz.numbers, pz.tags, pz.momenta, pz.masses, pz.magmoms,
                         pz.charges, pz.scaled_positions, pz.cell, pz.pbc, pz.celldisp, pz.constraint,
                         pz.calculator, pz.info, pz.velocities, pz.silent, pz.zeolite_type)

    def cap_atoms(self):
        ...

    @staticmethod
    def _get_cluster_indices(zeolite, index: int, size: int):
        nl = NeighborList(natural_cutoffs(zeolite), self_interaction=False, bothways=True)
        nl.update(zeolite)

        cluster_indices = set()
        new_cluster_indices = set([index])

        while True:
            current_cluster_indices = set()
            for cluster_index in new_cluster_indices:
                cluster_indices.add(cluster_index)
                if len(cluster_indices) >= size:
                    return list(cluster_indices)
                for new_index in nl.get_neighbors(cluster_index)[0]:
                    current_cluster_indices.add(new_index)
            new_cluster_indices = current_cluster_indices

    @staticmethod
    def _get_new_cluster_mapping(zeolite, cluster, indices: List[int]):
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

