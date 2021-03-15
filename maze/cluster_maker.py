from abc import ABC
from maze.perfect_zeotype import Zeotype
from maze.zeotypes import ImperfectZeotype
from typing import Iterable
from ase.neighborlist import natural_cutoffs, NeighborList
from typing import List


class ClusterMaker(ABC):
    def get_cluster_indices(self, zeotype: Zeotype, start_site: int, **kwargs):
        raise NotImplementedError

    @staticmethod
    def get_cluster(zeotype: Zeotype, indices: Iterable[int], name="Cluster"):
        cluster = ImperfectZeotype(zeotype, ztype=name)
        to_delete = cluster.get_indices_compliment(cluster, indices)
        cluster = cluster.delete_atoms(to_delete)
        return cluster

    @staticmethod
    def get_open_defect(zeotype: Zeotype, indices: Iterable[int], name="Open Defect"):
        new_od = ImperfectZeotype(zeotype, ztype=name)
        new_od = new_od.delete_atoms(indices_to_delete)
        return new_od


class DefaultClusterMaker(ClusterMaker):
    def __init__(self):
        # note you can change this function to be whatever you want
        self.get_cluster_indices = self.get_oh_cluster_indices

    def get_cluster_indices(self, zeotype: Zeotype, start_site: int, **kwargs):
        return self.get_cluster_indices(zeotype, start_site)

    @staticmethod
    def get_oh_cluster_multi_t_sites(zeolite: Zeotype, t_sites: Iterable[int]) -> List[int]:
        """
        get an OH cluster with multiple T sites
        :param zeolite: The MAZE-sim from which to extract the cluster
        :param t_sites: the central t site
        :return: A list of indices of the cluster
        """
        all_indces = set()
        for t_site in t_sites:
            all_indces.update(Cluster.get_oh_cluster_indices(zeolite, t_site))

        return list(all_indces)

    @staticmethod
    def get_oh_cluster_indices(zeolite: Zeotype, t_site: int) -> List[int]:
        """
        Create a cluster that only includes one central T site and then Oxygen
        and Hydrogen atoms. This is different than the other cluster selection
        methods that take in other
        :param zeolite: The zeolite from which the cluster indices will be drawn
        :param t_site: The index of the T site around which the cluster will be built
        :return: The indices of the new cluster
        """

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
