from source.zeotype import Zeotype, ImperfectZeotype, Cluster
from ase.visualize import view


def main():
    cif_dir = "/Users/dda/Code/zeotype/data/BEA.cif"
    zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    view(zeolite)
    c_in = Cluster.get_oh_cluster_indices(zeolite, 185)
    cluster, open_framework = zeolite.get_cluster(cluster_indices=c_in)
    cluster = cluster.cap_atoms()
    open_framework = open_framework.cap_atoms()
    # print(cluster.additions)
    # cluster = cluster.remove_caps('h_caps', 'h_caps_6_')
    # print(cluster.name)
    # print(cluster.index_mapper.main_index)
    # open_framework = open_framework.cap_atoms()
    # #view(open_framework)
    # view(cluster)
    #
    # pass
    # #cluster = cluster.cap_atoms()

    # view(cluster)
    # view(open_framework.cap_atoms())
    # view(open_framework)

    open_framework = open_framework.remove_caps()
    # view(open_framework)
    # open_framework = open_framework.integrate_other_zeotype(cluster)
    # view(open_framework)


if __name__ == '__main__':
    main()
