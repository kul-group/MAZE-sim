from source.zeotype import Zeotype, ImperfectZeotype, Cluster
from ase.visualize import view


def main():
    cif_dir = "/Users/dda/Code/zeotype/data/BEA.cif"
    zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    c_in = Cluster.get_oh_cluster_indices(zeolite, 185)
    print(c_in)
    cluster, open_framework = zeolite.get_cluster(0,0,0,cluster_indices=c_in)

    view(cluster)

    #cluster = cluster.cap_atoms()

    # view(cluster)
    # view(open_framework.cap_atoms())
    #view(open_framework)

    # open_framework = open_framework.remove_caps()
    # view(open_framework)
    # open_framework = open_framework.integrate_other_zeotype(cluster)
    # view(open_framework)

if __name__ == '__main__':
    main()
