from source.zeotype import Zeotype
from ase.visualize import view


def main():
    cif_dir = "/Users/dda/Code/zeotype/data/GOO.cif"
    zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    # make cluster
    ci = zeolite.add_cluster(95, 6)
    cluster = zeolite.clusters[ci]
    view(zeolite)
    view(cluster)
    # cap atoms
    cluster.cap_atoms(verbose=True)
    view(cluster)

    # change all cluster atoms to Hg
    for i in range(len(cluster)):
        cluster[i].symbol = 'Hg'
    view(cluster)
    view(zeolite)

    # now integrate in cluster
    zeolite.integrate_cluster(ci)
    view(zeolite)


if __name__ == '__main__':
    main()
