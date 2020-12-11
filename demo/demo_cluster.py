from source.zeotype import Zeotype, ImperfectZeotype
from ase.visualize import view


def main():
    cif_dir = "/Users/dda/Code/zeotype/data/GOO.cif"
    zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    #iz = ImperfectZeotype(zeolite)
    #iz.delete_atoms([i for i in range(0, 90)])
    # make cluster
    cluster, open_framework = zeolite.add_cluster(84, 400, 3)
    view(cluster)
    view(open_framework)


if __name__ == '__main__':
    main()
