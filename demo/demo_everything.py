from source.zeotype import Zeotype, ImperfectZeotype, Cluster
from ase.visualize import view

def main():
    cif_dir = "/Users/dda/Code/zeotype/data/BEA.cif"
    zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    c_in = Cluster.get_oh_cluster_indices(zeolite, 185)
    cluster, od = zeolite.get_cluster(cluster_indices=c_in)
    view(cluster)
    cluster = cluster.cap_atoms()
    print(cluster.additions)
    view(cluster)
    cluster = cluster.remove_caps('h_caps', 'h_caps_6_')
    print(cluster.additions)
    view(cluster)
    for a in cluster:
        cluster[a.index].symbol = 'Xe'
    view(cluster)
    repair = od.integrate_other_zeotype(cluster)
    view(repair)
    view(od.create_silanol_defect(164))


if __name__ == '__main__':
    main()
