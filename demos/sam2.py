from ase.io import read
from maze import Zeotype, Cluster
import os
from ase.visualize import view
from pathlib import Path

def main():
    data_dir = os.path.join(Path(os.getcwd()).parent, 'data', 'sam')

    #sn_ti_traj_filepath = os.path.join(data_dir, "sn_ti-periodic.traj")
    #cluster_traj_filepath = os.path.join(data_dir, "sn_ti-8T-cluster.traj")
    vasp = '/Users/dda/Downloads/sn_bea_tmpo_structures/02_open2/03_Sn/opt_from_vasp.traj'
    #pc = Zeotype(read(cluster_traj_filepath)) #perfect cluster
    z = Zeotype(read(vasp))
    view(z)
    #cluster_indices = list(set(Cluster.get_oh_cluster_indices(z, 191)).union(Cluster.get_oh_cluster_indices(z, 152)))#Cluster.get_cluster_indices_multi_T_site(z, [191, 152], 1000, 3)
    cluster_indices = Cluster.get_oh_cluster_multi_t_sites(z, [191])
    cluster, od = z.get_cluster(0,0,0, cluster_indices=cluster_indices)
    cluster = cluster.cap_atoms()
    view(cluster)
    # view(pc)
    #cluster = cluster.cap_atoms()
    #od = od.cap_atoms()
    #od = od.cap_atoms()
    #view(od)
    #pass
    #view(z)
    #view(cluster)

    #view(od)
    #view(cluster)
    #view(pc)
    # cluster_indices = [191,
    #                    168, 120, 84, 96, 112,
    #                    159, 79, 63, 87, 35,
    #                     175, 119, 127, 82, 103,
    #                    72,
    #                    152,
    #                    136, 28, 40, 36, 20,
    #                    172, 116, 80, 124, 100,
    #                    144, 64, 0, 48, 64]
    #
    # cluster_index = z.add_custom_cluster(cluster_indices)
    # z.clusters[cluster_index].cap_atoms(hcap_type='si')
    # view(z.clusters[cluster_index])
    # #view(z)
    # write(os.path.join(data_dir, 'sam_t8_cluster.traj'), z.clusters[cluster_index])
    # write(os.path.join(data_dir, 'overlay.traj'), z)
    # for a in z.clusters[cluster_index]:
    #     if a.index not in z.clusters[cluster_index].zeotype_to_cluster_index_map:
    #         z.append(a)
    # view(z)data/sam/
if __name__ == "__main__":
    main()

