from ase.io import read, write
from source.zeotype import Zeotype
import os
from ase.visualize import view

def main():
    #data_dir = os.path.join(os.getcwd(), 'data', 'sam')
    #filename = "sn_ti-8T-cluster.traj"
    # filepath = os.path.join(data_dir, filename)
    data_dir = '/Users/dda/Code/zeotype/data'
    filepath1 = "/Users/dda/Code/zeotype/data/sam/sn_ti-8T-cluster.traj"
    filepath2 = "/Users/dda/Code/zeotype/data/sam/sn_ti-periodic.traj"
    view(read(filepath1))
    z = Zeotype(read(filepath2))
    cluster_indices = [191,
                       168, 120, 84, 96, 112,
                       159, 79, 63, 87, 35,
                       175, 119, 127, 82, 103,
                       72,
                       152,
                       136, 28, 40, 36, 20,
                       172, 116, 80, 124, 100,
                       144, 64, 0, 48, 64]

    cluster_index = z.add_custom_cluster(cluster_indices)
    z.clusters[cluster_index].cap_atoms()
    view(z.clusters[cluster_index])
    view(z)
    print(z.clusters[cluster_index].zeotype_to_cluster_index_map)
    write(os.path.join(data_dir, 'dexter_t8_cluster.traj'), z.clusters[cluster_index])

    # z = Zeotype(a)
    # z2 = Zeotype(read(filepath2))
    # view(z2)
    # view(z)


if __name__ == "__main__":
    main()