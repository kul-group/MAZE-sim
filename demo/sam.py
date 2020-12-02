from ase.io import read, write
from source.zeotype import Zeotype
import os
from ase.visualize import view
from pathlib import Path

def main():
    data_dir = os.path.join(Path(os.getcwd()).parent, 'data', 'sam')

    sn_ti_traj_filepath = os.path.join(data_dir, "sn_ti-periodic.traj")
    z = Zeotype(read(sn_ti_traj_filepath))

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
    write(os.path.join(data_dir, 'sam_t8_cluster.traj'), z.clusters[cluster_index])



if __name__ == "__main__":
    main()