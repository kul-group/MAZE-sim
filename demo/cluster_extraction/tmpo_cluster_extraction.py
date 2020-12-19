from ase.io import read, write
from source.zeotype import Zeotype, Cluster
import os
from ase.visualize import view
from pathlib import Path
from ase.neighborlist import natural_cutoffs, NeighborList
import glob

def find_tmpo(atoms):
    """
    Args:
        atoms:
    """
    tmpo_indices = []
    p_index = None
    for a in atoms:
        if a.symbol == 'P':
            p_index = a.index
            break
    tmpo_indices.append(p_index)

    nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
    nl.update(atoms)
    p_nl = nl.get_neighbors(p_index)[0]
    tmpo_indices.extend(p_nl)
    for i in p_nl:
        if atoms[i].symbol == 'C':
            tmpo_indices.extend(nl.get_neighbors(i)[0])

    return tmpo_indices


def main():
    input_data_dir = '/Users/dda/Box/Kulkarni/data/sn_bea_tmpo_structures/'
    output_data_dir = '/Users/dda/Box/Kulkarni/data/sn_bea_tmpo_structures/output'
    glob_cmd = os.path.join(input_data_dir, '**/**/*.traj')
    traj_files = glob.glob(glob_cmd)
    for traj_file in traj_files:
        z = Zeotype(read(traj_file))
        tmpo = find_tmpo(z)
        tin_index = [a.index for a in z if a.symbol == 'Sn'][0]
        cluster_indices = Cluster.get_oh_cluster_multi_t_sites(z, [tin_index])
        cluster_indices = list(set(cluster_indices).union(set(tmpo)))
        cluster, od = z.get_cluster(0,0,0, cluster_indices=cluster_indices)
        cluster = cluster.cap_atoms()
        output_dir_path = os.path.join(output_data_dir,
                                       traj_file.split(input_data_dir)[-1].split(os.path.basename(traj_file))[0])
        output_filepath = os.path.join(output_dir_path, 'cluster.traj')
        Path(output_dir_path).mkdir(exist_ok=True, parents=True)
        write(output_filepath, cluster)
        print("wrote cluster to ", output_filepath)

if __name__ == "__main__":
    main()

