import os
from maze.zeolite import PerfectZeolite
from maze.io_zeolite import *

if __name__ == "__main__":
    data_path = '/Users/dda/Code/MAZE-sim/data/'
    cif_path = os.path.join(data_path, 'BEA.cif')
    output_path = os.path.join(data_path, 'my_first_zeotype')
    my_zeotype = PerfectZeolite.build_from_cif_with_labels(cif_path)
    print(my_zeotype.index_mapper.main_index)
    cluster, od = my_zeotype.get_cluster(1,10,10)
    # data_dir = os.path.join(Path(os.getcwd()).parent, 'data')
    # output_traj = os.path.join(data_dir, 'test.traj')
    # cif_dir = os.path.join(data_dir, 'BEA.cif')
    # zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    save_zeotypes(output_path, [my_zeotype, cluster, od])
    new_dict = read_zeotypes(output_path + '.zeo')
    print(new_dict)
    print(new_dict['parent'])
    print(new_dict['parent'].index_mapper.main_index)
    # Trajectory(output_traj,'w', zeolite)
