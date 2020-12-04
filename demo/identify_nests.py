from source.zeotype import Zeotype
from ase.visualize import view
from ase.io import read
import os
from pathlib import Path

def main():
    data_dir = os.path.join(Path(os.getcwd()).parent, 'data', 'sam')
    sn_8T_cluster = os.path.join(data_dir, "sn_ti-8T-cluster.traj")
    z = Zeotype(read(sn_8T_cluster))
    view(z)
    si_list = z.find_silanol_groups()
    print(si_list)
    print(z.find_silanol_nest_T_sites())

if __name__ == "__main__":
    main()
