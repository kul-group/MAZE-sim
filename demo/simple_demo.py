import os
from source.zeotype import Zeotype, ImperfectZeotype
from ase.visualize import view
from pathlib import Path



def main():
    data_dir = os.path.join(Path(os.getcwd()).parent, 'data')
    cif_dir = os.path.join(data_dir, 'GOO.cif')
    zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    print(isinstance(  zeolite, Zeotype))


if __name__ == "__main__":
    main()
