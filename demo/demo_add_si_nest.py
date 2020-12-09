import os
from source.zeotype import Zeotype
from ase.visualize import view
from patlib import path
def main():
    data_dir = os.path.join(Path(os.getcwd()).parent, 'data')
    cif_dir = os.path.join(data_dir, 'GOO.cif')
    zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    T1_index = 95
    zeolite[T1_index].symbol = 'Hg'
    view(zeolite)
    zeolite.create_silanol_defect(T1_index)
    view(zeolite)


if __name__ == '__main__':
    main()
