from source.zeotype import Zeotype
from ase.visualize import view

def main():
    cif_dir = "/Users/dda/Code/zeotype/data/GOO.cif"
    zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    T1_index = 95
    zeolite[T1_index].symbol = 'Hg'
    view(zeolite)
    zeolite.create_silanol_defect(T1_index)
    view(zeolite)


if __name__ == '__main__':
    main()
