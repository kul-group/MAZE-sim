import os
from maze import Zeotype
from pathlib import Path

def main():
    data_dir = os.path.join(Path(os.getcwd()).parent, 'data')
    cif_dir = os.path.join(data_dir, 'GOO.cif')
    zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    print(issubclass(Zeotype, zeolite))
    #iz = Zeolite(zeolite)
    # iz.name = "fred"
    # iz.parent_zeotype = zeolite
    # iz.index_mapper = iz.parent_zeotype.index_mapper
    # iz.add_iz_to_index_mapper()
    #
    # T1_index = 95
    # zeolite[T1_index].symbol = 'Hg'
    # view(zeolite)
    # iz.create_silanol_defect(T1_index)
    # view(zeolite)


if __name__ == '__main__':
    main()
