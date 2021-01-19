
from maze import Zeotype

def main():
    cif_dir = "//data/GOO.cif"
    zeolite = Zeotype.build_from_cif_with_labels(cif_dir)
    z = zeolite[[i for i in range(0,30)]]
    print(id(zeolite.index_mapper), id(z.index_mapper), zeolite.index_mapper is z.index_mapper)


if __name__ == "__main__":
    main()